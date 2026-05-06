#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <errno.h>
#include <math.h>
#include <stdint.h>
#include "mmpriv.h"
#include "khash.h"
#include "ksw2.h"

typedef struct mm_ra_ins_s {
	char *seq;
	uint32_t count;
	double weight;
	struct mm_ra_ins_s *next;
} mm_ra_ins_t;

KHASH_MAP_INIT_INT(ra_gap, mm_ra_ins_t*)

typedef struct {
	uint32_t base_count[5];
	double base_weight[5];
	uint32_t del_count;
	double del_weight;
	uint32_t depth;
} mm_ra_pos_t;

typedef struct {
	uint32_t ref_len;
	uint8_t *ref_seq;
	mm_ra_pos_t *pos;
	khash_t(ra_gap) *ins;

	uint32_t primary_reads_total;
	uint32_t primary_reads_used;
	uint32_t supplementary_pieces;
	uint32_t secondary_alignments;
	uint32_t forward_primary;
	uint32_t reverse_primary;

	uint64_t sum_read_len;
	uint64_t sum_mapq;
	double sum_identity;
	uint64_t sum_softclip5;
	uint64_t sum_softclip3;
	uint64_t qual_sum;
	uint64_t qual_bases;

	uint64_t mismatch_bases;
	uint64_t insertion_bases;
	uint64_t deletion_bases;
	uint64_t skipped_bases;

	uint32_t *read_lens;
	uint32_t n_read_lens, m_read_lens;
	uint32_t *mapqs;
	uint32_t n_mapqs, m_mapqs;
	double *identities;
	uint32_t n_identities, m_identities;
} mm_ra_ref_t;

struct mm_ref_analysis_collect_s {
	const mm_idx_t *mi;
	mm_ref_analysis_opt_t opt;
	mm_ra_ref_t *refs;
	int64_t total_reads;
	int64_t processed_reads;
	int64_t report_every;
	int64_t next_report;
	double start_time;
};

static void ra_format_duration(double seconds, char *buf, size_t buf_len)
{
	if (seconds < 0.0) seconds = 0.0;
	if (seconds < 60.0) {
		snprintf(buf, buf_len, "%.1fs", seconds);
	} else {
		int total = (int)(seconds + 0.5);
		int sec = total % 60;
		int min = (total / 60) % 60;
		int hrs = total / 3600;
		if (hrs > 0) snprintf(buf, buf_len, "%dh%02dm%02ds", hrs, min, sec);
		else snprintf(buf, buf_len, "%dm%02ds", min, sec);
	}
}

static void ra_emit_progress(const mm_ref_analysis_collect_t *c, int final)
{
	double elapsed, throughput;
	char elapsed_buf[32], eta_buf[32];
	if (c == 0) return;
	elapsed = realtime() - c->start_time;
	if (elapsed < 0.0) elapsed = 0.0;
	throughput = elapsed > 0.0? c->processed_reads / elapsed : 0.0;
	ra_format_duration(elapsed, elapsed_buf, sizeof(elapsed_buf));
	if (c->total_reads > 0 && throughput > 0.0) {
		double remaining = (c->total_reads - c->processed_reads) / throughput;
		if (remaining < 0.0) remaining = 0.0;
		ra_format_duration(remaining, eta_buf, sizeof(eta_buf));
		fprintf(stderr, "[M::ref_analysis] %s %lld/%lld reads; elapsed %s; throughput %.2f reads/s; ETA %s\n",
			final? "completed after analyzing" : "analyzed",
			(long long)c->processed_reads, (long long)c->total_reads, elapsed_buf, throughput, eta_buf);
	} else {
		strcpy(eta_buf, "unknown");
		fprintf(stderr, "[M::ref_analysis] %s %lld reads; elapsed %s; throughput %.2f reads/s; ETA %s\n",
			final? "completed after analyzing" : "analyzed",
			(long long)c->processed_reads, elapsed_buf, throughput, eta_buf);
	}
	fflush(stderr);
}

static void ra_free_bseq_records(mm_bseq1_t *seqs, int n)
{
	int i;
	if (seqs == 0) return;
	for (i = 0; i < n; ++i) {
		free(seqs[i].name);
		free(seqs[i].seq);
		if (seqs[i].qual) free(seqs[i].qual);
		if (seqs[i].comment) free(seqs[i].comment);
	}
	free(seqs);
}

int64_t mm_ref_analysis_count_reads(int n_segs, const char **fn)
{
	int i, n = 0;
	int64_t total = 0;
	mm_bseq1_t *seqs = 0;
	if (n_segs < 1 || fn == 0) return -1;
	for (i = 0; i < n_segs; ++i)
		if (fn[i] == 0 || strcmp(fn[i], "-") == 0)
			return -1;
	if (n_segs == 1) {
		mm_bseq_file_t *fp = mm_bseq_open(fn[0]);
		if (fp == 0) return -1;
		while ((seqs = mm_bseq_read3(fp, 1<<20, 0, 0, 0, &n)) != 0 && n > 0) {
			total += n;
			ra_free_bseq_records(seqs, n);
			seqs = 0;
		}
		mm_bseq_close(fp);
	} else {
		mm_bseq_file_t **fps = CALLOC(mm_bseq_file_t*, n_segs);
		if (fps == 0) return -1;
		for (i = 0; i < n_segs; ++i) {
			fps[i] = mm_bseq_open(fn[i]);
			if (fps[i] == 0) break;
		}
		if (i < n_segs) {
			while (--i >= 0) mm_bseq_close(fps[i]);
			free(fps);
			return -1;
		}
		while ((seqs = mm_bseq_read_frag2(n_segs, fps, 1<<20, 0, 0, &n)) != 0 && n > 0) {
			total += n / n_segs;
			ra_free_bseq_records(seqs, n);
			seqs = 0;
		}
		for (i = 0; i < n_segs; ++i) mm_bseq_close(fps[i]);
		free(fps);
	}
	return total;
}

static void ra_str_enlarge(kstring_t *s, int l)
{
	if (s->l + l + 1 > s->m) {
		s->m = s->l + l + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static void ra_str_putc(kstring_t *s, char c)
{
	ra_str_enlarge(s, 1);
	s->s[s->l++] = c;
	s->s[s->l] = 0;
}

static void ra_str_puts(kstring_t *s, const char *src)
{
	int l = strlen(src);
	ra_str_enlarge(s, l);
	memcpy(s->s + s->l, src, l);
	s->l += l;
	s->s[s->l] = 0;
}

static void ra_str_putw(kstring_t *s, uint32_t x)
{
	char buf[32];
	snprintf(buf, sizeof(buf), "%u", x);
	ra_str_puts(s, buf);
}

static int cmp_u32(const void *a, const void *b)
{
	uint32_t x = *(const uint32_t*)a, y = *(const uint32_t*)b;
	return x < y? -1 : x > y;
}

static int cmp_dbl(const void *a, const void *b)
{
	double x = *(const double*)a, y = *(const double*)b;
	return x < y? -1 : x > y;
}

static double median_u32(const uint32_t *a, uint32_t n)
{
	uint32_t *tmp;
	double med;
	if (n == 0) return -1.0;
	tmp = MALLOC(uint32_t, n);
	memcpy(tmp, a, n * sizeof(uint32_t));
	qsort(tmp, n, sizeof(uint32_t), cmp_u32);
	med = (n & 1)? tmp[n>>1] : (tmp[(n>>1) - 1] + tmp[n>>1]) * 0.5;
	free(tmp);
	return med;
}

static double median_dbl(const double *a, uint32_t n)
{
	double *tmp, med;
	if (n == 0) return -1.0;
	tmp = MALLOC(double, n);
	memcpy(tmp, a, n * sizeof(double));
	qsort(tmp, n, sizeof(double), cmp_dbl);
	med = (n & 1)? tmp[n>>1] : (tmp[(n>>1) - 1] + tmp[n>>1]) * 0.5;
	free(tmp);
	return med;
}

static void push_u32(uint32_t **a, uint32_t *n, uint32_t *m, uint32_t v)
{
	if (*n == *m) {
		*m = *m? *m << 1 : 8;
		*a = REALLOC(uint32_t, *a, *m);
	}
	(*a)[(*n)++] = v;
}

static void push_dbl(double **a, uint32_t *n, uint32_t *m, double v)
{
	if (*n == *m) {
		*m = *m? *m << 1 : 8;
		*a = REALLOC(double, *a, *m);
	}
	(*a)[(*n)++] = v;
}

static double read_mean_qual(const mm_bseq1_t *t)
{
	int i;
	double s = 0.0;
	if (t->qual == 0 || t->l_seq == 0) return -1.0;
	for (i = 0; i < t->l_seq; ++i) {
		int q = (uint8_t)t->qual[i] - 33;
		s += q > 0? q : 0;
	}
	return s / t->l_seq;
}

static void append_cigar_text(kstring_t *s, uint32_t op, uint32_t len)
{
	if (len == 0) return;
	ra_str_putw(s, len);
	ra_str_putc(s, "MIDNSHP=XB"[op]);
}

static void free_ins_list(mm_ra_ins_t *p)
{
	while (p) {
		mm_ra_ins_t *q = p->next;
		free(p->seq);
		free(p);
		p = q;
	}
}

static void mm_ra_ref_destroy(mm_ra_ref_t *r)
{
	khiter_t k;
	if (r->ins) {
		for (k = 0; k < kh_end(r->ins); ++k)
			if (kh_exist(r->ins, k))
				free_ins_list(kh_val(r->ins, k));
		kh_destroy(ra_gap, r->ins);
	}
	free(r->ref_seq);
	free(r->pos);
	free(r->read_lens);
	free(r->mapqs);
	free(r->identities);
	memset(r, 0, sizeof(*r));
}

static int mm_ra_ref_init(mm_ra_ref_t *r, const mm_idx_t *mi, uint32_t rid)
{
	if (r->ref_seq) return 0;
	r->ref_len = mi->seq[rid].len;
	r->ref_seq = MALLOC(uint8_t, r->ref_len);
	r->pos = CALLOC(mm_ra_pos_t, r->ref_len);
	r->ins = kh_init(ra_gap);
	if (r->ref_seq == 0 || r->pos == 0 || r->ins == 0) return -1;
	mm_idx_getseq(mi, rid, 0, r->ref_len, r->ref_seq);
	return 0;
}

static int mm_ra_is_primary(const mm_reg1_t *r)
{
	return r && r->p && r->parent == r->id && r->sam_pri;
}

static int mm_ra_is_supp(const mm_reg1_t *r)
{
	return r && r->p && r->parent == r->id && !r->sam_pri;
}

static int mm_ra_is_secondary(const mm_reg1_t *r)
{
	return r && r->p && r->parent != r->id;
}

static double ins_weight(const uint8_t *qqual, int qoff, int len, int has_qual)
{
	int i;
	double w = 0.0;
	if (!has_qual || len <= 0) return 1.0;
	for (i = 0; i < len; ++i) w += qqual[qoff + i] + 1.0;
	return w / len;
}

static void add_insertion(mm_ra_ref_t *rr, uint32_t gap, const uint8_t *qseq, const uint8_t *qqual, int qoff, int len, int has_qual)
{
	int ret;
	khiter_t k;
	mm_ra_ins_t *p;
	char *seq = MALLOC(char, len + 1);
	int i;
	for (i = 0; i < len; ++i)
		seq[i] = "ACGTN"[qseq[qoff + i] < 5? qseq[qoff + i] : 4];
	seq[len] = 0;
	k = kh_put(ra_gap, rr->ins, gap, &ret);
	if (ret) kh_val(rr->ins, k) = 0;
	for (p = kh_val(rr->ins, k); p; p = p->next)
		if (strcmp(p->seq, seq) == 0)
			break;
	if (p == 0) {
		p = CALLOC(mm_ra_ins_t, 1);
		p->seq = seq;
		p->next = kh_val(rr->ins, k);
		kh_val(rr->ins, k) = p;
	} else free(seq);
	++p->count;
	p->weight += ins_weight(qqual, qoff, len, has_qual);
	rr->insertion_bases += len;
}

static void extract_oriented_query(const mm_bseq1_t *t, const mm_reg1_t *r, uint8_t **qseq_out, uint8_t **qqual_out, int *has_qual_out)
{
	int i, len = r->qe - r->qs;
	uint8_t *qseq = MALLOC(uint8_t, len);
	uint8_t *qqual = 0;
	int has_qual = t->qual != 0;
	if (has_qual) qqual = MALLOC(uint8_t, len);
	if (!r->rev) {
		for (i = r->qs; i < r->qe; ++i) {
			qseq[i - r->qs] = seq_nt4_table[(uint8_t)t->seq[i]];
			if (has_qual) {
				int q = (uint8_t)t->qual[i] - 33;
				qqual[i - r->qs] = q > 0? q : 0;
			}
		}
	} else {
		for (i = r->qs; i < r->qe; ++i) {
			uint8_t c = seq_nt4_table[(uint8_t)t->seq[i]];
			int j = r->qe - i - 1;
			qseq[j] = c >= 4? 4 : 3 - c;
			if (has_qual) {
				int q = (uint8_t)t->qual[i] - 33;
				qqual[j] = q > 0? q : 0;
			}
		}
	}
	*qseq_out = qseq;
	*qqual_out = qqual;
	*has_qual_out = has_qual;
}

static void accumulate_primary(mm_ref_analysis_collect_t *c, const mm_bseq1_t *t, const mm_reg1_t *r)
{
	mm_ra_ref_t *rr = &c->refs[r->rid];
	uint8_t *qseq, *qqual;
	int has_qual, qoff = 0, roff = r->rs, i, soft5, soft3;
	double identity, meanq;
	if (mm_ra_ref_init(rr, c->mi, r->rid) < 0) return;
	++rr->primary_reads_total;
	if (r->mapq < c->opt.min_mapq || r->qe - r->qs < c->opt.min_aln_len) return;
	++rr->primary_reads_used;
	if (r->rev) ++rr->reverse_primary;
	else ++rr->forward_primary;
	push_u32(&rr->read_lens, &rr->n_read_lens, &rr->m_read_lens, t->l_seq);
	push_u32(&rr->mapqs, &rr->n_mapqs, &rr->m_mapqs, r->mapq);
	identity = r->blen > 0? (double)r->mlen / r->blen : 0.0;
	push_dbl(&rr->identities, &rr->n_identities, &rr->m_identities, identity);
	rr->sum_read_len += t->l_seq;
	rr->sum_mapq += r->mapq;
	rr->sum_identity += identity;
	soft5 = r->rev? t->l_seq - r->qe : r->qs;
	soft3 = r->rev? r->qs : t->l_seq - r->qe;
	rr->sum_softclip5 += soft5;
	rr->sum_softclip3 += soft3;
	meanq = read_mean_qual(t);
	if (c->opt.use_qual && meanq >= 0.0) {
		int j;
		for (j = 0; j < t->l_seq; ++j) {
			int q = (uint8_t)t->qual[j] - 33;
			rr->qual_sum += q > 0? q : 0;
		}
		rr->qual_bases += t->l_seq;
	}
	extract_oriented_query(t, r, &qseq, &qqual, &has_qual);
	for (i = 0; i < (int)r->p->n_cigar; ++i) {
		uint32_t op = r->p->cigar[i] & 0xf, len = r->p->cigar[i] >> 4, j;
		if (op == MM_CIGAR_MATCH || op == MM_CIGAR_EQ_MATCH || op == MM_CIGAR_X_MISMATCH) {
			for (j = 0; j < len; ++j) {
				uint8_t qb = qseq[qoff + j] < 5? qseq[qoff + j] : 4;
				double qw = (c->opt.use_qual && has_qual)? qqual[qoff + j] + 1.0 : 1.0;
				mm_ra_pos_t *p = &rr->pos[roff + j];
				++p->base_count[qb];
				p->base_weight[qb] += qw;
				++p->depth;
				if (qb != rr->ref_seq[roff + j]) ++rr->mismatch_bases;
			}
			qoff += len;
			roff += len;
		} else if (op == MM_CIGAR_INS) {
			add_insertion(rr, roff, qseq, qqual, qoff, len, c->opt.use_qual && has_qual);
			qoff += len;
		} else if (op == MM_CIGAR_DEL || op == MM_CIGAR_N_SKIP) {
			for (j = 0; j < len; ++j) {
				mm_ra_pos_t *p = &rr->pos[roff + j];
				++p->del_count;
				p->del_weight += 1.0;
				++p->depth;
			}
			if (op == MM_CIGAR_DEL) rr->deletion_bases += len;
			else rr->skipped_bases += len;
			roff += len;
		}
	}
	free(qseq);
	free(qqual);
}

mm_ref_analysis_collect_t *mm_ref_analysis_collect_init(const mm_idx_t *mi, const mm_ref_analysis_opt_t *opt)
{
	mm_ref_analysis_collect_t *c = CALLOC(mm_ref_analysis_collect_t, 1);
	if (c == 0) return 0;
	c->mi = mi;
	if (opt) c->opt = *opt;
	else mm_ref_analysis_opt_init(&c->opt);
	c->refs = CALLOC(mm_ra_ref_t, mi->n_seq);
	if (c->refs == 0) {
		free(c);
		return 0;
	}
	c->report_every = 1000;
	c->next_report = 1000;
	c->start_time = realtime();
	return c;
}

void mm_ref_analysis_collect_destroy(mm_ref_analysis_collect_t *c)
{
	uint32_t i;
	if (c == 0) return;
	if (c->refs) {
		for (i = 0; i < c->mi->n_seq; ++i)
			mm_ra_ref_destroy(&c->refs[i]);
		free(c->refs);
	}
	free(c);
}

void mm_ref_analysis_collect_set_total(mm_ref_analysis_collect_t *c, int64_t total_reads)
{
	if (c == 0) return;
	c->total_reads = total_reads;
	c->report_every = total_reads > 0 && total_reads < 10000? total_reads / 10 : 1000;
	if (c->report_every < 1) c->report_every = 1;
	c->next_report = c->report_every;
	c->start_time = realtime();
	if (total_reads > 0)
		fprintf(stderr, "[M::ref_analysis] starting per-reference analysis for %lld reads; update interval %lld\n",
			(long long)total_reads, (long long)c->report_every);
	else
		fprintf(stderr, "[M::ref_analysis] starting per-reference analysis; total reads unknown\n");
	fflush(stderr);
}

void mm_ref_analysis_step(mm_ref_analysis_collect_t *c, const mm_bseq1_t *t, int n_reg, mm_reg1_t *regs)
{
	int i;
	const mm_reg1_t *primary = 0;
	if (c == 0 || t == 0) return;
	if (regs && n_reg > 0) {
		for (i = 0; i < n_reg; ++i) {
			mm_reg1_t *r = &regs[i];
			if (r->p == 0) continue;
			if (mm_ra_is_secondary(r)) {
				mm_ra_ref_t *rr = &c->refs[r->rid];
				if (mm_ra_ref_init(rr, c->mi, r->rid) == 0) ++rr->secondary_alignments;
			} else if (mm_ra_is_supp(r)) {
				mm_ra_ref_t *rr = &c->refs[r->rid];
				if (mm_ra_ref_init(rr, c->mi, r->rid) == 0) ++rr->supplementary_pieces;
			} else if (mm_ra_is_primary(r) && primary == 0) {
				primary = r;
			}
		}
		if (primary) accumulate_primary(c, t, primary);
	}
	++c->processed_reads;
	if (c->processed_reads >= c->next_report || (c->total_reads > 0 && c->processed_reads >= c->total_reads)) {
		ra_emit_progress(c, 0);
		c->next_report = c->processed_reads + c->report_every;
	}
}

static double position_support(const mm_ra_pos_t *p, int base, int use_weight)
{
	return use_weight? p->base_weight[base] : p->base_count[base];
}

static double deletion_support(const mm_ra_pos_t *p, int use_weight)
{
	return use_weight? p->del_weight : p->del_count;
}

static double insertion_score(const mm_ra_ins_t *p, int use_weight)
{
	return use_weight? p->weight : p->count;
}

static double gap_support(const mm_ra_ref_t *rr, uint32_t gap)
{
	uint32_t left = gap > 0? rr->pos[gap - 1].depth : 0;
	uint32_t right = gap < rr->ref_len? rr->pos[gap].depth : 0;
	return left > right? left : right;
}

static const mm_ra_ins_t *best_insertion_at_gap(const mm_ra_ref_t *rr, uint32_t gap, int use_weight, double *best_score_out, double *span_out)
{
	khiter_t k;
	const mm_ra_ins_t *best = 0;
	mm_ra_ins_t *p;
	double best_score = 0.0;
	double span = gap_support(rr, gap);
	if (best_score_out) *best_score_out = 0.0;
	if (span_out) *span_out = span;
	if (rr->ins == 0 || (k = kh_get(ra_gap, rr->ins, gap)) == kh_end(rr->ins))
		return 0;
	for (p = kh_val(rr->ins, k); p; p = p->next) {
		double sc = insertion_score(p, use_weight);
		if (best == 0 || sc > best_score || (sc == best_score && strcmp(p->seq, best->seq) < 0))
			best = p, best_score = sc;
	}
	if (best_score_out) *best_score_out = best_score;
	if (best && best->count > 0 && best_score * 2.0 > span)
		return best;
	return 0;
}

static void append_cigar_run(kstring_t *s, uint32_t *last_op, uint32_t *last_len, uint32_t op, uint32_t len)
{
	if (len == 0) return;
	if (*last_len > 0 && *last_op == op) {
		*last_len += len;
	} else {
		if (*last_len > 0) append_cigar_text(s, *last_op, *last_len);
		*last_op = op;
		*last_len = len;
	}
}

static void build_consensus_and_cigar(const mm_ra_ref_t *rr, int use_weight, char **cons_out, uint32_t *len_out, char **cigar_out, uint64_t *eq_out, uint64_t *x_out, uint64_t *i_out, uint64_t *d_out, uint64_t *n_out, double *edit_out)
{
	kstring_t cons = {0,0,0}, cigar = {0,0,0};
	uint32_t i, last_op = UINT32_MAX, last_len = 0;
	uint64_t eq = 0, x = 0, ins = 0, del = 0;
	for (i = 0; i <= rr->ref_len; ++i) {
		const mm_ra_ins_t *best_ins;
		double best_score = 0.0, span = 0.0;
		best_ins = best_insertion_at_gap(rr, i, use_weight, &best_score, &span);
		if (best_ins) {
			uint32_t ins_len = (uint32_t)strlen(best_ins->seq);
			if (ins_len > 0) {
				ra_str_puts(&cons, best_ins->seq);
				append_cigar_run(&cigar, &last_op, &last_len, MM_CIGAR_INS, ins_len);
				ins += ins_len;
			}
		}
		if (i == rr->ref_len) break;
		if (rr->pos[i].depth == 0) {
			ra_str_putc(&cons, "ACGTN"[rr->ref_seq[i] < 5? rr->ref_seq[i] : 4]);
			append_cigar_run(&cigar, &last_op, &last_len, MM_CIGAR_EQ_MATCH, 1);
			++eq;
		} else {
			int b, refb = rr->ref_seq[i] < 5? rr->ref_seq[i] : 4, bestb = refb;
			double best = position_support(&rr->pos[i], refb, use_weight);
			double del_support = deletion_support(&rr->pos[i], use_weight);
			for (b = 0; b < 5; ++b) {
				double sc = position_support(&rr->pos[i], b, use_weight);
				if (sc > best) best = sc, bestb = b;
			}
			if (del_support > best) {
				append_cigar_run(&cigar, &last_op, &last_len, MM_CIGAR_DEL, 1);
				++del;
			} else {
				ra_str_putc(&cons, "ACGTN"[bestb]);
				if (bestb == refb) {
					append_cigar_run(&cigar, &last_op, &last_len, MM_CIGAR_EQ_MATCH, 1);
					++eq;
				} else {
					append_cigar_run(&cigar, &last_op, &last_len, MM_CIGAR_X_MISMATCH, 1);
					++x;
				}
			}
		}
	}
	if (last_len > 0) append_cigar_text(&cigar, last_op, last_len);
	*cons_out = cons.s;
	*len_out = cons.l;
	*cigar_out = cigar.s;
	*eq_out = eq;
	*x_out = x;
	*i_out = ins;
	*d_out = del;
	*n_out = 0;
	*edit_out = cons.l > 0? (double)(x + ins + del) / cons.l : 0.0;
}

static void fill_row(const mm_idx_t *mi, const mm_ra_ref_t *rr, uint32_t rid, mm_ref_analysis_row_t *row)
{
	uint32_t i;
	uint64_t depth_sum = 0;
	uint32_t *depths;
	double use_weight = rr->qual_bases > 0;
	memset(row, 0, sizeof(*row));
	row->ref_name = strdup(mi->seq[rid].name);
	row->ref_len = mi->seq[rid].len;
	row->primary_reads_total = rr->primary_reads_total;
	row->primary_reads_used = rr->primary_reads_used;
	row->supplementary_pieces = rr->supplementary_pieces;
	row->secondary_alignments = rr->secondary_alignments;
	row->forward_primary = rr->forward_primary;
	row->reverse_primary = rr->reverse_primary;
	row->median_read_len = row->mean_read_len = row->avg_read_qual = -1.0;
	row->mean_mapq = row->median_mapq = -1.0;
	row->mean_identity = row->median_identity = -1.0;
	row->coverage_breadth = row->mean_depth = row->median_depth = -1.0;
	row->mean_softclip_5p = row->mean_softclip_3p = -1.0;
	row->consensus_edit_rate = -1.0;
	row->mismatch_bases = rr->mismatch_bases;
	row->insertion_bases = rr->insertion_bases;
	row->deletion_bases = rr->deletion_bases;
	row->skipped_bases = rr->skipped_bases;
	if (rr->primary_reads_used == 0) return;
	row->min_read_len = UINT32_MAX;
	for (i = 0; i < rr->n_read_lens; ++i) {
		if (rr->read_lens[i] < row->min_read_len) row->min_read_len = rr->read_lens[i];
		if (rr->read_lens[i] > row->max_read_len) row->max_read_len = rr->read_lens[i];
	}
	row->median_read_len = median_u32(rr->read_lens, rr->n_read_lens);
	row->mean_read_len = (double)rr->sum_read_len / rr->primary_reads_used;
	if (rr->qual_bases > 0) row->avg_read_qual = (double)rr->qual_sum / rr->qual_bases;
	row->mean_mapq = (double)rr->sum_mapq / rr->primary_reads_used;
	row->median_mapq = median_u32(rr->mapqs, rr->n_mapqs);
	row->mean_identity = rr->sum_identity / rr->primary_reads_used;
	row->median_identity = median_dbl(rr->identities, rr->n_identities);
	row->mean_softclip_5p = (double)rr->sum_softclip5 / rr->primary_reads_used;
	row->mean_softclip_3p = (double)rr->sum_softclip3 / rr->primary_reads_used;
	depths = MALLOC(uint32_t, rr->ref_len);
	for (i = 0; i < rr->ref_len; ++i) {
		depths[i] = rr->pos[i].depth;
		if (depths[i] > 0) ++row->covered_bases;
		depth_sum += depths[i];
	}
	row->coverage_breadth = rr->ref_len > 0? (double)row->covered_bases / rr->ref_len : 0.0;
	row->mean_depth = rr->ref_len > 0? (double)depth_sum / rr->ref_len : 0.0;
	row->median_depth = median_u32(depths, rr->ref_len);
	free(depths);
	build_consensus_and_cigar(rr, use_weight > 0.0, &row->consensus_seq, &row->consensus_len, &row->consensus_cigar, &row->cigar_eq, &row->cigar_x, &row->cigar_i, &row->cigar_d, &row->cigar_n, &row->consensus_edit_rate);
}

int mm_ref_analysis_collect_finish(mm_ref_analysis_collect_t *c, mm_ref_analysis_result_t **out)
{
	uint32_t i;
	mm_ref_analysis_result_t *res;
	if (out == 0 || c == 0) return -1;
	res = CALLOC(mm_ref_analysis_result_t, 1);
	if (res == 0) return -1;
	res->n_rows = c->mi->n_seq;
	res->rows = CALLOC(mm_ref_analysis_row_t, res->n_rows);
	if (res->rows == 0) {
		free(res);
		return -1;
	}
	for (i = 0; i < res->n_rows; ++i)
		fill_row(c->mi, &c->refs[i], i, &res->rows[i]);
	ra_emit_progress(c, 1);
	*out = res;
	return 0;
}

void mm_ref_analysis_result_destroy(mm_ref_analysis_result_t *res)
{
	uint32_t i;
	if (res == 0) return;
	for (i = 0; i < res->n_rows; ++i) {
		free(res->rows[i].ref_name);
		free(res->rows[i].consensus_seq);
		free(res->rows[i].consensus_cigar);
	}
	free(res->rows);
	free(res);
}

static void fputs_double(FILE *fp, double x)
{
	if (x < 0.0) fputs("NA", fp);
	else fprintf(fp, "%.6f", x);
}

int mm_ref_analysis_write(const mm_ref_analysis_result_t *res, const char *prefix)
{
	char *fa_fn, *stats_fn, *cigar_fn;
	FILE *fa, *stats, *cigar;
	uint32_t i;
	size_t l;
	if (res == 0 || prefix == 0) return -1;
	l = strlen(prefix);
	fa_fn = MALLOC(char, l + 32);
	stats_fn = MALLOC(char, l + 40);
	cigar_fn = MALLOC(char, l + 40);
	sprintf(fa_fn, "%s.consensus.fa", prefix);
	sprintf(stats_fn, "%s.reference_stats.tsv", prefix);
	sprintf(cigar_fn, "%s.consensus_cigar.tsv", prefix);
	fa = fopen(fa_fn, "w");
	stats = fopen(stats_fn, "w");
	cigar = fopen(cigar_fn, "w");
	if (!fa || !stats || !cigar) {
		if (fa) fclose(fa);
		if (stats) fclose(stats);
		if (cigar) fclose(cigar);
		free(fa_fn); free(stats_fn); free(cigar_fn);
		return -1;
	}
	fprintf(stats, "ref_name\tref_len\tprimary_reads_total\tprimary_reads_used\tsupplementary_pieces\tsecondary_alignments\tmin_read_len\tmax_read_len\tmedian_read_len\tmean_read_len\tavg_read_qual\tmean_mapq\tmedian_mapq\tmean_identity\tmedian_identity\tforward_primary\treverse_primary\tcovered_bases\tcoverage_breadth\tmean_depth\tmedian_depth\tmean_softclip_5p\tmean_softclip_3p\tmismatch_bases\tinsertion_bases\tdeletion_bases\tskipped_bases\tconsensus_len\tconsensus_edit_rate\n");
	fprintf(cigar, "ref_name\tref_len\tconsensus_len\tcigar\teq_count\tx_count\ti_count\td_count\tn_count\tconsensus_edit_rate\n");
	for (i = 0; i < res->n_rows; ++i) {
		const mm_ref_analysis_row_t *r = &res->rows[i];
		if (r->primary_reads_total > 0 && r->consensus_seq)
			fprintf(fa, ">%s\n%s\n", r->ref_name, r->consensus_seq);
		fprintf(stats, "%s\t%u\t%u\t%u\t%u\t%u\t", r->ref_name, r->ref_len, r->primary_reads_total, r->primary_reads_used, r->supplementary_pieces, r->secondary_alignments);
		if (r->primary_reads_used == 0) {
			fprintf(stats, "0\t0\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t0\t0\t0\tNA\tNA\tNA\t0\t0\t0\t0\t0\t0\tNA\n");
		} else {
			fprintf(stats, "%u\t%u\t", r->min_read_len, r->max_read_len);
			fputs_double(stats, r->median_read_len); fputc('\t', stats);
			fputs_double(stats, r->mean_read_len); fputc('\t', stats);
			fputs_double(stats, r->avg_read_qual); fputc('\t', stats);
			fputs_double(stats, r->mean_mapq); fputc('\t', stats);
			fputs_double(stats, r->median_mapq); fputc('\t', stats);
			fputs_double(stats, r->mean_identity); fputc('\t', stats);
			fputs_double(stats, r->median_identity); fputc('\t', stats);
			fprintf(stats, "%u\t%u\t%u\t", r->forward_primary, r->reverse_primary, r->covered_bases);
			fputs_double(stats, r->coverage_breadth); fputc('\t', stats);
			fputs_double(stats, r->mean_depth); fputc('\t', stats);
			fputs_double(stats, r->median_depth); fputc('\t', stats);
			fputs_double(stats, r->mean_softclip_5p); fputc('\t', stats);
			fputs_double(stats, r->mean_softclip_3p); fputc('\t', stats);
			fprintf(stats, "%llu\t%llu\t%llu\t%llu\t%u\t",
				(unsigned long long)r->mismatch_bases,
				(unsigned long long)r->insertion_bases,
				(unsigned long long)r->deletion_bases,
				(unsigned long long)r->skipped_bases,
				r->consensus_len);
			fputs_double(stats, r->consensus_edit_rate);
			fputc('\n', stats);
		}
		fprintf(cigar, "%s\t%u\t%u\t%s\t%llu\t%llu\t%llu\t%llu\t%llu\t",
			r->ref_name,
			r->ref_len,
			r->consensus_len,
			r->consensus_cigar? r->consensus_cigar : "*",
			(unsigned long long)r->cigar_eq,
			(unsigned long long)r->cigar_x,
			(unsigned long long)r->cigar_i,
			(unsigned long long)r->cigar_d,
			(unsigned long long)r->cigar_n);
		fputs_double(cigar, r->consensus_edit_rate);
		fputc('\n', cigar);
	}
	fclose(fa);
	fclose(stats);
	fclose(cigar);
	free(fa_fn); free(stats_fn); free(cigar_fn);
	return 0;
}

int mm_ref_analysis_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, const mm_ref_analysis_opt_t *aopt, mm_ref_analysis_result_t **out)
{
	mm_ref_analysis_collect_t *c;
	mm_mapopt_t opt2;
	mm_ref_analysis_opt_t local_opt;
	int ret;
	if (out == 0) return -1;
	*out = 0;
	if (opt && opt->split_prefix) return -2;
	if (aopt == 0) {
		mm_ref_analysis_opt_init(&local_opt);
		aopt = &local_opt;
	}
	opt2 = *opt;
	opt2.flag |= MM_F_CIGAR | MM_F_REF_ANALYSIS | MM_F_NO_PRINT;
	c = mm_ref_analysis_collect_init(idx, aopt);
	if (c == 0) return -1;
	if (aopt->count_reads_for_eta)
		mm_ref_analysis_collect_set_total(c, mm_ref_analysis_count_reads(1, &fn));
	else
		mm_ref_analysis_collect_set_total(c, 0);
	ret = mm_map_file_frag_collect(idx, 1, &fn, &opt2, n_threads, aopt, c);
	if (ret == 0)
		ret = mm_ref_analysis_collect_finish(c, out);
	mm_ref_analysis_collect_destroy(c);
	return ret;
}
