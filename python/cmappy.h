#ifndef CMAPPY_H
#define CMAPPY_H

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <pthread.h>
#include <zlib.h>
#include "minimap.h"
#include "mmpriv.h"
#include "kalloc.h"
#include "kseq.h"
KSEQ_DECLARE(gzFile)

typedef struct {
	const char *ctg;
	int32_t ctg_start, ctg_end;
	int32_t qry_start, qry_end;
	int32_t blen, mlen, NM, ctg_len;
	uint8_t mapq, is_primary;
	int8_t strand, trans_strand;
	int32_t seg_id;
	int32_t n_cigar32;
	uint32_t *cigar32;
} mm_hitpy_t;

typedef struct {
	int n_regs;
	int rep_len;
	mm_reg1_t *regs;
} mappy_map_result_t;

typedef struct {
	const mm_idx_t *mi;
	const mm_mapopt_t *opt;
	mm_bseq1_t *seqs;
	mappy_map_result_t *results;
	mm_tbuf_t **bufs;
	int n_seq;
	volatile int next_index;
} mappy_batch_map_t;

typedef struct {
	mappy_batch_map_t *shared;
	int tid;
} mappy_batch_worker_t;

static inline void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
{
	h->ctg = mi->seq[r->rid].name;
	h->ctg_len = mi->seq[r->rid].len;
	h->ctg_start = r->rs, h->ctg_end = r->re;
	h->qry_start = r->qs, h->qry_end = r->qe;
	h->strand = r->rev? -1 : 1;
	h->mapq = r->mapq;
	h->mlen = r->mlen;
	h->blen = r->blen;
	h->NM = r->blen - r->mlen + r->p->n_ambi;
	h->trans_strand = r->p->trans_strand == 1? 1 : r->p->trans_strand == 2? -1 : 0;
	h->is_primary = (r->id == r->parent);
	h->seg_id = r->seg_id;
	h->n_cigar32 = r->p->n_cigar;
	h->cigar32 = r->p->cigar;
}

static inline void mm_free_reg1(mm_reg1_t *r)
{
	free(r->p);
}

static inline kseq_t *mm_fastx_open(const char *fn)
{
	gzFile fp;
	fp = fn && strcmp(fn, "-") != 0? gzopen(fn, "r") : gzdopen(fileno(stdin), "r");
	return kseq_init(fp);
}

static inline void mm_fastx_close(kseq_t *ks)
{
	gzFile fp;
	fp = ks->f->f;
	kseq_destroy(ks);
	gzclose(fp);
}

static inline int mm_verbose_level(int v)
{
	if (v >= 0) mm_verbose = v;
	return mm_verbose;
}

static inline void mm_reset_timer(void)
{
	extern double realtime(void);
	mm_realtime0 = realtime();
}

static inline mm_reg1_t *mm_map_aux_nogil(const mm_idx_t *mi, const char* seqname, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt)
{
	mm_reg1_t *r;
	if (seq2 == 0) {
		r = mm_map(mi, strlen(seq1), seq1, n_regs, b, opt, seqname);
	} else {
		int _n_regs[2];
		mm_reg1_t *regs[2];
		char *seq[2];
		int i, len[2];

		len[0] = strlen(seq1);
		len[1] = strlen(seq2);
		seq[0] = (char*)seq1;
		seq[1] = strdup(seq2);
		for (i = 0; i < len[1]>>1; ++i) {
			int t = seq[1][len[1] - i - 1];
			seq[1][len[1] - i - 1] = seq_comp_table[(uint8_t)seq[1][i]];
			seq[1][i] = seq_comp_table[t];
		}
		if (len[1]&1) seq[1][len[1]>>1] = seq_comp_table[(uint8_t)seq[1][len[1]>>1]];
		mm_map_frag(mi, 2, len, (const char**)seq, _n_regs, regs, b, opt, seqname);
		for (i = 0; i < _n_regs[1]; ++i)
			regs[1][i].rev = !regs[1][i].rev;
		*n_regs = _n_regs[0] + _n_regs[1];
		regs[0] = (mm_reg1_t*)realloc(regs[0], sizeof(mm_reg1_t) * (*n_regs));
		memcpy(&regs[0][_n_regs[0]], regs[1], _n_regs[1] * sizeof(mm_reg1_t));
		free(regs[1]);
		free(seq[1]);
		r = regs[0];
	}
	return r;
}

extern unsigned char seq_comp_table[256];
static inline mm_reg1_t *mm_map_aux(const mm_idx_t *mi, const char* seqname, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt)
{
	mm_reg1_t *r;

	Py_BEGIN_ALLOW_THREADS
	r = mm_map_aux_nogil(mi, seqname, seq1, seq2, n_regs, b, opt);
	Py_END_ALLOW_THREADS

	return r;
}

static inline void mappy_bseq_free_records(mm_bseq1_t *seqs, int n_seq)
{
	int i;
	if (seqs == 0) return;
	for (i = 0; i < n_seq; ++i) {
		free(seqs[i].name);
		free(seqs[i].seq);
		free(seqs[i].qual);
		free(seqs[i].comment);
	}
	free(seqs);
}

static inline void mappy_batch_result_destroy(mappy_map_result_t *r)
{
	int i;
	if (r == 0 || r->regs == 0) return;
	for (i = 0; i < r->n_regs; ++i)
		mm_free_reg1(&r->regs[i]);
	free(r->regs);
	r->regs = 0;
	r->n_regs = 0;
	r->rep_len = 0;
}

static inline void mappy_batch_results_free(mappy_map_result_t *results, int n_seq)
{
	int i;
	if (results == 0) return;
	for (i = 0; i < n_seq; ++i)
		mappy_batch_result_destroy(&results[i]);
	free(results);
}

static void *mappy_batch_worker(void *data)
{
	mappy_batch_worker_t *w = (mappy_batch_worker_t*)data;
	mappy_batch_map_t *s = w->shared;
	for (;;) {
		int i = __sync_fetch_and_add(&s->next_index, 1);
		mappy_map_result_t *res;
		if (i >= s->n_seq) break;
		res = &s->results[i];
		res->n_regs = 0;
		res->rep_len = 0;
		res->regs = 0;
		res->regs = mm_map(s->mi, s->seqs[i].l_seq, s->seqs[i].seq, &res->n_regs, s->bufs[w->tid], s->opt, s->seqs[i].name);
		res->rep_len = s->bufs[w->tid]->rep_len;
	}
	return 0;
}

static inline int mappy_map_batch(const mm_idx_t *mi, const mm_mapopt_t *opt, mm_bseq1_t *seqs, int n_seq, int n_threads, mappy_map_result_t *results)
{
	int i;
	pthread_t *threads;
	mappy_batch_worker_t *workers;
	mappy_batch_map_t shared;
	if (n_seq <= 0) return 0;
	if (n_threads < 1) n_threads = 1;
	memset(&shared, 0, sizeof(shared));
	shared.mi = mi;
	shared.opt = opt;
	shared.seqs = seqs;
	shared.results = results;
	shared.n_seq = n_seq;
	shared.next_index = 0;
	shared.bufs = (mm_tbuf_t**)calloc(n_threads, sizeof(mm_tbuf_t*));
	for (i = 0; i < n_threads; ++i)
		shared.bufs[i] = mm_tbuf_init();
	if (n_threads == 1) {
		mappy_batch_worker_t worker;
		worker.shared = &shared;
		worker.tid = 0;
		mappy_batch_worker(&worker);
	} else {
		threads = (pthread_t*)calloc(n_threads, sizeof(pthread_t));
		workers = (mappy_batch_worker_t*)calloc(n_threads, sizeof(mappy_batch_worker_t));
		for (i = 0; i < n_threads; ++i) {
			workers[i].shared = &shared;
			workers[i].tid = i;
			pthread_create(&threads[i], 0, mappy_batch_worker, &workers[i]);
		}
		for (i = 0; i < n_threads; ++i)
			pthread_join(threads[i], 0);
		free(threads);
		free(workers);
	}
	for (i = 0; i < n_threads; ++i)
		mm_tbuf_destroy(shared.bufs[i]);
	free(shared.bufs);
	return 0;
}

static inline void mappy_kstring_reserve(kstring_t *s, size_t add)
{
	if (s->l + add + 1 > s->m) {
		s->m = s->l + add + 1;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline void mappy_kstring_clear(kstring_t *s)
{
	if (s == 0) return;
	s->l = 0;
	if (s->s) s->s[0] = 0;
}

static inline void mappy_kputsn(kstring_t *s, const char *src, size_t n)
{
	mappy_kstring_reserve(s, n);
	memcpy(s->s + s->l, src, n);
	s->l += n;
	s->s[s->l] = 0;
}

static inline void mappy_kputs(kstring_t *s, const char *src)
{
	mappy_kputsn(s, src, strlen(src));
}

static inline void mappy_kputc(kstring_t *s, char c)
{
	mappy_kstring_reserve(s, 1);
	s->s[s->l++] = c;
	s->s[s->l] = 0;
}

static inline void mappy_kputw(kstring_t *s, long long x)
{
	char buf[32];
	int i = 0;
	unsigned long long u;
	if (x < 0) {
		mappy_kputc(s, '-');
		u = (unsigned long long)(-x);
	} else u = (unsigned long long)x;
	do {
		buf[i++] = (char)('0' + (u % 10));
		u /= 10;
	} while (u > 0);
	while (i-- > 0)
		mappy_kputc(s, buf[i]);
}

static inline void mappy_kwrite_legacy_cigar(kstring_t *s, const mm_reg1_t *r)
{
	uint32_t i;
	for (i = 0; i < r->p->n_cigar; ++i) {
		mappy_kputw(s, r->p->cigar[i] >> 4);
		mappy_kputc(s, MM_CIGAR_STR[r->p->cigar[i] & 0xf]);
	}
}

static inline void mappy_write_legacy_record(kstring_t *s, const mm_idx_t *mi, const mm_bseq1_t *t, const mm_reg1_t *r, void *km, int want_cs, int want_md, char **tag_buf, int *tag_cap)
{
	mm_hitpy_t h;
	int l_tag;
	mm_reg2hitpy(mi, (mm_reg1_t*)r, &h);
	mappy_kstring_clear(s);
	mappy_kputs(s, t->name);
	mappy_kputc(s, '\t');
	mappy_kputw(s, t->l_seq);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.qry_start);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.qry_end);
	mappy_kputc(s, '\t');
	mappy_kputc(s, h.strand > 0? '+' : h.strand < 0? '-' : '?');
	mappy_kputc(s, '\t');
	mappy_kputs(s, h.ctg);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.ctg_len);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.ctg_start);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.ctg_end);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.mlen);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.blen);
	mappy_kputc(s, '\t');
	mappy_kputw(s, h.mapq);
	mappy_kputc(s, '\t');
	mappy_kputs(s, h.is_primary? "tp:A:P" : "tp:A:S");
	mappy_kputc(s, '\t');
	if (h.trans_strand > 0) mappy_kputs(s, "ts:A:+");
	else if (h.trans_strand < 0) mappy_kputs(s, "ts:A:-");
	else mappy_kputs(s, "ts:A:.");
	mappy_kputc(s, '\t');
	mappy_kputs(s, "cg:Z:");
	mappy_kwrite_legacy_cigar(s, r);
	if (want_cs) {
		l_tag = mm_gen_cs(km, tag_buf, tag_cap, mi, r, t->seq, 1);
		mappy_kputs(s, "\tcs:Z:");
		mappy_kputsn(s, *tag_buf, l_tag);
	}
	if (want_md) {
		l_tag = mm_gen_MD(km, tag_buf, tag_cap, mi, r, t->seq);
		mappy_kputs(s, "\tMD:Z:");
		mappy_kputsn(s, *tag_buf, l_tag);
	}
	mappy_kputc(s, '\n');
}

#define MAPPY_OUTPUT_LEGACY 0
#define MAPPY_OUTPUT_PAF 1
#define MAPPY_OUTPUT_SAM 2

static inline int mappy_write_sam_header(FILE *fp, const mm_idx_t *idx)
{
	uint32_t i;
	if (fp == 0) return -1;
	fprintf(fp, "@HD\tVN:1.6\tSO:unsorted\tGO:query\n");
	if (idx) {
		for (i = 0; i < idx->n_seq; ++i)
			fprintf(fp, "@SQ\tSN:%s\tLN:%u\n", idx->seq[i].name, idx->seq[i].len);
	}
	fprintf(fp, "@PG\tID:minimap2.py\tPN:minimap2.py\tVN:%s\n", MM_VERSION);
	return 0;
}

static inline int mappy_write_batch(FILE *fp, const mm_idx_t *mi, const mm_mapopt_t *opt, mm_bseq1_t *seqs, mappy_map_result_t *results, int n_seq, int output_format)
{
	int i, j, mapped = 0;
	kstring_t out = {0,0,0};
	char *tag_buf = 0;
	int tag_cap = 0;
	void *km = 0;
	int want_cs = !!(opt->flag & MM_F_OUT_CS);
	int want_md = !!(opt->flag & MM_F_OUT_MD);
	const mm_reg1_t *regss_arr[1];
	if (fp == 0) return 0;
	km = km_init();
	for (i = 0; i < n_seq; ++i) {
		mappy_map_result_t *res = &results[i];
		if (res->n_regs > 0) ++mapped;
		if (output_format == MAPPY_OUTPUT_LEGACY) {
			for (j = 0; j < res->n_regs; ++j) {
				mappy_write_legacy_record(&out, mi, &seqs[i], &res->regs[j], km, want_cs, want_md, &tag_buf, &tag_cap);
				fwrite(out.s, 1, out.l, fp);
			}
		} else if (output_format == MAPPY_OUTPUT_PAF) {
			for (j = 0; j < res->n_regs; ++j) {
				mappy_kstring_clear(&out);
				mm_write_paf4(&out, mi, &seqs[i], &res->regs[j], km, opt->flag, res->rep_len, 0, 0);
				if (out.l > 0 && out.s[out.l - 1] != '\n') mappy_kputc(&out, '\n');
				if (out.l > 0) fwrite(out.s, 1, out.l, fp);
			}
		} else if (output_format == MAPPY_OUTPUT_SAM) {
			regss_arr[0] = res->regs;
			if (res->n_regs > 0) {
				for (j = 0; j < res->n_regs; ++j) {
					mappy_kstring_clear(&out);
					mm_write_sam3(&out, mi, &seqs[i], 0, j, 1, &res->n_regs, regss_arr, km, opt->flag, res->rep_len);
					if (out.l > 0 && out.s[out.l - 1] != '\n') mappy_kputc(&out, '\n');
					if (out.l > 0) fwrite(out.s, 1, out.l, fp);
				}
			} else {
				mappy_kstring_clear(&out);
				mm_write_sam3(&out, mi, &seqs[i], 0, -1, 1, &res->n_regs, regss_arr, km, opt->flag, res->rep_len);
				if (out.l > 0 && out.s[out.l - 1] != '\n') mappy_kputc(&out, '\n');
				if (out.l > 0) fwrite(out.s, 1, out.l, fp);
			}
		}
	}
	if (fp != stdout) fflush(fp);
	free(out.s);
	free(tag_buf);
	km_destroy(km);
	return mapped;
}

static inline char *mappy_revcomp(int len, const uint8_t *seq)
{
	int i;
	char *rev;
	rev = (char*)malloc(len + 1);
	for (i = 0; i < len; ++i)
		rev[len - i - 1] = seq_comp_table[seq[i]];
	rev[len] = 0;
	return rev;
}

static char *mappy_fetch_seq(const mm_idx_t *mi, const char *name, int st, int en, int *len)
{
	int i, rid;
	char *s;
	*len = 0;
	rid = mm_idx_name2id(mi, name);
	if (rid < 0) return 0;
	if ((uint32_t)st >= mi->seq[rid].len || st >= en) return 0;
	if (en < 0 || (uint32_t)en > mi->seq[rid].len)
		en = mi->seq[rid].len;
	s = (char*)malloc(en - st + 1);
	*len = mm_idx_getseq(mi, rid, st, en, (uint8_t*)s);
	for (i = 0; i < *len; ++i)
		s[i] = "ACGTN"[(uint8_t)s[i]];
	s[*len] = 0;
	return s;
}

static mm_idx_t *mappy_idx_seq(int w, int k, int is_hpc, int bucket_bits, const char *seq, int len)
{
	const char *fake_name = "N/A";
	char *s;
	mm_idx_t *mi;
	s = (char*)calloc(len + 1, 1);
	memcpy(s, seq, len);
	mi = mm_idx_str(w, k, is_hpc, bucket_bits, 1, (const char**)&s, (const char**)&fake_name);
	free(s);
	return mi;
}

#endif
