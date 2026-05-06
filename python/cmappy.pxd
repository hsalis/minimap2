from libc.stdint cimport int8_t, uint8_t, int32_t, int64_t, uint32_t, uint64_t
from libc.stdio cimport FILE

cdef extern from "minimap.h":
	#
	# Options
	#
	ctypedef struct mm_idxopt_t:
		short k, w, flag, bucket_bits
		int64_t mini_batch_size
		uint64_t batch_size
		int32_t compact_k
		float compact_ratio

	ctypedef struct mm_mapopt_t:
		int64_t flag
		int seed
		int sdust_thres

		int max_qlen

		int bw, bw_long
		int max_gap, max_gap_ref
		int max_frag_len
		int max_chain_skip, max_chain_iter
		int min_cnt
		int min_chain_score
		float chain_gap_scale
		float chain_skip_scale
		int rmq_size_cap, rmq_inner_dist
		int rmq_rescue_size
		float rmq_rescue_ratio

		float mask_level
		int mask_len
		float pri_ratio
		int best_n

		float alt_drop

		int a, b, q, e, q2, e2
		int transition
		int sc_ambi
		int noncan
		int junc_bonus, junc_pen
		int zdrop, zdrop_inv
		int end_bonus
		int min_dp_max
		int min_ksw_len
		int anchor_ext_len, anchor_ext_shift
		float max_clip_ratio

		int rank_min_len
		float rank_frac

		int pe_ori, pe_bonus

		int jump_min_match

		float mid_occ_frac
		float q_occ_frac
		int32_t min_mid_occ
		int32_t mid_occ
		int32_t max_occ
		int64_t mini_batch_size
		int64_t max_sw_mat
		int64_t cap_kalloc

		const char *split_prefix

	ctypedef struct mm_ref_analysis_opt_t:
		int32_t min_mapq
		int32_t min_aln_len
		int32_t use_qual
		int32_t count_reads_for_eta

	ctypedef struct mm_ref_analysis_row_t:
		char *ref_name
		uint32_t ref_len
		uint32_t primary_reads_total
		uint32_t primary_reads_used
		uint32_t supplementary_pieces
		uint32_t secondary_alignments
		uint32_t min_read_len
		uint32_t max_read_len
		double median_read_len
		double mean_read_len
		double avg_read_qual
		double mean_mapq
		double median_mapq
		double mean_identity
		double median_identity
		uint32_t forward_primary
		uint32_t reverse_primary
		uint32_t covered_bases
		double coverage_breadth
		double mean_depth
		double median_depth
		double mean_softclip_5p
		double mean_softclip_3p
		uint64_t mismatch_bases
		uint64_t insertion_bases
		uint64_t deletion_bases
		uint64_t skipped_bases
		char *consensus_seq
		uint32_t consensus_len
		char *consensus_cigar
		uint64_t cigar_eq
		uint64_t cigar_x
		uint64_t cigar_i
		uint64_t cigar_d
		uint64_t cigar_n
		double consensus_edit_rate

	ctypedef struct mm_ref_analysis_result_t:
		uint32_t n_rows
		mm_ref_analysis_row_t *rows

	int mm_set_opt(char *preset, mm_idxopt_t *io, mm_mapopt_t *mo)
	int mm_verbose
	void mm_ref_analysis_opt_init(mm_ref_analysis_opt_t *opt)

	#
	# Indexing
	#
	ctypedef struct mm_idx_seq_t:
		char *name
		uint64_t offset
		uint32_t len

	ctypedef struct mm_idx_bucket_t:
		pass

	ctypedef struct mm_idx_t:
		int32_t b, w, k, flag
		uint32_t n_seq
		mm_idx_seq_t *seq
		uint32_t *S
		mm_idx_bucket_t *B
		void *km
		void *h

	ctypedef struct mm_idx_reader_t:
		pass

	mm_idx_reader_t *mm_idx_reader_open(const char *fn, const mm_idxopt_t *opt, const char *fn_out)
	mm_idx_t *mm_idx_reader_read(mm_idx_reader_t *r, int n_threads)
	void mm_idx_reader_close(mm_idx_reader_t *r)
	void mm_idx_destroy(mm_idx_t *mi)
	void mm_mapopt_update(mm_mapopt_t *opt, const mm_idx_t *mi)

	int mm_idx_index_name(mm_idx_t *mi)
	int mm_ref_analysis_file(const mm_idx_t *idx, const char *fn, const mm_mapopt_t *opt, int n_threads, const mm_ref_analysis_opt_t *aopt, mm_ref_analysis_result_t **out)
	int mm_ref_analysis_write(const mm_ref_analysis_result_t *res, const char *prefix)
	void mm_ref_analysis_result_destroy(mm_ref_analysis_result_t *res)

	#
	# Mapping (key struct defined in cmappy.h below)
	#
	ctypedef struct mm_reg1_t:
		pass

	ctypedef struct mm_tbuf_t:
		pass

	mm_tbuf_t *mm_tbuf_init()
	void mm_tbuf_destroy(mm_tbuf_t *b)
	void *mm_tbuf_get_km(mm_tbuf_t *b)
	int mm_gen_cs(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq, int no_iden)
	int mm_gen_MD(void *km, char **buf, int *max_len, const mm_idx_t *mi, const mm_reg1_t *r, const char *seq)

cdef extern from "bseq.h":
	ctypedef struct mm_bseq_file_t:
		pass

	ctypedef struct mm_bseq1_t:
		int l_seq, rid
		char *name, *seq, *qual, *comment

	mm_bseq_file_t *mm_bseq_open(const char *fn)
	void mm_bseq_close(mm_bseq_file_t *fp)
	mm_bseq1_t *mm_bseq_read3(mm_bseq_file_t *fp, int64_t chunk_size, int with_qual, int with_comment, int frag_mode, int *n_)
	int mm_bseq_eof(mm_bseq_file_t *fp)

#
# Helper header (because it is hard to expose mm_reg1_t with Cython)
#
cdef extern from "cmappy.h":
	ctypedef struct mm_hitpy_t:
		const char *ctg
		int32_t ctg_start, ctg_end
		int32_t qry_start, qry_end
		int32_t blen, mlen, NM, ctg_len
		uint8_t mapq, is_primary
		int8_t strand, trans_strand
		int32_t seg_id
		int32_t n_cigar32
		uint32_t *cigar32

	ctypedef struct mappy_map_result_t:
		mm_reg1_t *regs
		int n_regs

	void mm_reg2hitpy(const mm_idx_t *mi, mm_reg1_t *r, mm_hitpy_t *h)
	void mm_free_reg1(mm_reg1_t *r)
	mm_reg1_t *mm_map_aux(const mm_idx_t *mi, const char* seqname, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt)
	mm_reg1_t *mm_map_aux_nogil(const mm_idx_t *mi, const char* seqname, const char *seq1, const char *seq2, int *n_regs, mm_tbuf_t *b, const mm_mapopt_t *opt) nogil
	char *mappy_fetch_seq(const mm_idx_t *mi, const char *name, int st, int en, int *l)
	mm_idx_t *mappy_idx_seq(int w, int k, int is_hpc, int bucket_bits, const char *seq, int l)

	ctypedef struct kstring_t:
		unsigned l, m
		char *s

	ctypedef struct kstream_t:
		pass

	ctypedef struct kseq_t:
		kstring_t name, comment, seq, qual
		int last_char
		kstream_t *f

	kseq_t *mm_fastx_open(const char *fn)
	void mm_fastx_close(kseq_t *ks)
	int kseq_read(kseq_t *seq)

	void mappy_bseq_free_records(mm_bseq1_t *seqs, int n) nogil
	void mappy_batch_result_destroy(mappy_map_result_t *result) nogil
	void mappy_batch_results_free(mappy_map_result_t *results, int n) nogil
	int mappy_map_batch(const mm_idx_t *mi, const mm_mapopt_t *opt, mm_bseq1_t *seqs, int n_seq, int n_threads, mappy_map_result_t *results) nogil
	int mappy_write_sam_header(FILE *fp, const mm_idx_t *idx)
	int mappy_write_batch(FILE *fp, const mm_idx_t *mi, const mm_mapopt_t *opt, mm_bseq1_t *seqs, mappy_map_result_t *results, int n_seq, int output_format) nogil

	char *mappy_revcomp(int l, const uint8_t *seq)
	int mm_verbose_level(int v)
	void mm_reset_timer()

cdef extern from "mmpriv.h":
	ctypedef struct mm_ref_analysis_collect_t:
		pass

	mm_ref_analysis_collect_t *mm_ref_analysis_collect_init(const mm_idx_t *mi, const mm_ref_analysis_opt_t *opt)
	void mm_ref_analysis_collect_destroy(mm_ref_analysis_collect_t *c)
	void mm_ref_analysis_collect_set_total(mm_ref_analysis_collect_t *c, int64_t total_reads)
	void mm_ref_analysis_step(mm_ref_analysis_collect_t *c, const mm_bseq1_t *t, int n_reg, mm_reg1_t *regs)
	int mm_ref_analysis_collect_finish(mm_ref_analysis_collect_t *c, mm_ref_analysis_result_t **out)
