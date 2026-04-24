from libc.stdint cimport uint8_t, int8_t
from libc.stdlib cimport free, calloc
from libc.stdio cimport FILE, fopen, fclose, fflush, stdout
cimport cmappy
import sys
import time

__version__ = '2.30'

cdef int MM_I_COMPACT_REFS = 0x8
cdef long long MM_F_CIGAR = 0x004
cdef long long MM_F_OUT_SAM = 0x008
cdef long long MM_F_OUT_CG = 0x020
cdef long long MM_F_OUT_CS = 0x040
cdef long long MM_F_OUT_MD = 0x1000000

cdef int MAPPY_OUTPUT_LEGACY = 0
cdef int MAPPY_OUTPUT_PAF = 1
cdef int MAPPY_OUTPUT_SAM = 2

cmappy.mm_reset_timer()


def _format_duration(seconds):
    seconds = max(0.0, float(seconds))
    if seconds < 60:
        return f"{seconds:.1f}s"
    minutes, sec = divmod(int(round(seconds)), 60)
    if minutes < 60:
        return f"{minutes}m{sec:02d}s"
    hours, minutes = divmod(minutes, 60)
    return f"{hours}h{minutes:02d}m{sec:02d}s"


def _emit_progress(processed, mapped, total_reads, start_time, batch_no=None, verbose=True, final=False):
    elapsed = time.perf_counter() - start_time
    throughput = processed / elapsed if elapsed > 0 else 0.0
    if total_reads and throughput > 0:
        remaining = max(0, total_reads - processed)
        eta = _format_duration(remaining / throughput)
        progress = f"{processed}/{total_reads}"
    else:
        eta = "unknown"
        progress = str(processed)
    if final:
        prefix = "[mappy] completed"
    elif verbose and batch_no is not None:
        prefix = f"[mappy] batch {batch_no}"
    else:
        prefix = "[mappy] progress"
    print(
        f"{prefix}: processed {progress} reads; mapped {mapped}; elapsed {_format_duration(elapsed)}; "
        f"throughput {throughput:.2f} reads/s; ETA {eta}",
        file=sys.stderr,
        flush=True,
    )


cdef inline int _output_format_code(object output_format):
    cdef str fmt
    if output_format is None:
        return MAPPY_OUTPUT_LEGACY
    fmt = output_format.decode() if isinstance(output_format, bytes) else str(output_format)
    fmt = fmt.lower()
    if fmt == 'legacy':
        return MAPPY_OUTPUT_LEGACY
    if fmt == 'paf':
        return MAPPY_OUTPUT_PAF
    if fmt == 'sam':
        return MAPPY_OUTPUT_SAM
    raise ValueError(f"unsupported output format: {output_format}")


cdef class Alignment:
    cdef int _ctg_len, _r_st, _r_en
    cdef int _q_st, _q_en
    cdef int _NM, _mlen, _blen
    cdef int8_t _strand, _trans_strand
    cdef uint8_t _mapq, _is_primary
    cdef int _seg_id
    cdef _ctg, _cigar, _cs, _MD

    def __cinit__(self, ctg, cl, cs, ce, strand, qs, qe, mapq, cigar, is_primary, mlen, blen, NM, trans_strand, seg_id, cs_str, MD_str):
        self._ctg = ctg if isinstance(ctg, str) else ctg.decode()
        self._ctg_len, self._r_st, self._r_en = cl, cs, ce
        self._strand, self._q_st, self._q_en = strand, qs, qe
        self._NM, self._mlen, self._blen = NM, mlen, blen
        self._mapq = mapq
        self._cigar = cigar
        self._is_primary = is_primary
        self._trans_strand = trans_strand
        self._seg_id = seg_id
        self._cs = cs_str
        self._MD = MD_str

    @property
    def ctg(self): return self._ctg

    @property
    def ctg_len(self): return self._ctg_len

    @property
    def r_st(self): return self._r_st

    @property
    def r_en(self): return self._r_en

    @property
    def strand(self): return self._strand

    @property
    def trans_strand(self): return self._trans_strand

    @property
    def blen(self): return self._blen

    @property
    def mlen(self): return self._mlen

    @property
    def NM(self): return self._NM

    @property
    def is_primary(self): return (self._is_primary != 0)

    @property
    def q_st(self): return self._q_st

    @property
    def q_en(self): return self._q_en

    @property
    def mapq(self): return self._mapq

    @property
    def cigar(self): return self._cigar

    @property
    def read_num(self): return self._seg_id + 1

    @property
    def cs(self): return self._cs

    @property
    def MD(self): return self._MD

    @property
    def cigar_str(self):
        return "".join(map(lambda x: str(x[0]) + 'MIDNSHP=XB'[x[1]], self._cigar))

    def __str__(self):
        if self._strand > 0: strand = '+'
        elif self._strand < 0: strand = '-'
        else: strand = '?'
        if self._is_primary != 0: tp = 'tp:A:P'
        else: tp = 'tp:A:S'
        if self._trans_strand > 0: ts = 'ts:A:+'
        elif self._trans_strand < 0: ts = 'ts:A:-'
        else: ts = 'ts:A:.'
        a = [str(self._q_st), str(self._q_en), strand, self._ctg, str(self._ctg_len), str(self._r_st), str(self._r_en),
            str(self._mlen), str(self._blen), str(self._mapq), tp, ts, "cg:Z:" + self.cigar_str]
        if self._cs != "": a.append("cs:Z:" + self._cs)
        if self._MD != "": a.append("MD:Z:" + self._MD)
        return "\t".join(a)


cdef class ThreadBuffer:
    cdef cmappy.mm_tbuf_t *_b

    def __cinit__(self):
        self._b = cmappy.mm_tbuf_init()

    def __dealloc__(self):
        cmappy.mm_tbuf_destroy(self._b)


cdef class Aligner:
    cdef cmappy.mm_idx_t *_idx
    cdef cmappy.mm_idxopt_t idx_opt
    cdef cmappy.mm_mapopt_t map_opt

    def __cinit__(self, fn_idx_in=None, preset=None, k=None, w=None, min_cnt=None, min_chain_score=None, min_dp_score=None, bw=None, bw_long=None, best_n=None, n_threads=3, fn_idx_out=None, max_frag_len=None, extra_flags=None, seq=None, scoring=None, sc_ambi=None, max_chain_skip=None, compact_repeats=False, compact_k=27, compact_ratio=0.20):
        self._idx = NULL
        cmappy.mm_set_opt(NULL, &self.idx_opt, &self.map_opt)
        if preset is not None:
            cmappy.mm_set_opt(str.encode(preset), &self.idx_opt, &self.map_opt)
        self.map_opt.flag |= MM_F_CIGAR
        self.idx_opt.batch_size = 0x7fffffffffffffff
        if k is not None: self.idx_opt.k = k
        if w is not None: self.idx_opt.w = w
        if min_cnt is not None: self.map_opt.min_cnt = min_cnt
        if min_chain_score is not None: self.map_opt.min_chain_score = min_chain_score
        if min_dp_score is not None: self.map_opt.min_dp_max = min_dp_score
        if bw is not None: self.map_opt.bw = bw
        if bw_long is not None: self.map_opt.bw_long = bw_long
        if best_n is not None: self.map_opt.best_n = best_n
        if max_frag_len is not None: self.map_opt.max_frag_len = max_frag_len
        if extra_flags is not None: self.map_opt.flag |= extra_flags
        if scoring is not None and len(scoring) >= 4:
            self.map_opt.a, self.map_opt.b = scoring[0], scoring[1]
            self.map_opt.q, self.map_opt.e = scoring[2], scoring[3]
            self.map_opt.q2, self.map_opt.e2 = self.map_opt.q, self.map_opt.e
            if len(scoring) >= 6:
                self.map_opt.q2, self.map_opt.e2 = scoring[4], scoring[5]
                if len(scoring) >= 7:
                    self.map_opt.sc_ambi = scoring[6]
        if sc_ambi is not None: self.map_opt.sc_ambi = sc_ambi
        if max_chain_skip is not None: self.map_opt.max_chain_skip = max_chain_skip
        self.idx_opt.compact_k = compact_k
        self.idx_opt.compact_ratio = compact_ratio
        if compact_repeats: self.idx_opt.flag |= MM_I_COMPACT_REFS

        cdef cmappy.mm_idx_reader_t *r

        if seq is None:
            if fn_idx_out is None:
                r = cmappy.mm_idx_reader_open(str.encode(fn_idx_in), &self.idx_opt, NULL)
            else:
                r = cmappy.mm_idx_reader_open(str.encode(fn_idx_in), &self.idx_opt, str.encode(fn_idx_out))
            if r is not NULL:
                self._idx = cmappy.mm_idx_reader_read(r, n_threads)
                cmappy.mm_idx_reader_close(r)
                cmappy.mm_mapopt_update(&self.map_opt, self._idx)
                cmappy.mm_idx_index_name(self._idx)
        else:
            self._idx = cmappy.mappy_idx_seq(self.idx_opt.w, self.idx_opt.k, self.idx_opt.flag&1, self.idx_opt.bucket_bits, str.encode(seq), len(seq))
            cmappy.mm_mapopt_update(&self.map_opt, self._idx)
            self.map_opt.mid_occ = 1000

    def __dealloc__(self):
        if self._idx is not NULL:
            cmappy.mm_idx_destroy(self._idx)

    def __bool__(self):
        return (self._idx != NULL)

    cdef list _collect_alignments(self, cmappy.mm_reg1_t *regs, int n_regs, bytes seq1, bytes seq2, bint has_seq2, bint cs, bint MD, cmappy.mm_tbuf_t *buf):
        cdef cmappy.mm_hitpy_t h
        cdef list out = []
        cdef list cigar
        cdef int i, k
        cdef int l_cs_str = 0
        cdef int m_cs_str = 0
        cdef char *cs_str = NULL
        cdef char *seq1_ptr = seq1
        cdef char *seq2_ptr = NULL
        cdef char *cur_seq = NULL
        cdef void *km = cmappy.mm_tbuf_get_km(buf)
        cdef str _cs, _MD

        if has_seq2:
            seq2_ptr = seq2

        try:
            for i in range(n_regs):
                cmappy.mm_reg2hitpy(self._idx, &regs[i], &h)
                cigar = []
                _cs = ''
                _MD = ''
                for k in range(h.n_cigar32):
                    c = h.cigar32[k]
                    cigar.append([c >> 4, c & 0xf])
                if cs or MD:
                    cur_seq = seq2_ptr if h.seg_id > 0 and has_seq2 else seq1_ptr
                    if cs:
                        l_cs_str = cmappy.mm_gen_cs(km, &cs_str, &m_cs_str, self._idx, &regs[i], cur_seq, 1)
                        _cs = (<bytes>cs_str[:l_cs_str]).decode()
                    if MD:
                        l_cs_str = cmappy.mm_gen_MD(km, &cs_str, &m_cs_str, self._idx, &regs[i], cur_seq)
                        _MD = (<bytes>cs_str[:l_cs_str]).decode()
                out.append(Alignment(h.ctg, h.ctg_len, h.ctg_start, h.ctg_end, h.strand, h.qry_start, h.qry_end, h.mapq, cigar, h.is_primary, h.mlen, h.blen, h.NM, h.trans_strand, h.seg_id, _cs, _MD))
            return out
        finally:
            free(cs_str)

    def map(self, seq, seq2=None, name=None, buf=None, cs=False, MD=False, max_frag_len=None, extra_flags=None):
        cdef cmappy.mm_reg1_t *regs = NULL
        cdef ThreadBuffer b
        cdef int n_regs = 0
        cdef cmappy.mm_mapopt_t map_opt
        cdef list alignments = []
        cdef bytes _seq
        cdef bytes _seq2 = b''
        cdef bytes _name

        if self._idx == NULL:
            return
        if ((self.map_opt.flag & MM_F_CIGAR) and (self._idx.flag & 2)):
            return
        map_opt = self.map_opt
        if max_frag_len is not None: map_opt.max_frag_len = max_frag_len
        if extra_flags is not None: map_opt.flag |= extra_flags

        if buf is None: b = ThreadBuffer()
        else: b = buf

        _seq = seq if isinstance(seq, bytes) else seq.encode()
        if name is not None:
            _name = name if isinstance(name, bytes) else name.encode()

        if seq2 is None:
            if name is None:
                regs = cmappy.mm_map_aux(self._idx, NULL, _seq, NULL, &n_regs, b._b, &map_opt)
            else:
                regs = cmappy.mm_map_aux(self._idx, _name, _seq, NULL, &n_regs, b._b, &map_opt)
        else:
            _seq2 = seq2 if isinstance(seq2, bytes) else seq2.encode()
            if name is None:
                regs = cmappy.mm_map_aux(self._idx, NULL, _seq, _seq2, &n_regs, b._b, &map_opt)
            else:
                regs = cmappy.mm_map_aux(self._idx, _name, _seq, _seq2, &n_regs, b._b, &map_opt)

        try:
            alignments = self._collect_alignments(regs, n_regs, _seq, _seq2, seq2 is not None, cs, MD, b._b)
        finally:
            i = 0
            while i < n_regs:
                cmappy.mm_free_reg1(&regs[i])
                i += 1
            free(regs)

        for aln in alignments:
            yield aln

    def map_batch(self, records, cs=False, MD=False, max_frag_len=None, extra_flags=None, n_threads=1):
        cdef list items = list(records)
        cdef Py_ssize_t n_seq = len(items)
        cdef cmappy.mm_bseq1_t *seqs = NULL
        cdef cmappy.mappy_map_result_t *results = NULL
        cdef cmappy.mm_mapopt_t map_opt
        cdef list names = []
        cdef list seq_bufs = []
        cdef ThreadBuffer buf = ThreadBuffer()
        cdef list out = []
        cdef Py_ssize_t i
        cdef object rec
        cdef object name
        cdef object seq
        cdef bytes name_b
        cdef bytes seq_b
        cdef int rc
        cdef int thread_count

        if self._idx == NULL:
            return []
        if ((self.map_opt.flag & MM_F_CIGAR) and (self._idx.flag & 2)):
            return []
        if n_seq == 0:
            return []

        map_opt = self.map_opt
        if max_frag_len is not None: map_opt.max_frag_len = max_frag_len
        if extra_flags is not None: map_opt.flag |= extra_flags

        seqs = <cmappy.mm_bseq1_t*>calloc(n_seq, sizeof(cmappy.mm_bseq1_t))
        results = <cmappy.mappy_map_result_t*>calloc(n_seq, sizeof(cmappy.mappy_map_result_t))
        if seqs == NULL or results == NULL:
            free(seqs)
            free(results)
            raise MemoryError('unable to allocate batch mapping buffers')

        try:
            names = [None] * n_seq
            seq_bufs = [None] * n_seq
            for i in range(n_seq):
                rec = items[i]
                if not isinstance(rec, (tuple, list)) or len(rec) < 2:
                    raise ValueError('each batch record must be a tuple/list of at least (name, seq)')
                name = rec[0]
                seq = rec[1]
                name_b = name if isinstance(name, bytes) else str(name).encode()
                seq_b = seq if isinstance(seq, bytes) else str(seq).encode()
                names[i] = name_b
                seq_bufs[i] = seq_b
                seqs[i].name = name_b
                seqs[i].seq = seq_b
                seqs[i].l_seq = len(seq_b)
                seqs[i].qual = NULL
                seqs[i].comment = NULL
                seqs[i].rid = i

            if n_threads is None or n_threads < 1:
                n_threads = 1
            thread_count = <int>n_threads
            with nogil:
                rc = cmappy.mappy_map_batch(self._idx, &map_opt, seqs, <int>n_seq, thread_count, results)
            if rc != 0:
                raise RuntimeError('threaded batch mapping failed')

            out = []
            for i in range(n_seq):
                out.append(self._collect_alignments(results[i].regs, results[i].n_regs, seq_bufs[i], b'', 0, cs, MD, buf._b))
            return out
        finally:
            if results != NULL:
                with nogil:
                    cmappy.mappy_batch_results_free(results, <int>n_seq)
                results = NULL
            free(seqs)

    def map_file(self, query_path, output_path='-', output_format='legacy', n_threads=1, cs=False, MD=False, max_frag_len=None, extra_flags=None, verbose=True, total_reads=0, progress_every=1000, batch_bases=33554432, batch_reads=4096):
        cdef cmappy.mm_bseq_file_t *reader = NULL
        cdef cmappy.mm_bseq1_t *seqs = NULL
        cdef cmappy.mm_bseq1_t *chunk_seqs
        cdef cmappy.mappy_map_result_t *results = NULL
        cdef cmappy.mm_mapopt_t map_opt
        cdef FILE *out_fp = stdout
        cdef bytes query_b
        cdef bytes output_b
        cdef int output_code = _output_format_code(output_format)
        cdef int n_batch = 0
        cdef int chunk_n
        cdef int offset
        cdef int rc
        cdef int thread_count
        cdef int write_rc
        cdef int with_qual
        cdef int chunk_mapped
        cdef int i
        cdef long processed = 0
        cdef long mapped = 0
        cdef long last_reported = 0
        cdef long batch_no = 0
        cdef object start_time
        cdef object stats

        if self._idx == NULL:
            raise ValueError('index is not loaded')
        if ((self.map_opt.flag & MM_F_CIGAR) and (self._idx.flag & 2)):
            raise ValueError('alignment output requires reference sequences in the index')

        if n_threads is None or n_threads < 1:
            n_threads = 1
        thread_count = <int>n_threads
        if batch_reads is None or batch_reads < 1:
            batch_reads = 4096
        if batch_bases is None or batch_bases < 1:
            batch_bases = 33554432
        if progress_every is None or progress_every < 1:
            progress_every = 1000

        map_opt = self.map_opt
        if max_frag_len is not None: map_opt.max_frag_len = max_frag_len
        if extra_flags is not None: map_opt.flag |= extra_flags
        if output_code == MAPPY_OUTPUT_SAM:
            map_opt.flag |= MM_F_OUT_SAM
        elif output_code == MAPPY_OUTPUT_PAF:
            map_opt.flag |= MM_F_OUT_CG
        if cs:
            map_opt.flag |= MM_F_OUT_CS
        if MD:
            map_opt.flag |= MM_F_OUT_MD

        query_b = query_path if isinstance(query_path, bytes) else str(query_path).encode()
        reader = cmappy.mm_bseq_open(query_b)
        if reader is NULL:
            raise OSError(f"failed to open query file: {query_path}")

        if output_path not in (None, '-'):
            output_b = output_path if isinstance(output_path, bytes) else str(output_path).encode()
            out_fp = fopen(output_b, b"wb")
            if out_fp == NULL:
                cmappy.mm_bseq_close(reader)
                raise OSError(f"failed to open output file: {output_path}")

        with_qual = 1 if output_code == MAPPY_OUTPUT_SAM else 0
        start_time = time.perf_counter()

        try:
            if output_code == MAPPY_OUTPUT_SAM:
                rc = cmappy.mappy_write_sam_header(out_fp, self._idx)
                if rc != 0:
                    raise RuntimeError('failed to write SAM header')

            while True:
                n_batch = 0
                seqs = cmappy.mm_bseq_read3(reader, batch_bases, with_qual, 0, 0, &n_batch)
                if seqs == NULL or n_batch <= 0:
                    break
                try:
                    offset = 0
                    while offset < n_batch:
                        chunk_n = batch_reads if n_batch - offset > batch_reads else n_batch - offset
                        chunk_seqs = seqs + offset
                        results = <cmappy.mappy_map_result_t*>calloc(chunk_n, sizeof(cmappy.mappy_map_result_t))
                        if results == NULL:
                            raise MemoryError('unable to allocate threaded mapping results')
                        try:
                            with nogil:
                                rc = cmappy.mappy_map_batch(self._idx, &map_opt, chunk_seqs, chunk_n, thread_count, results)
                            if rc != 0:
                                raise RuntimeError('threaded file mapping failed')
                            with nogil:
                                write_rc = cmappy.mappy_write_batch(out_fp, self._idx, &map_opt, chunk_seqs, results, chunk_n, output_code)
                            if write_rc < 0:
                                raise RuntimeError('failed to write mapped output batch')

                            chunk_mapped = 0
                            for i in range(chunk_n):
                                if results[i].n_regs > 0:
                                    chunk_mapped += 1
                            processed += chunk_n
                            mapped += chunk_mapped
                            batch_no += 1

                            if verbose:
                                if processed - last_reported >= progress_every or (total_reads and processed >= total_reads):
                                    _emit_progress(processed, mapped, total_reads, start_time, batch_no=batch_no, verbose=True)
                                    last_reported = processed
                            else:
                                if processed - last_reported >= progress_every or (total_reads and processed >= total_reads):
                                    _emit_progress(processed, mapped, total_reads, start_time, verbose=False)
                                    last_reported = processed
                        finally:
                            if results != NULL:
                                with nogil:
                                    cmappy.mappy_batch_results_free(results, chunk_n)
                                results = NULL
                        offset += chunk_n
                finally:
                    with nogil:
                        cmappy.mappy_bseq_free_records(seqs, n_batch)
                    seqs = NULL
        finally:
            if seqs != NULL:
                with nogil:
                    cmappy.mappy_bseq_free_records(seqs, n_batch)
            if out_fp != NULL:
                fflush(out_fp)
                if out_fp != stdout:
                    fclose(out_fp)
            if reader != NULL:
                cmappy.mm_bseq_close(reader)

        if processed != last_reported or processed == 0:
            _emit_progress(processed, mapped, total_reads, start_time, verbose=verbose, final=True)
        else:
            _emit_progress(processed, mapped, total_reads, start_time, verbose=verbose, final=True)

        stats = {
            'processed': processed,
            'mapped': mapped,
            'elapsed_sec': time.perf_counter() - start_time,
            'reads_per_sec': (processed / (time.perf_counter() - start_time)) if processed > 0 and (time.perf_counter() - start_time) > 0 else 0.0,
            'output_format': output_format,
            'threads': n_threads,
        }
        return stats

    def seq(self, str name, int start=0, int end=0x7fffffff):
        cdef int l
        cdef char *s
        if self._idx == NULL: return
        if ((self.map_opt.flag & MM_F_CIGAR) and (self._idx.flag & 2)): return
        s = cmappy.mappy_fetch_seq(self._idx, name.encode(), start, end, &l)
        if l == 0: return None
        r = s[:l] if isinstance(s, str) else s[:l].decode()
        free(s)
        return r

    @property
    def k(self): return self._idx.k

    @property
    def w(self): return self._idx.w

    @property
    def n_seq(self): return self._idx.n_seq

    @property
    def seq_names(self):
        cdef char *p
        if self._idx == NULL: return
        sn = []
        for i in range(self._idx.n_seq):
            p = self._idx.seq[i].name
            s = p if isinstance(p, str) else p.decode()
            sn.append(s)
        return sn


def fastx_read(fn, read_comment=False):
    cdef cmappy.kseq_t *ks
    ks = cmappy.mm_fastx_open(str.encode(fn))
    if ks is NULL: return None
    while cmappy.kseq_read(ks) >= 0:
        if ks.qual.l > 0: qual = ks.qual.s if isinstance(ks.qual.s, str) else ks.qual.s.decode()
        else: qual = None
        name = ks.name.s if isinstance(ks.name.s, str) else ks.name.s.decode()
        seq = ks.seq.s if isinstance(ks.seq.s, str) else ks.seq.s.decode()
        if read_comment:
            if ks.comment.l > 0: comment = ks.comment.s if isinstance(ks.comment.s, str) else ks.comment.s.decode()
            else: comment = None
            yield name, seq, qual, comment
        else:
            yield name, seq, qual
    cmappy.mm_fastx_close(ks)


def revcomp(seq):
    l = len(seq)
    bseq = seq if isinstance(seq, bytes) else seq.encode()
    cdef char *s = cmappy.mappy_revcomp(l, bseq)
    r = s[:l] if isinstance(s, str) else s[:l].decode()
    free(s)
    return r


def verbose(v=None):
    if v is None: v = -1
    return cmappy.mm_verbose_level(v)
