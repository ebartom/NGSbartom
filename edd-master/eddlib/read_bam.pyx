from libc.stdlib cimport calloc, free
import numpy as np
cimport numpy as np
from pysam.libchtslib cimport bam1_t, bam1_core_t
from pysam.libcsamfile cimport Samfile


cdef np.ndarray[np.float64_t] \
    agg_n(np.ndarray[np.float64_t] xs, int n):
  cdef:
    np.ndarray[np.float64_t] ys
    int new_len
    int i, k
  new_len = len(xs) / n
  if len(xs) % n > 0:
    new_len += 1
  ys = np.zeros(new_len, dtype=np.float64)
  k = 0
  for i in range(len(xs)):
    ys[k] += xs[i]
    if i % n == n - 1:
      k += 1
  return ys

def aggregate_every_n_bins(bins_by_chrom, n):
  return {k:agg_n(xs, n) for k, xs in bins_by_chrom.items()}

def read_bam_into_bins(chrom_sizes, bin_size, bam_filename):
  b = BamCounter(chrom_sizes, bam_filename, bin_size)
  b.process_bam()
  return b.chrom_bins

cdef class BamCounter:
  cdef:
    object chrom_sizes
    public object chrom_bins
    double **cb
    size_t cb_len
    size_t *cnbins
    Samfile fp
    size_t bin_size

  def __cinit__(self, chrom_sizes, bam_filename, bin_size):
    self.fp = Samfile(bam_filename)
    self.cb_len = self.fp.nreferences
    self.cb = <double**> calloc(self.cb_len, sizeof(double*));
    self.cnbins = <size_t*> calloc(self.cb_len, sizeof(size_t))

  def __init__(self, chrom_sizes, bam_filename, bin_size):
    self.chrom_sizes = chrom_sizes
    self.bin_size = bin_size
    self.chrom_bins = {}
    cdef np.ndarray[np.float64_t] a

    # add array bin pointers to fast lookup array
    for i in range(self.cb_len):
      chrom_name = self.fp.getrname(i)
      if not chrom_name in self.chrom_sizes:
        self.cb[i] = NULL
        self.cnbins[i] = 0
      else:
        chrom_size = self.chrom_sizes[chrom_name]
        nbins = int(chrom_size / bin_size)
        if chrom_size % bin_size != 0:
          nbins += 1
        a = np.zeros(nbins, dtype=np.float64)
        self.chrom_bins[chrom_name] = a
        self.cb[i] = <double *> a.data
        self.cnbins[i] = nbins

  def process_bam(self):
    cdef:
      bam1_t *b
      size_t start, cov, sidx
      size_t nreads = 0
      size_t start_to_bin_end
      double r
    while self.fp.cnext() >= 0:
      b = self.fp.getCurrent()
      if (b.core.flag & 0x4 or
          self.cb[b.core.tid] == NULL or
          b.core.pos == 0):
        continue
      nreads += 1
      start = b.core.pos - 1 # bam 1-based, bed is 0-based
      sidx = start / self.bin_size
      start_to_bin_end = self.bin_size - (start % self.bin_size)
      if b.core.l_qseq < start_to_bin_end:
        # the whole sequence fits within a bin
        assert(self.cnbins[b.core.tid] > sidx)
        self.cb[b.core.tid][sidx] += 1
      else:
        # the sequence spans two bins
        r = start_to_bin_end / <double>b.core.l_qseq;
        self.cb[b.core.tid][sidx] += r
        if (self.cnbins[b.core.tid] > (sidx + 1)):
          # only add to next bin if present (ok to drop if not)
          self.cb[b.core.tid][sidx + 1] += (1 - r)
    return nreads;

  def __dealloc__(self):
    free(self.cb)
    free(self.cnbins)
