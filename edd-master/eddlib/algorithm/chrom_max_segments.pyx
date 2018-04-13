from libc.stdlib cimport malloc, free, realloc
cimport cython
cimport numpy as np

cdef class Segment:
    cdef readonly int from_idx, to_idx
    cdef readonly float score

    def __init__(self, s, fi, ti):
        self.score = s
        self.from_idx = fi
        self.to_idx = ti


cdef struct segment:
    float I
    float L
    float R
    int lidx
    int ridx

cdef int consider(segment *s, segment *buf, int k):
    cdef int j = k - 1
    while (j >= 0 and buf[j].L >= s.L):
        j -= 1
    if (j == -1 or buf[j].R >= s.R):
        buf[k] = s[0]
        return k + 1
    else:
        s.I = s.R - buf[j].L
        s.L = buf[j].L
        s.lidx = buf[j].lidx
        s.ridx = s.ridx
        return consider(s, buf, j)

cdef enum:
    N = 1000

cdef int max_segments_impl(np.ndarray[np.float64_t] xs,
                           segment **buf, int len_buf):
    cdef float csum = 0
    cdef int k = 0
    cdef int i
    cdef segment s
    cdef float x, Lk, Rk
    for i in range(len(xs)):
        x = xs[i]
        Lk = csum
        Rk = csum + x
        csum += x;
        if (k == len_buf):
            len_buf *= 2
            if (len_buf < N):
                # if len_xs is initially zero
              len_buf = N
            buf[0] = <segment *> realloc(buf[0], len_buf * sizeof(segment))
        if (x > 0):
            s = segment(x, Lk, Rk, i, i)
            k = consider(&s, buf[0], k)
    return k

def max_segments(np.ndarray[np.float64_t] xs):
    cdef:
      segment *buf = NULL
      int buf_len, i
    buf_len = max_segments_impl(xs, &buf, 0)
    res = []
    for i in range(buf_len):
        y = Segment(buf[i].I, buf[i].lidx, buf[i].ridx)
        res.append(y)
    free(buf)
    return res

def maximum_segment(np.ndarray[np.float64_t] xs):
  cdef:
    segment *buf = NULL
    int i, buf_len
    float cur_score, max_score = 0
  buf_len = max_segments_impl(xs, &buf, 0)
  for i in range(buf_len):
    cur_score = buf[i].I
    if cur_score > max_score:
      max_score = cur_score
  free(buf)
  return max_score
