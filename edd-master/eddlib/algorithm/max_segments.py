import operator
import numpy as np
import collections
import itertools
from eddlib import util 
from chrom_max_segments import max_segments
import unalignable_regions
import logbook
from statsmodels.sandbox.stats.multicomp import multipletests

log = logbook.Logger(__name__)

class GenomeBins(object):
    '''
    Genome wide bins with scores.
    '''

    def __init__(self, chrom_bins):
        '''
        Keyword arguments:
        chrom_bins -- dict of names->bins where a name may span a 
                        complete chromosome or just a part
        '''
        self.chrom_bins = chrom_bins
        self.chrom_scores = {k:np.array([x.score for x in v], dtype=float) 
                for k, v in chrom_bins.items()}

    @classmethod
    def with_unalignable_regions_masked(cls, chrom_bins, unalignable_regions_file):
        if unalignable_regions_file is not None:
            g = unalignable_regions.read_file(unalignable_regions_file)
            chrom_bins, _ = unalignable_regions.split_on_regions(chrom_bins, g)
        return cls(chrom_bins)

    def max_segments(self, filter_trivial=0):
        def chrom_max(chrom):
            xs = max_segments(self.chrom_scores[chrom])
            bins = self.chrom_bins[chrom]
            return [util.bed(bins[x.from_idx].chrom, bins[x.from_idx].start,
                        bins[x.to_idx].end, x.score) for x in xs]
        segments_per_chrom = collections.defaultdict(list)
        for x in itertools.chain.from_iterable(chrom_max(k) for k in self.chrom_bins):
            if x.score > filter_trivial:
                segments_per_chrom[x.chrom].append(x)
        num_intervals = 0
        for xs in segments_per_chrom.values():
            xs.sort(key=operator.attrgetter('start'))
            num_intervals += len(xs)
        if filter_trivial > 0:
            log.notice('Removed trivial intervals with score less than %.4f.' % filter_trivial)
            log.notice('%d intervals (potential peaks) remaining.' % num_intervals)
        return segments_per_chrom

    @classmethod
    def df_as_bins(cls, df, gap_file):
        '''
        converts a already scored df to an object
        containing a dict of bins per chromosome, separated by gaps.
        '''
        chromd = collections.defaultdict(list)
        for _, x in df.iterrows():
            b = util.bed(x['chrom'], x['start'], x['end'],
                    x['score'])
            chromd[b.chrom].append(b)
        return cls.with_unalignable_regions_masked(chromd, gap_file)

    def scale_neg_scores(self, scale):
        return GenomeBins(
            {k:[util.bed(x.chrom, x.start, x.end,
                         x.score if x.score > 0 else x.score * scale)
                for x in xs]
             for k, xs in self.chrom_bins.items()})
    

class IntervalTest(object):

    def __init__(self, max_intervals, mc_res):
        self.max_intervals = max_intervals
        self.mc_res = np.sort(mc_res)
        self._pvalues = None
        self._pvalues = None

    def pvalues(self):
        res = []
        for xs in self.max_intervals.values():
            pos = len(self.mc_res) - np.searchsorted(self.mc_res,
                    [x.score for x in xs])
            pvals = (pos + 1) / float(len(self.mc_res) + 1)
            for x, pval in zip(xs, pvals):
                res.append((x, pval))
        self._pvalues = res
        return res

    def qvalues(self, below=0.1):
        if not self._pvalues:
            self.pvalues()
        pvals = [x[1] for x in self._pvalues]
        qvals = list(multipletests(pvals, method='fdr_bh')[1])
        res = [(q,p,x) for (q, (x,p)) in zip(qvals, self._pvalues)
            if q < below]
        self._qvalues = res
        log.notice('got %d peaks with qvalue below %.2f. From %d possible.' % (
            len(res), below, len(pvals)))
        return res

    @classmethod
    def segments_to_bedstream(cls, segments, bedstream):
        for segment in segments:
            bedstream.write('\t'.join((segment.chrom, str(segment.start),
                                str(segment.end),
                                str(segment.score))))
            bedstream.write('\n')

    def as_bed(self, output_file):
        if self._qvalues is None:
            self.qvalues()
        log.notice('Saving significant peaks.')
        segments = [x for (_, _, x) in self._qvalues]
        segments.sort(key=operator.itemgetter(0,1))
        with open(output_file, 'w') as of:
            self.segments_to_bedstream(segments, of)
        log.notice('Done. Wrote %d peaks.' % len(self._qvalues))

