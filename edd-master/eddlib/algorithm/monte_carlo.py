import numpy as np
from chrom_max_segments import maximum_segment
from bisect import bisect_left
import multiprocessing
from statsmodels.sandbox.stats.multicomp import multipletests

class MonteCarlo(object):
    """
    Initialized with a dict of chrom_name and observed scores per chromosome.
    A trial performs a random shuffling and a subsequent count of segment sizes.

    Shuffles across chromosomes 
    """

    def __init__(self, observed_data):
        """
        Arguments:
        - `observed_data`: {<chrom-name>: array of scores}
        """

        self.names = observed_data.keys()
        self.chrom_lengths = [len(x) for x in observed_data.values()]
        self._mc_arr = np.concatenate(observed_data.values())

    def __generate_permutation(self):
        a = self._mc_arr#.copy()
        np.random.shuffle(a)
        # last interval is implicitly known (rest of array), so it's left out
        split_idx = np.cumsum(self.chrom_lengths[:-1])
        xs = np.split(a, split_idx)
        return dict(zip(self.names, xs))

    def trial(self):
        'Runs a trial and returns the score of the maximum segment'
        p = self.__generate_permutation()
        return max(maximum_segment(xs) for xs in p.values())
        
    def __call__(self, i):
        '''same as self.trial(), but prints a "." to stderr when done to
        indicate progress...'''
        np.random.seed()
        res = self.trial()
        return res
    workers = None 
    @classmethod
    def run_simulation(cls, observed_data, niter=4, nprocs=4):
        mc = cls(observed_data)
        if nprocs > 1:
            if cls.workers is None:
                cls.workers = multiprocessing.Pool(nprocs)
                # async makes keyboard interrupt work
            xs = cls.workers.map_async(mc, range(niter)).get(99999999)
        else:
            xs = [mc(i) for i in range(niter)]
        return np.sort(xs)


def fdr_qvals(obs, mc):
    '''
    compute pvalues and fdr correct them.
    '''
    def compute_pvalues(obs, mc):
        for o in obs:
            h0_greaterequal = len(mc) - bisect_left(mc, o)
            pvalue = float(h0_greaterequal + 1) / (len(mc) + 1)
            yield pvalue

    obs = np.sort(obs)[::-1]
    mc = np.sort(mc)
    pvals = list(compute_pvalues(obs, mc))
    qvals = list(multipletests(pvals, method='fdr_bh')[1])
    return dict(pvals=pvals, qvals=qvals)

