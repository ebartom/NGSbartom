from collections import namedtuple
import math
from logbook import Logger

log = Logger(__name__)
bed = namedtuple('BedGraph', 'chrom start end score')
bin = namedtuple('Bin', 'chrom start end ip_count input_count')

def save_bin_score_file(df, ratio_file):
    df['chrom start end score'.split()].sort(['chrom', 'start']).to_csv(
        ratio_file, sep='\t', index=False, header=False)

def ci_lower_bound(pos, neg):
    '''
    computes lower bound of 95% confidence interval for true binomial proportion.

    sources:
    http://www.evanmiller.org/how-not-to-sort-by-average-rating.html
    http://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
    '''
    n = pos + neg
    if n == 0:
        return 0
    z = 1.96 # for a 95% confidence interval
    pratio = float(pos) / n
    return (pratio + z*z/(2*n) - z * math.sqrt((pratio*(1-pratio)+z*z/(4*n))/n))/(1+z*z/n)

