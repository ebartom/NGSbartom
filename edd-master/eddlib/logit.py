import numpy as np
from logbook import Logger
from statsmodels.stats.proportion import proportion_confint

log = Logger(__name__)

def logit(xs):
    return np.log(xs) - np.log(1 - xs)

def get_medians(df):
    neg = np.median(df.ix[df.score < 0].score)
    pos = np.median(df.ix[df.score > 0].score)
    return neg, pos

def ci_for_df(odf, ci_method, ci_min=0.25):
    df = odf.copy()
    df['tot_reads'] = df.ip + df.input
    df['avg'] = df.ip / df.tot_reads.astype(float)
    assert ci_method in ('normal', 'agresti_coull',
                         'beta', 'wilson', 'jeffrey')
    df['ci_low'], df['ci_high'] = proportion_confint(df.ip, df.tot_reads,
                                                     method=ci_method)
    df['ci_diff'] = df.ci_high - df.ci_low
    has_ip_and_ctr_reads = np.logical_and(df.ip > 0, df.input > 0)
    small_CI = np.logical_and(df.ci_diff < ci_min, has_ip_and_ctr_reads)
    df['score'] = logit(df.ix[small_CI].avg)
    return df

def get_nib_ratio(df):
    nbins_with_reads = (df.tot_reads > 0).sum()
    nbins_ok = len(df.score.dropna())
    return float(nbins_with_reads - nbins_ok) / nbins_with_reads

def extrapolate_low_info_bins(odf):
    df = odf.copy()
    median_neg, median_pos = get_medians(df.dropna())
    df.ix[np.isnan(df.score), 'score'] = median_neg
    return df
    
def neg_score_scale(odf, scale):
    df = odf.copy()
    df.ix[df.score < 0, 'score'] *= scale
    return df
