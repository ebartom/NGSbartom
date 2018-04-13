import collections
import operator
import os
import logbook
import pybedtools


log = logbook.Logger(__name__)

class UnalignableRegions(object):

    def __init__(self, chrom, start, end):
        self.chrom = chrom
        self.start = start
        self.end = end

    def bigger_than(self, x):
        '''
        returns true if gap comes after param x (no overlap)
        '''
        return x.end <= self.start

    def overlaps(self, x):
        '''
        returns true if the two intervals overlaps
        '''
        return (self.start <= x.start < self.end or
                self.start < x.end <= self.end)

    def __repr__(self):
        return 'UnalignableRegions(%s, %d, %d)' % (self.chrom, self.start, self.end)

def read_file(path):
    res = []
    if os.stat(path).st_size > 0:
        for e in pybedtools.BedTool(path).sort().merge():
            x = UnalignableRegions(e.chrom, e.start, e.end)
            res.append(x)
    else:
        log.notice('unalignable regions file is empty, skipping.')
    tot = sum(x.end - x.start for x in res)
    log.notice('Unalignable regions file read. Got %d regions. Total coverage: %.2fMB' % (
        len(res) ,tot / 1e6))
    return res

def split_on_regions(scores_per_chrom, regions):
    dg = collections.defaultdict(list)
    revdict = {}
    d = {}
    for g in regions:
        dg[g.chrom].append(g)

    for chrom, regions in dg.items():
        if not chrom in scores_per_chrom:
            continue
        regions.sort(key=operator.attrgetter('start'))
        groups = []
        cur_grp = []
        cur_gap = 0
        ndropped = 0
        bins = iter(scores_per_chrom[chrom])
        try:
            while True:
                x = bins.next()
                if regions[cur_gap].bigger_than(x):
                    cur_grp.append(x)
                elif regions[cur_gap].overlaps(x):
                    if len(cur_grp) > 0:
                        groups.append(cur_grp)
                        cur_grp = []
                    ndropped += 1
                else:
                    cur_gap += 1
                    if len(regions) == cur_gap:
                        cur_grp.append(x)
                        cur_grp.extend(bins)
                        break
                    if not regions[cur_gap].overlaps(x):
                        cur_grp.append(x)
        except StopIteration:
            pass
        if len(cur_grp) > 0:
            groups.append(cur_grp)
        names = ['%s_%d' % (chrom, i) for i in range(len(groups))]
        d.update(zip(names, groups))
        for name in names:
            revdict[name] = chrom
    for chrom in (scores_per_chrom.viewkeys() - dg.viewkeys()):
        d[chrom] = scores_per_chrom[chrom]
    return d, revdict
