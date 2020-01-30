'''
Copy number variation (CNV) correction module
Author: Alexander Wu, Wei Li
'''
import numpy as np
import scipy

def read_CNVdata(CN_file,cell_list):
    '''reads a file contaning a matrix of copy number data and filters out
    copy number data for inputted set of desired cell lines'''

    ndarr = np.genfromtxt(CN_file,delimiter='\t',names=True,dtype=None)

    # dictionary of cell line indices from copy number data matrix
    cell_dict = {key:val for (val,key) in enumerate(ndarr.dtype.names)}
    # dictionary of gene symbol indices from copy number data matrix
    gene_dict = {key:val for (val,key) in enumerate(ndarr['SYMBOL'])}

    # identify matches in list of desired cell lines and cell lines in CN data
    inds = [] 
    matches = []
    for cell in cell_list:
        for name in cell_dict:
            if cell.upper() == name.upper():
                inds.append(cell_dict[name])
                matches.append(cell)

    # convert ndarray into array with float values (instead of string)
    # NOTE: also adjusting log2(CN) to CN by exponentiating
    arr = np.asarray([list(row) for row in ndarr])[:,inds]
    arr = arr.astype(np.float)
    arr = 2**arr

    # dictionary of cell line indices from filtered array of CN data
    new_cell_dict = {key:val for (val,key) in enumerate(matches)}

    return (arr,new_cell_dict,gene_dict)

def match_sgrnaCN(score_list,gene_list,CN_arr,CN_genedict):
    '''filter out scores and CN data for genes that are represented by both the
    sgRNA dataset and CNV dataset'''

    # isolate data for regression model
    score_vals = []
    CN_vals = []
    for i in range(len(score_list)):
        score = score_list[i]
        gene = gene_list[i]
        # retain sgRNA score if its target gene is represented in CNV data
        if gene in CN_genedict:
            score_vals.append(score)
            CN_ind = CN_genedict[gene]
            CN_vals.append(CN_arr[CN_ind,0])

    return score_vals,CN_vals
    
###############################################################################
### PIECEWISE NORMALIZATION ###################################################
###############################################################################

def betascore_piecewisenorm(allgenedict,betalabels,CN_arr,CN_celldict,CN_genedict):
    from scipy import stats
    from scipy import optimize

    # identify cell lines on which to CNV normalize 
    # (exclude cell lines not represented in the CNV dataset)
    cell = betalabels[1:]
    cellinds = [i for i in range(len(cell)) if cell[i] in CN_celldict]

    # for each cell line/sample represented
    for i in cellinds:
        # collect BetaScore values and gene list
        score_vals = []
        CN_vals = []
        gene_list = []
        for (tgid,tginst) in allgenedict.iteritems():
            gene = tginst.prefix
            gene_list.append(gene)
            # retain sgRNA score if its target gene is represented in CNV data
            if gene in CN_genedict:
                # only select BetaScore, CNV data of specific cell line
                start = len(tginst.beta_estimate)-len(cell)
                score_vals.append(tginst.beta_estimate[start+i])
                CN_ind = CN_genedict[gene]
                CN_vals.append(CN_arr[CN_ind,CN_celldict[cell[i]]])

        # convert lists to arrays
        score_vals = np.array(score_vals)
        CN_vals = np.array(CN_vals)

        # breakpoint linear model function
        def linreg_bp(bp):
        # regression model for genes w/ CN <= breakpoint
            slope,intercept,r_value,p_value,std_err = \
                stats.linregress(CN_vals[CN_vals<=bp],score_vals[CN_vals<=bp])
            return (slope,intercept)

        # least squares function
        def leastsq_bp(bp):
            (slope,intercept) = linreg_bp(bp)
            # before breakpoint
            sse = sum((score_vals[CN_vals<=bp]-(intercept+\
                    slope*CN_vals[CN_vals<=bp]))**2)
            # after breakpoint
            sse += sum((score_vals[CN_vals>bp]-(intercept+slope*bp))**2)

            return sse

        # find optimal breakpoint
        opt_bp = optimize.minimize(leastsq_bp,2,bounds=((1,max(CN_vals)),))
        opt_bp = opt_bp.x[0]
        # create piecewise linear regression model
        (slope,intercept) = linreg_bp(opt_bp)

        for gene in allgenedict:

            # normalize if sgRNA target gene represented in CNV data
            if gene in CN_genedict:

                # identify approx. score given CN of gene
                # (by applying piecewise linear regression model)
                CN = CN_arr[CN_genedict[gene],CN_celldict[cell[i]]]
                if CN>=opt_bp:
                    est_score = intercept + slope*opt_bp
                else:
                    est_score = intercept + slope*CN

                # estimated score for CN=1
                mu_score = intercept + slope

                # normalize and update scores ***
                start = len(allgenedict[gene].beta_estimate)-len(cell)
                allgenedict[gene].beta_estimate[start+i] += (mu_score-est_score)

def sgRNAscore_piecewisenorm(score_list,gene_list,CN_arr,CN_genedict):
    '''normalizes sgRNA scores for each sgRNA with respect to CN using a 
    piecewise linear regression model'''

    from scipy import optimize
    from scipy import stats

    norm_list = score_list[:]

    # isolate data for regression model
    score_vals,CN_vals = match_sgrnaCN(score_list,gene_list,CN_arr,CN_genedict)
    # convert lists to arrays
    score_vals = np.array(score_vals)
    CN_vals = np.array(CN_vals)

    # breakpoint linear model function
    def linreg_bp(bp):
    # regression model for genes w/ CN <= breakpoint
        slope,intercept,r_value,p_value,std_err = \
            stats.linregress(CN_vals[CN_vals<=bp],score_vals[CN_vals<=bp])
        return (slope,intercept)

    # least squares function
    def leastsq_bp(bp):
        (slope,intercept) = linreg_bp(bp)
        # before breakpoint
        sse = sum((score_vals[CN_vals<=bp]-(intercept+\
                slope*CN_vals[CN_vals<=bp]))**2)
        # after breakpoint
        sse += sum((score_vals[CN_vals>bp]-(intercept+slope*bp))**2)

        return sse

    # find optimal breakpoint
    opt_bp = optimize.minimize(leastsq_bp,2,bounds=((1,max(CN_vals)),))
    opt_bp = opt_bp.x[0]
    # create piecewise linear regression model
    (slope,intercept) = linreg_bp(opt_bp)

    for i in range(len(score_list)):
        gene = gene_list[i]

        # normalize if sgRNA target gene represented in CNV data
        if gene in CN_genedict:

            # identify approx. score given CN of gene
            # (by applying piecewise linear regression model)
            CN = CN_arr[CN_genedict[gene],0]
            if CN>=opt_bp:
                est_score = intercept + slope*opt_bp
            else:
                est_score = intercept + slope*CN

            # estimated score for CN=1
            mu_score = intercept + slope

            # normalize scores ***
            norm_list[i] += (mu_score-est_score)

    return norm_list
