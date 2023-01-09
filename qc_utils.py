import pandas as pd
import numpy as np
import scanpy as sc
import matplotlib.pylab as plt
from scipy import stats
import skimage.filters as skif
from sklearn.preprocessing import normalize
import glob

def reorder_AnnData(AnnData, descending = True):
    AnnData.obs['total_counts'] = np.array(AnnData.X.sum(axis=1) )
    if(descending==True):
        AnnData = AnnData[np.argsort(AnnData.obs['total_counts'])[::-1]].copy()
        AnnData.uns['Order'] = 'Descending'
    elif(descending==False):
        AnnData = AnnData[np.argsort(AnnData.obs['total_counts'])[:]].copy()
        AnnData.uns['Order'] = 'Ascending'
    return(AnnData)

def read_dropest(dir_path,reorder=True):
    data_matrix = glob.glob(dir_path+"/*.mtx")[0]
    data_genes = glob.glob(dir_path+"/*features.*")[0]
    data_barcodes = glob.glob(dir_path+"/*barcodes.*")[0]
    adata = sc.read_mtx(data_matrix).T
    adata.var.index = pd.read_csv(data_genes,header=None)[0].values
    adata.obs.index = pd.read_csv(data_barcodes,header=None)[0].values
    adata.obs.index.name = 'Cells'
    adata.var.index.name = 'Genes'
    adata = reorder_AnnData(adata,descending=True)
    adata.raw = adata
    return(adata)

def find_inflection(adata_in,mito_tag = "MT-",run_qc=True):
    if(run_qc):
        print("Calculating QC Metrics")
        adata_in.var['Mitochondrial'] = adata_in.var.index.str.startswith(mito_tag)
        sc.pp.calculate_qc_metrics(adata_in,qc_vars=['Mitochondrial'],use_raw=True,inplace=True)
    tc_ngbc_ratio = adata_in.obs['total_counts']/adata_in.obs['n_genes_by_counts']
    data_cumsum = np.cumsum(adata_in.obs['total_counts'])
    x_vals=np.arange(0,adata_in.n_obs)
    secant_coef=data_cumsum[adata_in.obs.shape[0]-1]/adata_in.obs.shape[0]
    secant_line=secant_coef*x_vals
    secant_dist=abs(data_cumsum-secant_line)
    inflection_points = secant_dist.argsort()[::-1]
    inflection_percentiles = [0,15,30]
    inflection_percentiles_inds = np.percentile(x_vals[x_vals>inflection_points[0]],inflection_percentiles).astype(int)
    color=plt.cm.tab10(np.linspace(0,1,adata_in.n_obs))

    fig = plt.figure( figsize=(20,6))
    ax1 = fig.add_subplot(121)
    ax2 = fig.add_subplot(122)
    ax1.plot(np.array(data_cumsum), label="Cumulative Sum")
    ax1.plot(np.array(secant_dist), label="Secant Distance")
    for percentile in inflection_percentiles_inds:
        ax1.axvline(x=percentile,ymin=0,c=color[percentile],linestyle='--',linewidth=2,label="Inflection point {}".format(percentile))
    ax1.legend()
    ax1.set_xlabel("Cell Rank")
    ax1.set_ylabel("Total Counts per Cell")
    ax1.set_title('Inflection Curve')
        
    ax2.plot(np.sort(tc_ngbc_ratio.values)[::-1],label='Total Counts/N Genes By Counts')
    #ax2.plot(np.sort(adata_in.obs['pct_counts_Mitochondrial'])[::-1]/100)
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_title('Total Counts/N Genes By Counts')
    ax2.set_xlabel("Cell Rank")
    ax2.set_ylabel("Total Counts/N Genes By Counts")
    ax2.legend()
        
    print("Inflection point at {} for {} percentiles of greatest secant distances".format(inflection_percentiles_inds,inflection_percentiles))
    return(inflection_percentiles_inds)

def relative_diversity(adata_in,zscore_minimum = 0,log_scale=True):
    adata_in.obs['n_genes_by_counts_zscore'] = stats.zscore(adata_in.obs['n_genes_by_counts'])
    diversity_cutoff = skif.threshold_minimum(adata_in.obs['n_genes_by_counts_zscore'][adata_in.obs['n_genes_by_counts_zscore'] >= zscore_minimum])
        
    fig = plt.figure(figsize=(4,4))
    ax0 = plt.subplot(111)
    
    hist_plot_div = ax0.hist(adata_in.obs['n_genes_by_counts_zscore'],bins=np.sqrt(adata_in.n_obs).astype(int),log=log_scale)
    #ax0.set_ylim([0,max(hist_plot_div[0])])
    #ax0.set_xlim([0,max(hist_plot_div[1])])
    ax0.axvspan(diversity_cutoff,max(hist_plot_div[1]), alpha=0.25, color='green')
    ax0.vlines(diversity_cutoff,0,max(hist_plot_div[0]),color='green')
    ax0.set_title("Relative Transcript Diversity")
    adata_in.obs['relative_transcript_diversity_threshold'] = adata_in.obs['n_genes_by_counts_zscore']>=diversity_cutoff
    return(adata_in)

def subset_cleanup(adata_in,selection='Cell_Selection'):
    adata_in.X = adata_in.raw.X.copy()
    adata_in = adata_in[adata_in.obs[selection]].copy()
    adata_in = sc.AnnData(adata_in.to_df())
    return(adata_in)
