import omicverse as ov
from omicverse.utils import mde
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd

#Annotation of urinary aTEC using marker genes of adaptive states
# import urine scRNA-seq data (aTEC) of DKD patients generated from this study
atec_adata=ov.read('urine_aTEC_DKD.h5ad')
atec_atec_adata.var_names_make_unique()
atec_adata.obs_names_make_unique()
atec_adata



#prepare marker gene dictionary
res_marker_dict={
    'Epithelial':['SLC12A1','UMOD','PROM2','MUC1','SCNN1A','CDH1','KRT7','AQP3','SLC8A1',],
    'Mesenchymal':['MGP','NID2','THY1','VIM'],
    'Fibrotic':['HAVCR1','VCAM1','AQP1', 'ANPEP','LRP2','ZEB2','SERPINA1','NCAM1','BMP2','BMPR1A','IGFBP6', 'PTTG1','CENPF', 'TOP2A','MKI67'], 

}
	
# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(atec_adata,'leiden')
sc.pl.dotplot(atec_adata, res_marker_dict_adj, 'leiden', 
              dendrogram=True,standard_scale='var',cmap='coolwarm')





#annotate urinary aTEC according to the generated dictionary
cluster2annotation = {
     '1': 'aTEC_LOH_like',
     '2': 'aTEC_PT_like',
	 '6': 'aTEC_DCT_like', 
	 '8': 'aTEC_Mes',
	 '11': 'aTEC_PT_prolf',
     '13': 'aTEC_PT_like',

}

ov.single.scanpy_cellanno_from_dict(atec_adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')									   
ov.utils.embedding(atec_adata,
                   basis='X_umap', #X_mde
                   color=['major_celltype'],  
                   legend_loc='right margin',frameon='small', legend_fontoutline=1, 
                    palette=['#E069A6','#E6B7D2','#941456','#FFDEAD','#D2B48C','#FCBC10','#E069A6',]) 





#compute cell type ratio
ov.utils.plot_cellproportion(adata=atec_adata,celltype_clusters='major_celltype',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))





#GO analysis
#prepare GO library
pathway_dict=ov.utils.geneset_prepare('GO_2021.txt',organism='Human') # Downloaded from https://maayanlab.cloud/Enrichr/#libraries





# calculate and visualize all GO terms enriched in early DKD & advanced DKD 
adata_mer.uns['log1p']['base']=None
sc.tl.rank_genes_groups(atec_adata, 'group', method='t-test',n_genes=100)
res=ov.single.pathway_enrichment(atec_adata,pathways_dict=pathway_dict,organism='Human',cutoff=0.05, logfc_threshold=2, pvalue_type='adjust',
                                     group_by='group',plot=True)





#compute AUCelll score based on all GO terms
adata_aucs=ov.single.pathway_aucell_enrichment(atec_adata,
                                                  pathways_dict=pathway_dict,
                                                  num_workers=8)
adata_aucs.obs=atec_adata[adata_aucs.obs.index].obs
adata_aucs.obsm=atec_adata[adata_aucs.obs.index].obsm
adata_aucs.obsp=atec_adata[adata_aucs.obs.index].obsp
adata_aucs





# calculate all GO terms enriched in each urinary aTEC cluster
adata_mer.uns['log1p']['base']=None
sc.tl.rank_genes_groups(atec_adata, 'major_celltype', method='t-test',n_genes=100)
res=ov.single.pathway_enrichment(atec_adata,pathways_dict=pathway_dict,organism='Human',cutoff=0.05, logfc_threshold=2, pvalue_type='adjust',
                                     group_by='major_celltype',plot=True)





#visualize the enriched GO terms in each urinary aTEC clusters with heatmap colored by AUCell score 
ax=ov.single.pathway_enrichment_plot(res,plot_title='Enrichment',cmap='PuRd',term_num=4,
                                         xticklabels=True,cbar=False,square=True,vmax=10,
                                         yticklabels=True,cbar_kws={'label': '-log10(qvalue)','shrink': 0.5,})

