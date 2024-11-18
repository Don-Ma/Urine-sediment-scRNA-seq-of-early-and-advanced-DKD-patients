import omicverse as ov
from omicverse.utils import mde
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd

#Annotation of urinary TEC using marker genes of adaptive states
# import urine scRNA-seq data (TEC) of DKD patients generated from this study
atec_adata=ov.read('urine_aTEC_DKD.h5ad')
atec_atec_adata.var_names_make_unique()
atec_adata.obs_names_make_unique()
atec_adata



#prepare marker gene dictionary
res_marker_dict={
    'aTAL':['PROM1','DCDC2'],
    'aPT':['HAVCR1','VCAM1'],
    'norPT':['AQP1', 'ANPEP','LRP2',], 
    'Mes':['MGP','NID2','THY1',],
    'norLOH':['SLC12A1','UMOD',"TACSTD2",],
    'norDCT/CT':["FGF13",'GATA3','AQP3'],
}
	
# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(atec_adata,'leiden')
sc.pl.dotplot(atec_adata, res_marker_dict, 'leiden', 
              dendrogram=True,standard_scale='var',cmap='coolwarm')

#prepare marker gene dictionary
res_marker_dict_adj={
    'Progenitor':['SOX4','SOX9','PAX2','HES1',],
    'Fibrosis':['MGP', 'BGN', 'FN1','ACTA2',],
    'EMT':['ZEB2','SERPINA1','NCAM1', 'BMP2',],
    'Degenerative':['APOE','CLU','S100A6','FTL',],
    'Proliferative':['PTTG1','TOP2A','MKI67','CENPF'], 
}

# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(atec_adata,'leiden')
sc.pl.dotplot(atec_adata, res_marker_dict_adj, 'leiden', 
              dendrogram=True,standard_scale='var',cmap='coolwarm')



#annotate urinary aTEC according to the generated dictionary
cluster2annotation = {
         '1': 'aTEC_LOH_like',
         '2': 'aTEC_PT_like',
	 '6': 'aTEC_DCT/CT_like', 
	 '8': 'aTEC_Mes_like',
	 '11': 'aTEC_PT_prolf',
	 '13': 'aTEC_PT_like',
}

ov.single.scanpy_cellanno_from_dict(atec_adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')									   
ov.utils.embedding(atec_adata,
                   basis='X_umap', #X_mde
                   color=['major_celltype'],  
                   legend_loc='right margin',frameon='small', legend_fontoutline=1, 
                   palette=['#E069A6','#E6B7D2','#941456','#FFDEAD','#D2B48C',]) 



#compute cell type ratio
ov.utils.plot_cellproportion(adata=atec_adata,celltype_clusters='major_celltype',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))




#GO analysis
#prepare GO library
pathway_dict=ov.utils.geneset_prepare('GO_2021.txt',organism='Human') # Downloaded from https://maayanlab.cloud/Enrichr/#libraries





# calculate and visualize all GO terms enriched in early DKD & advanced DKD 
#ely=atec_adata[atec_adata.obs['group'].isin(['Early_DKD'])]
adv=atec_adata[atec_adata.obs['group'].isin(['Adv_DKD'])]
adv.uns['log1p']['base']=None
sc.tl.rank_genes_groups(adv, 'major_celltype', method='t-test',n_genes=200)
res=ov.single.pathway_enrichment(adv,pathways_dict=pathway_dict,organism='Human',cutoff=0.05, logfc_threshold=2, pvalue_type='adjust',
                                     group_by='major_celltype',plot=True)




#compute AUCelll score based on all GO terms
adata_aucs=ov.single.pathway_aucell_enrichment(atec_adata,
                                                  pathways_dict=pathway_dict,
                                                  num_workers=8)
adata_aucs.obs=atec_adata[adata_aucs.obs.index].obs
adata_aucs.obsm=atec_adata[adata_aucs.obs.index].obsm
adata_aucs.obsp=atec_adata[adata_aucs.obs.index].obsp
adata_aucs




#visualize the enriched GO terms in each urinary aTEC clusters with heatmap colored by AUCell score 
ax=ov.single.pathway_enrichment_plot(res,plot_title='Enrichment',cmap='OrRd',
                                         xticklabels=True,cbar=False,square=True,vmax=10,
                                         yticklabels=True,cbar_kws={'label': '-log10(qvalue)','shrink': 0.5,})





