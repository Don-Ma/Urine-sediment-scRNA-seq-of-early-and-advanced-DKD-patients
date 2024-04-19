import omicverse as ov
from omicverse.utils import mde
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd





#Define aTEC in kidney scRNA-seq data of DKD patients
#import adata
kid_adata=ov.read('kidney_scRNA_DKD.h5ad') #raw data downloaded under GEO accession:GSE131882 by Wilson PC et al Proc Natl Acad Sci U S A., 2019
kid_adata.var_names_make_unique()
kid_adata.obs_names_make_unique()
kid_adata




#visualize expression of PROM1 & DCDC2 
sc.pl.umap(kid_adata,color=['PROM1','DCDC2',])




#calculate the PROM1+/DCDC2+ cells in the PT, LOH, DCT, CD and Mes clusters (presumbly contribute to the urinary aTEC)
gene1 = 'PROM1'
gene2 = 'DCDC2'

kid_adata.obs['aTEC'] = (kid_adata.raw[:,'{}'.format(gene1)].X.todense() > 0) & (kid_adata.raw[:,'{}'.format(gene2)].X.todense() > 0)
kid_adata.obs['aTEC'] = kid_adata.obs['aTEC'].astype(str).astype('category')
kid_adata.obs['aTEC']





#visualize the PROM1+/DCDC2+ cells
ov.utils.plot_embedding_celltype(kid_adata,figsize=(8,4.5),basis='X_umap',
                            celltype_key='aTEC',
                            title='            Cell Counts',
                            celltype_range=(2,4),
                            embedding_range=(4,10),)





#compute the ratio of PROM1+/DCDC2+ cells in the PT, LOH, DCT, CD and Mes clusters
pt_adata=kid_adata[kid_adata.obs['major_celltype'].isin(['PT'])]
loh_adata=kid_adata[kid_adata.obs['major_celltype'].isin(['LOH'])]
dct_adata=kid_adata[kid_adata.obs['major_celltype'].isin(['DCT'])]
ics_adata=kid_adata[kid_adata.obs['major_celltype'].isin(['CD-ICA','CD-ICB'])]
pc_adata=kid_adata[kid_adata.obs['major_celltype'].isin(['CD-PC'])]
mes_adata=kid_adata[kid_adata.obs['major_celltype'].isin(['Mes'])]

ov.utils.plot_cellproportion(adata=pt_adata,celltype_clusters='aTEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))

ov.utils.plot_cellproportion(adata=loh_adata,celltype_clusters='aTEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))

ov.utils.plot_cellproportion(adata=dct_adata,celltype_clusters='aTEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))

ov.utils.plot_cellproportion(adata=ics_adata,celltype_clusters='aTEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))

ov.utils.plot_cellproportion(adata=pc_adata,celltype_clusters='aTEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))

ov.utils.plot_cellproportion(adata=mes_adata,celltype_clusters='aTEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))





#merge urine scRNA-seq data with kidney scRNA-seq data of DKD patients
# obtain cells from the control group of kidney scRNA-seq data
ct_kid_adata=kid_adata[kid_adata.obs['group'].isin(['Control'])]
ct_kid_adata





# obtain cells clusters of interest
adata1=ct_kid_adata[ct_kid_adata.obs['major_celltype'].isin(['CD-PC','CD-ICA','CD-ICB','DCT','DCT/CT','LOH','Mes'])] 
adata1





#import urine scRNA-seq data
adata2=ov.read('urine_aTEC_DKD.h5ad')
adata2.var_names_make_unique()
adata2.obs_names_make_unique()
adata2





# merger data
adata1.obs['type']='Kidney_sn'
adata2.obs['type']='Urine_sc'

adata=sc.concat([adata1,adata2],merge='same')

adata.obs['type'].unique()
adata





#visualize the merged data before perform Harmony to remove batch effect
ov.utils.embedding(adata,
                basis='X_umap',frameon='small',
                color=['type'],show=False)





#perform Harmony to remove batch effect
adata_harmony=ov.single.batch_correction(adata,batch_key='type',
                                        methods='harmony',n_pcs=50) #50
adata	





# embedding & clustering the neighborhood graph again after removing batch effect
sc.pp.neighbors(adata_harmony, n_neighbors=15, n_pcs=50,
               use_rep='X_pca_harmony')	
sc.tl.umap(adata_harmony)
sc.tl.leiden(adata_harmony)





#visualize after removing batch effect
ov.utils.embedding(adata_harmony,
                basis='X_umap',frameon='small',
                color=['type'],show=False)
ov.utils.embedding(adata_harmony,
                basis='X_umap',frameon='small',
                color=['leiden'],show=False)





#Cell type annotation
#prepare marker gene dictionary
res_marker_dict2={
    'PT':['SLC5A12','CUBN','LRP2',],
	'Mes':['PDGFRB','ITGA8','CFH'],
	'CD-ICB':['SLC26A4','ATP6V0D2',], 
	'CD-ICA':['KIT','ATP6V0D2','AQP6'],
	'CD-PC':['SLC8A1','CDH1','GATA3','AQP3',],
	'DCT':['SLC12A3','ATP6V0D2','SLC8A1'],
	'LOH_I':['SLC12A1','UMOD'],
    'LOH_II':['SLC12A2','PROM1','DEFB1',], 
    'aTEC_DCT/LOH_like(urine)':['SLC12A1','UMOD','PROM2','MUC1','SCNN1A','CDH1','KRT7','AQP3','SLC8A1'],
    'aTEC_Mes':['MGP','NID2','THY1','VIM'],
    'aTEC_PT_like(urine)':['HAVCR1','VCAM1','AQP1', 'ANPEP','LRP2','ZEB2','SERPINA1','NCAM1','BMP2','BMPR1A','IGFBP6', 'PTTG1','CENPF', 'TOP2A','MKI67'], 
    
}
# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata,'leiden')
sc.pl.dotplot(adata, res_marker_dict2, 'leiden', 
              dendrogram=True,standard_scale='var')





#annotate urinary aTEC according to the generated dictionary
cluster2annotation2 = {
     '1': 'aTEC_DCT/LOH_like(urine)', 
     '2': 'aTEC_PT_like(urine)', 
     '3': 'CD-PC',
     '4': 'LOH',
     '5': 'DCT',
	 '6': 'aTEC_DCT/LOH_like(urine)',
	 '7': 'CD-PC',
	 '8': 'DCT',
	 '10': 'LOH',
	 '11': 'aTEC_PT_like(urine)',
	 '12': 'CD-ICA',
	 '13': 'DCT',
     '16':'CD-ICA', 
     '17': 'LOH', 
     '18':'Mes', 
     '19':'CD-PC', 
     '21':'CD-ICB',
     '22':'LOH',
     '24':'Mes',
     '26':'CD-PC', 
    '27':'CD-ICA',
    '31':'aTEC_Mes(urine)',

}

ov.single.scanpy_cellanno_from_dict(adata_harmony,anno_dict=cluster2annotation2,
                                       clustertype='leiden',)	
#anno on umap									   
ov.utils.embedding(adata_harmony,
                   basis='X_umap', #X_mde
                   color=['major_celltype'],  
                   legend_loc='right margin',frameon='small', legend_fontoutline=1,
                   palette=ov.utils.palette()[11:],)





#perform trajectory inference by VIA, initiate at aTEC_Mes(urine)
v0 = ov.single.pyVIA(adata=adata_harmony,adata_key='X_pca_harmony',adata_ncomps=80, basis='X_umap',
                         clusters='major_celltype',knn=30,random_seed=42,root_user=["aTEC_Mes(urine)"])
v0.run()





# assign color to cell clusters
fig, ax = plt.subplots(1,1,figsize=(5,4.5))
sc.pl.embedding(
    adata_harmony,
    basis="X_umap",
    color=['major_celltype'],
    frameon=False,
    ncols=1,
    wspace=0.5,
    show=False,
    ax=ax
)





#visualize the inferred trajectory projected onto UMAP
fig,ax1,ax2=v0.plot_trajectory_gams(basis='X_umap',clusters='major_celltype',draw_all_curves=False,figsize=(11,4.5))





#draw a vector stream plot of the more fine-grained directionality of cells along the trajectory projected onto UMAP
fig,ax=v0.plot_stream(basis='X_umap',density_grid=0.4, scatter_size=30, color_scheme='time', linewidth=0.5,
                             min_mass = 1, cutoff_perc = 5, scatter_alpha=0.3, marker_edgewidth=0.1,
                             density_stream = 2, smooth_transition=1, smooth_grid=0.5)





# visualize the probabilistic pathways from root to terminal state 
fig,axs=v0.plot_lineage_probability(figsize=(10,6)) #marker_lineages = [9,17,23,26,29]





# draw heatmaps showing gene dynamics along the pseudotimetime
mark=[ 'SLC12A1','UMOD','ATP6V0D2','PTPRC','GATA3','NPNT','SLC8A1','SCNN1B',]
fig,ax=v0.plot_gene_trend_heatmap(gene_list=mark,figsize=(3,15),
                          marker_lineages=[9,17,23,26,29],)

