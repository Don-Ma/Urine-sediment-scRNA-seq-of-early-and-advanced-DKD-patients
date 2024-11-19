import omicverse as ov
from omicverse.utils import mde
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd


#Data preprocess
# import urine scRNA-seq data from Batch 1
adata1 = sc.read_10x_mtx( 
    '/Batch1/',  # the directory with the `.mtx` file 
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata1.var.index = adata1.var.index.astype(str)
adata1.obs.index = adata1.obs.index.astype(str)

adata1.var_names_make_unique()
adata1.obs_names_make_unique()
adata1=ov.pp.qc(adata1,
              tresh={'mito_perc': 0.15, 'nUMIs': 500, 'detected_genes': 250}) # cells with >15% mtRNA, < 500 UMIs, < 250 detected_genes were filtered out
adata1



# import urine scRNA-seq data from Batch 2
adata2 = sc.read_10x_mtx( 
    '/Batch1/',  # the directory with the `.mtx` file 
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading
adata2.var.index = adata2.var.index.astype(str)
adata2.obs.index = adata2.obs.index.astype(str)

adata2.var_names_make_unique()
adata2.obs_names_make_unique()
adata2=ov.pp.qc(adata2,
              tresh={'mito_perc': 0.15, 'nUMIs': 500, 'detected_genes': 250}) # cells with >15% mtRNA, < 500 UMIs, < 250 detected_genes were filtered out
adata2



#import control data downloaded under GEO accession: GSE157640 by Abedini A et al 2019, JASN
adata3 = sc.read_10x_mtx( 
    '/ctrl/filtered_feature_bc_matrix/',  # the directory with the `.mtx` file 
    var_names='gene_symbols',                # use gene symbols for the variable names (variables-axis index)
    cache=True)                              # write a cache file for faster subsequent reading

adata2.var.index = adata2.var.index.astype(str)
adata2.obs.index = adata2.obs.index.astype(str)
adata3.var_names_make_unique()
adata3.obs_names_make_unique()
adata3=ov.pp.qc(adata3,
              tresh={'mito_perc': 0.15, 'nUMIs': 500, 'detected_genes': 250}) # cells with >15% mtRNA, < 500 UMIs, < 250 detected_genes were filtered out
adata3



#merge different datasets
adata1.obs['batch']='s1mix1' #control
adata2.obs['batch']='s2mix2' #batch1
adata3.obs['batch']='s3mix3' #batch2

adata=sc.concat([adata1,adata2,adata3],merge='same')

adata.obs['batch'].unique()
adata


#ov.utils.store_layers(adata,layers='counts')

#preprocess data
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',
                       n_HVGs=3000,batch_key='batch')
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata



#PCA analysis
ov.pp.scale(adata)
ov.pp.pca(adata,layer='scaled',n_pcs=50)

# Embedding the neighborhood graph
sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50,
               use_rep='scaled|original|X_pca')	
sc.tl.umap(adata)



#clustering the neighborhood graph
sc.tl.leiden(adata)



#visualize before removing batch effect
ov.utils.embedding(adata,
                basis='X_umap',frameon='small',
                color=['batch'],show=False)
ov.utils.embedding(adata,
                basis='X_umap',frameon='small',
                color=['leiden'],show=False)



#Batch effect removal
#perform Harmony to remove batch effect
adata_harmony=ov.single.batch_correction(adata,batch_key='batch',
                                        methods='harmony',n_pcs=50)	
adata



adata_harmony



# embedding & clustering the neighborhood graph again after removing batch effect
sc.pp.neighbors(adata_harmony, n_neighbors=15, n_pcs=50,
               use_rep='X_pca_harmony')	
sc.tl.umap(adata_harmony)
sc.tl.leiden(adata_harmony)


#visualize after removing batch effect
ov.utils.embedding(adata_harmony,
                basis='X_umap',frameon='small',
                color=['batch'],show=False)
ov.utils.embedding(adata_harmony,
                basis='X_umap',frameon='small',
                color=['leiden'],show=False)


#Cell type annotation
#prepare marker gene dictionary
res_marker_dict2={
        'UGEC':['KRT13','PSCA'], 
        'Leuko_kidney':['CD74','APOE','LST1',],
    	'Kidney':[ 'OCIAD2', 'CRYAB','EPCAM'],
	
}
# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata_harmony,'leiden')
sc.pl.dotplot(adata_harmony, res_marker_dict2, 'leiden', 
              dendrogram=True,standard_scale='var')


#annotate cell clusters according to the generated dictionary
cluster2annotation = {
     '0': 'Urogenital_tract(UGEC)', 
     '1': 'Kidney', 
     '2': 'Urogenital_tract(UGEC)',
     '3': 'Leuko_kidney',
     '4': 'Kidney',
     '5': 'Urogenital_tract(UGEC)', 
	 '6': 'Urogenital_tract(UGEC)', 
	 '7': 'Urogenital_tract(UGEC)',
	 '8': 'Kidney',
	 '9': 'Leuko_kidney',
	 '10': 'Kidney',
	 '11': 'Urogenital_tract(UGEC)', 
	 '12': 'Kidney',
	 '13': 'Urogenital_tract(UGEC)', 
	 '14': 'Kidney',
	 '15': 'Kidney',
     '16': 'Leuko_kidney',
     '17': 'Urogenital_tract(UGEC)'
}
ov.single.scanpy_cellanno_from_dict(adata_harmony,anno_dict=cluster2annotation,
                                       clustertype='leiden')	



# visualization with clusters' labels
from matplotlib import patheffects
fig, ax = plt.subplots(figsize=(5,5),)
ov.utils.embedding(adata_harmony,
                   basis='X_umap', #X_mde
                   color=['leiden'],  
                    legend_loc=None, 
                   legend_fontoutline=1,
                   show=False,
                   #frameon='small',
                   ax=ax,
                   palette=ov.utils.palette()[14:],
                  )


ov.utils.gen_mpl_labels(
    adata_harmony,
    'leiden',
    exclude=("None",),  
    basis='X_umap',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 13 ,weight='bold',
                     path_effects=[patheffects.withStroke(linewidth=1, foreground='w')] ),
)



# visualization with celly type labels & colored by group    
from matplotlib import patheffects
fig, ax = plt.subplots(figsize=(5,5),)
ov.utils.embedding(adata_harmony,
                   basis='X_umap', #X_mde
                   color=['group'],  
                    legend_loc=None, 
                   legend_fontoutline=1,
                   show=False,
                   #frameon='small',
                   ax=ax,
                   palette=ov.utils.palette()[16:],
                  )

ov.utils.gen_mpl_labels(
    adata_harmony,
    'major_celltype',
    exclude=("None",),  
    basis='X_umap',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 13 ,weight='bold',
                     path_effects=[patheffects.withStroke(linewidth=1, foreground='w')] ),
)


#save the adata as h5ad file
adata_harmony.write_h5ad('urine_DKD.h5ad',compression='gzip')  



#Cell type variation between uring from early DKD & advanced DKD patients
# obtain cells from early and advanced DKD patients
adata=adata_harmony[adata_harmony.obs['group'].isin(['Early_DKD','Adv_DKD'])]
adata

#prepare marker gene dictionary
res_marker_dict_2nd={
    'UGEC':['KRT13','PSCA','UPK2'], 
    'aTEC':['DCDC2','PROM1','PAX2'], 
    'Podocytes':['NPHS1','NPHS2','PODXL'], 
    'M2_Macro':['C1QA','C1QB','C1QC'], 
    'MDC':['LST1','CLEC10A','CD1C'],
}

# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata,'major_celltype')
sc.pl.dotplot(adata, res_marker_dict_2nd, 'major_celltype', 
              dendrogram=True,standard_scale='var',cmap='Greens')



#annotate cell clusters according to the generated dictionary
cluster2annotation = {
     '0': 'UGEC', 
     '1': 'aTEC_LOH_like',
     '2': 'aTEC_PT_like',
     '3': 'M2_Macro',
     '4': 'UGEC',
     '5': 'UGEC',
	 '6': 'aTEC_DCT_like', 
	 '7': 'MDC',
	 '8': 'aTEC_Mes',
	 '9': 'UGEC',
	 '10': 'UGEC',
	 '11': 'aTEC_PT_prolf',
	 '12': 'Podocyte', 
	 '13': 'aTEC_PT_like',
	 '14': 'M2_Macro',
	 '15': 'UGEC',
}

ov.single.scanpy_cellanno_from_dict(adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')									   
#visualization with cell type labels 
ov.utils.embedding(adata,
                   basis='X_umap', #X_mde
                   color=['major_celltype'],  
                   legend_loc='right margin',frameon='small', legend_fontoutline=1, 
                   palette=ov.utils.palette()[13:],)



#compute cell type ratio
ov.utils.plot_cellproportion(adata=adata,celltype_clusters='major_celltype',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))



#extract the TEC cluster
gene1 = 'OCIAD2'
gene2 = 'CRYAB'
gene3 = 'EPCAM'

adata.obs['TEC'] = (adata.raw[:,'{}'.format(gene1)].X.todense() > 0) & (adata.raw[:,'{}'.format(gene2)].X.todense() > 0) & (adata.raw[:,'{}'.format(gene3)].X.todense() > 0)



#compute cell type ratio					   
ov.utils.plot_cellproportion(adata=adata,celltype_clusters='TEC',
                    visual_clusters='group',
                    visual_name='group',figsize=(2,4))
