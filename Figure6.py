import omicverse as ov
from omicverse.utils import mde
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd
from matplotlib import patheffects
import matplotlib.pyplot as plt





#Examine aTEC in urine of COVID-19 patients with AKI
#import adata
adata=ov.read('urine_covid19_AKI.h5ad') #raw data downloaded under GEO accession: GSE199321 by Klocke J, Kim SJ et al Kidney Int, 2022
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata





#prepare marker gene dictionary
res_marker_dict={
    '0-UGEC':['KRT13','PSCA','UPK2'],
    '1-Immu':['S100A4','LST1','CD74'],
    '2-aTEC':['DCDC2','PAX2','PROM1'], 
    '3-Podo':['NPHS1','NPHS2','PODXL',], 
	}

# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata,'major_celltype')
sc.pl.dotplot(adata, res_marker_dict, 'major_celltype', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# calculated the expression of marker genes in each cluster and the fraction	
#sc.tl.dendrogram(adata,'leiden')
#sc.pl.dotplot(adata, res_marker_dict, 'leiden', 
#              dendrogram=True,standard_scale='var')





#annotate urinary aTEC according to the generated dictionary
cluster2annotation = {
     '0': 'Immu', 
     '1': 'Immu',
     '2': 'Immu',
     '3': 'UGEC',
     '4': 'UGEC',
     '5': 'aTEC',
	 '6': 'Immu', 
	 '7': 'aTEC',
	 '8': 'Immu',
	 '9': 'Immu',
	 '10': 'Immu',
	 '11': 'UGEC',
	 '12': 'Immu', 
	 '13': 'aTEC',
	 '14': 'UGEC',
	 '15': 'Immu',
     '16':'UGEC',
     '17':'Immu',
     '18':'Immu',
     '19':'aTEC',
     '20':'Podo',
     '21':'Immu',
     '22':'aTEC',
     '23':'UGEC',
    
}

ov.single.scanpy_cellanno_from_dict(adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')	
    





# visualization
fig, ax = plt.subplots(figsize=(5,5))

ov.utils.embedding(adata,
                  basis='X_umap',
                  color=['major_celltype'],
                   show=False, legend_loc=None, add_outline=False, 
                   legend_fontoutline=2,ax=ax,#frameon='small',#
                   palette=['#941456','#FCBC10','#01A0A7','#EAEFC5',]#'#E0A7C8',
                   #palette=ov.utils.palette()[15:]
                   #palette=ov.utils.green_color,
                 )
ov.utils.gen_mpl_labels(
    adata,
    'major_celltype',
    exclude=("None",),  
    basis='X_umap',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 15 ,weight='bold',
    path_effects=[patheffects.withStroke(linewidth=2, foreground='w')] ),
)

ov.utils.plot_ConvexHull(adata,
                basis='X_umap',
                cluster_key='major_celltype',
                color='grey',
                hull_cluster='aTEC',
                ax=ax)





#prepare marker gene dictionary for different state
tec=adata[adata.obs['major_celltype'].isin(['aTEC'])]
res_marker_dict2={
    'Progenitor':['SOX4','SOX9','PAX2','HES1',],
    'Fibrosis':['MGP', 'BGN', 'FN1','ACTA2',],
    'Degenerative':['APOE','CLU','S100A6','FTL',],
    'Proliferative':['PTTG1','TOP2A','MKI67','CENPF'], 
}


# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(tec,'manual_annotations')
sc.pl.dotplot(tec, res_marker_dict2, 'manual_annotations', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# visualization
ov.pl.embedding(tec,
                basis='X_umap',
                color=['manual_annotations'],palette=['#EAEFC5','#01A0A7','#75C8CC','#1F577B'],
                frameon='small')


# compute ratio
ov.utils.plot_cellproportion(adata=tec,celltype_clusters='manual_annotations',
                    visual_clusters='group',           
                    visual_name='group',figsize=(2,6))





#Examine aTEC in urine of cardiac surgery patients with AKI
#import adata
adata=ov.read('urine_cardi_AKI.h5ad') #raw data downloaded under GEO accession: GSE199321 by Klocke J, Kim SJ et al Kidney Int, 2022
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata





#prepare marker gene dictionary
res_marker_dict={
    '0-UGEC':['KRT13','PSCA','UPK2'],
    '1-Immu':['S100A4','LZY','CD74'],
    '2-aTEC':['DCDC2','PAX2','PROM1'], 
    '3-PT':['MT1G','MT1H','MME'], 
	}

# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata,'major_celltype')
sc.pl.dotplot(adata, res_marker_dict, 'major_celltype', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# calculated the expression of marker genes in each cluster and the fraction	
#sc.tl.dendrogram(adata,'leiden')
#sc.pl.dotplot(adata, res_marker_dict, 'leiden', 
#              dendrogram=True,standard_scale='var')




#annotate urinary aTEC according to the generated dictionary
cluster2annotation = {
     '0': 'UGEC', 
     '1': 'aTEC',
     '2': 'UGEC',
     '3': 'aTEC',
     '4': 'Immu',
     '5': 'UGEC',
	 '6': 'UGEC', 
	 '7': 'Immu',
	 '8': 'UGEC',
	 '9': 'aTEC',
	 '10': 'UGEC',
	 '11': 'aTEC',
	 '12': 'UGEC', 
	 '13': 'Immu',
	 '14': 'aTEC',
	 '15': 'aTEC',
     '16':'aTEC',
     '17':'Immu',
     '18':'Immu',
     '19':'PT',
     '20':'aTEC',
     '21':'Immu',
     '22':'UGEC',
     '23':'UGEC',
     '24':'UGEC',
     '25':'UGEC',
     '26':'aTEC', 
     '27':'aTEC',
    
}

ov.single.scanpy_cellanno_from_dict(adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')	
    





# visualization
fig, ax = plt.subplots(figsize=(5,5))

ov.utils.embedding(adata,
                  basis='X_umap',
                  color=['major_celltype'],
                   show=False, legend_loc=None, add_outline=False, 
                   legend_fontoutline=2,ax=ax,
                    palette=['#941456','#FCBC10','#01A0A7','#EAEFC5',]
                 )
ov.utils.gen_mpl_labels(
    adata,
    'major_celltype',
    exclude=("None",),  
    basis='X_umap',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 15 ,weight='bold',
    path_effects=[patheffects.withStroke(linewidth=2, foreground='w')] ),
)

ov.utils.plot_ConvexHull(adata,
                basis='X_umap',
                cluster_key='major_celltype',
                color='grey',
                hull_cluster='aTEC',
                ax=ax)





#prepare marker gene dictionary for different state
tec=adata[adata.obs['major_celltype'].isin(['aTEC'])]
res_marker_dict2={
    'Progenitor':['SOX4','SOX9','PAX2','HES1'],
    'Fibrosis':['MGP', 'BGN', 'FN1',],
    'Degenerative':['APOE','CLU','S100A6','FTL'],
    'Proliferative':['PTTG1','TOP2A','MKI67','CENPF'], 
}


# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(tec,'manual_annotations')
sc.pl.dotplot(tec, res_marker_dict2, 'manual_annotations', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# visualization
ov.pl.embedding(tec,
                basis='X_umap',
                color=['manual_annotations'],palette=['#EAEFC5','#01A0A7','#75C8CC','#1F577B'],
                frameon='small')


# compute ratio
ov.utils.plot_cellproportion(adata=tec,celltype_clusters='manual_annotations',
                    visual_clusters='group',           
                    visual_name='group',figsize=(2,6))




#Examine aTEC in urine of DKD patients
#import adata
adata=ov.read('urine_DKD.h5ad') #raw data downloaded under GEO accession: GSE157640 by Abedini A et al JASN, 2019
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata





#prepare marker gene dictionary
res_marker_dict={
    '0-UGEC':['KRT13','PSCA','UPK2'],
    '1-Immu':['S100A4','LST1','CD74'],
    '2-aTEC':['DCDC2','PAX2','PROM1'], 
    '3-PT':['MT1G','MT1H','MME'], 
	}

# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata,'major_celltype')
sc.pl.dotplot(adata, res_marker_dict, 'major_celltype', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# calculated the expression of marker genes in each cluster and the fraction	
#sc.tl.dendrogram(adata,'leiden')
#sc.pl.dotplot(adata, res_marker_dict, 'leiden', 
#              dendrogram=True,standard_scale='var')




#annotate urinary aTEC according to the generated dictionary
cluster2annotation = {
     '0': 'UGEC', 
     '1': 'UGEC',
     '2': 'UGEC',
     '3': 'UGEC',
     '4': 'UGEC',
     '5': 'UGEC',
	 '6': 'UGEC', 
	 '7': 'UGEC',
	 '8': 'UGEC',
	 '9': 'UGEC',
	 '10': 'UGEC',
	 '11': 'aTEC',
	 '12': 'UGEC', 
	 '13': 'UGEC',
	 '14': 'Immu',
	 '15': 'Immu',
     '16':'UGEC',
     '17':'aTEC',
     '18':'PT',
     '19':'UGEC',
     '20':'Immu',
     '21':'UGEC',
     '22':'aTEC',
     '23':'Immu',
    
}

ov.single.scanpy_cellanno_from_dict(adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')	
    





# visualization
fig, ax = plt.subplots(figsize=(5,5))

ov.utils.embedding(adata,
                  basis='X_umap',
                  color=['major_celltype'],
                   show=False, legend_loc=None, add_outline=False, 
                   legend_fontoutline=2,ax=ax,
                    palette=['#941456','#FCBC10','#01A0A7','#EAEFC5',]
                 )
ov.utils.gen_mpl_labels(
    adata,
    'major_celltype',
    exclude=("None",),  
    basis='X_umap',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 15 ,weight='bold',
    path_effects=[patheffects.withStroke(linewidth=2, foreground='w')] ),
)

ov.utils.plot_ConvexHull(adata,
                basis='X_umap',
                cluster_key='major_celltype',
                color='grey',
                hull_cluster='aTEC',
                ax=ax)





#prepare marker gene dictionary for different state
tec=adata[adata.obs['major_celltype'].isin(['aTEC'])]
res_marker_dict2={
    'Progenitor':['SOX4','SOX9','PAX2','HES1'],
    'Fibrosis':['MGP', 'BGN', 'FN1',],
    'Degenerative':['APOE','CLU','S100A6','FTL'],
    'Proliferative':['PTTG1','TOP2A','MKI67','CENPF'], 
}


# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(tec,'manual_annotations')
sc.pl.dotplot(tec, res_marker_dict2, 'manual_annotations', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# visualization
ov.pl.embedding(tec,
                basis='X_umap',
                color=['manual_annotations'],palette=['#EAEFC5','#01A0A7','#75C8CC','#1F577B'],
                frameon='small')


# compute ratio
ov.utils.plot_cellproportion(adata=tec,celltype_clusters='manual_annotations',
                    visual_clusters='group',           
                    visual_name='group',figsize=(2,6))




#Examine aTEC in urine of FSGS patients
#import adata
adata=ov.read('urine_FSGS.h5ad') #raw data downloaded under GEO accession: GSE176465 by Latt KZ et al Kidney Int Rep, 2022
adata.var_names_make_unique()
adata.obs_names_make_unique()
adata



#prepare marker gene dictionary
res_marker_dict={
    '0-UGEC':['KRT13','PSCA','UPK2'],
    '1-Immu':['S100A4','LST1','CD74'],
    '2-aTEC':['DCDC2','PAX2','PROM1'], 
    '3-PT':['MT1G','MT1H','MME'], 
	}

# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(adata,'major_celltype')
sc.pl.dotplot(adata, res_marker_dict, 'major_celltype', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# calculated the expression of marker genes in each cluster and the fraction	
#sc.tl.dendrogram(adata,'leiden')
#sc.pl.dotplot(adata, res_marker_dict, 'leiden', 
#              dendrogram=True,standard_scale='var')




#annotate urinary aTEC according to the generated dictionary
cluster2annotation = {
     '0': 'UGEC', 
     '1': 'UGEC',
     '2': 'UGEC',
     '3': 'UGEC',
     '4': 'UGEC',
     '5': 'UGEC',
	 '6': 'UGEC', 
	 '7': 'UGEC',
	 '8': 'UGEC',
	 '9': 'UGEC',
	 '10': 'UGEC',
	 '11': 'aTEC',
	 '12': 'UGEC', 
	 '13': 'UGEC',
	 '14': 'Immu',
	 '15': 'Immu',
     '16':'UGEC',
     '17':'aTEC',
     '18':'PT',
     '19':'UGEC',
     '20':'Immu',
     '21':'UGEC',
     '22':'aTEC',
     '23':'Immu',
    
}

ov.single.scanpy_cellanno_from_dict(adata,anno_dict=cluster2annotation,
                                       clustertype='leiden')	
    





# visualization
fig, ax = plt.subplots(figsize=(5,5))

ov.utils.embedding(adata,
                  basis='X_umap',
                  color=['major_celltype'],
                   show=False, legend_loc=None, add_outline=False, 
                   legend_fontoutline=2,ax=ax,
                    palette=['#941456','#FCBC10','#01A0A7','#EAEFC5',]
                 )
ov.utils.gen_mpl_labels(
    adata,
    'major_celltype',
    exclude=("None",),  
    basis='X_umap',
    ax=ax,
    adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
    text_kwargs=dict(fontsize= 15 ,weight='bold',
    path_effects=[patheffects.withStroke(linewidth=2, foreground='w')] ),
)

ov.utils.plot_ConvexHull(adata,
                basis='X_umap',
                cluster_key='major_celltype',
                color='grey',
                hull_cluster='aTEC',
                ax=ax)




#prepare marker gene dictionary for different state
tec=adata[adata.obs['major_celltype'].isin(['aTEC'])]
res_marker_dict2={
    'Progenitor':['SOX4','SOX9','PAX2','HES1'],
    'Fibrosis':['MGP', 'BGN', 'FN1','ACTA2',],
    'Degenerative':['APOE','CLU','S100A6','FTL'],
    'Proliferative':['PTTG1','TOP2A','MKI67','CENPF'], 
}


# calculated the expression of marker genes in each cluster and the fraction	
sc.tl.dendrogram(tec,'manual_annotations')
sc.pl.dotplot(tec, res_marker_dict2, 'manual_annotations', 
              dendrogram=False,standard_scale='var',cmap='Greens')

# visualization
ov.pl.embedding(tec,
                basis='X_umap',
                color=['manual_annotations'],palette=['#EAEFC5','#01A0A7','#75C8CC','#1F577B'],
                frameon='small')


# compute ratio
ov.utils.plot_cellproportion(adata=tec,celltype_clusters='manual_annotations',
                    visual_clusters='group',           
                    visual_name='group',figsize=(2,6))
