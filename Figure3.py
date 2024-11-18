import omicverse as ov
from omicverse.utils import mde
import scanpy as sc
import matplotlib.pyplot as plt
import os
import pandas as pd






#import merged adata containing kidney and urine scRNA-seq data
kid_adata=ov.read('kidney&urinePT_scRNA_DKD.h5ad') # processed kidney scRNA-seqdata is downloaded under GEO accession:GSE211785 
kid_adata.var_names_make_unique()
kid_adata.obs_names_make_unique()
kid_adata



#perform trajectory inference 
Traj=ov.single.TrajInfer(kid_adata,basis='X_umap',groupby='Cluster_Idents',
                         use_rep='scaled|original|X_pca',n_comps=50)
Traj.set_origin_cells('iPT',)
Traj.set_terminal_cells(['urine_aTEC_PT_like','urine_aTEC_PT_prolf','urine_aTEC_Mes_like',"Fibroblast_1","Fibroblast_2",'MyoFib/VSMC'])
Traj.inference(method='palantir',num_waypoints=500)



#visualize on UMAP
Traj.palantir_plot_pseudotime(embedding_basis='X_umap',cmap='RdBu_r',s=3)



#visualize the selection of cells
Traj.palantir_cal_branch(eps=0,q=.5)


#import merged adata containing kidney and urine scRNA-seq data
kid_adata=ov.read('kidney&urineDT_scRNA_DKD.h5ad') # processed kidney scRNA-seqdata is downloaded under GEO accession:GSE211785 
kid_adata.var_names_make_unique()
kid_adata.obs_names_make_unique()
kid_adata



#perform trajectory inference 
Traj=ov.single.TrajInfer(kid_adata,basis='X_umap',groupby='Cluster_Idents',
                         use_rep='scaled|original|X_pca',n_comps=50)
Traj.set_origin_cells('urine_aTEC_DCT/CT/LOH_like',)
Traj.set_terminal_cells(['C_TAL','M_TAL','Des-Thin_Limb','Ascending_Thin_LOH','DCT1','DCT2','CNT','IC_A','IC_B','PC',])
Traj.inference(method='palantir',num_waypoints=500)



#visualize on UMAP
Traj.palantir_plot_pseudotime(embedding_basis='X_umap',cmap='RdBu_r',s=3)



#visualize the selection of cells
Traj.palantir_cal_branch(eps=0,q=.005)





