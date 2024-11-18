import omicverse as ov
#print(f"omicverse version: {ov.__version__}")
import scanpy as sc
#print(f"scanpy version: {sc.__version__}")
ov.utils.ov_plot_set()





#import urine scRNA-seq data of DKD patients generated from this study
adata_sc=ov.read('urine_DKD.h5ad')
adata_sc.var_names_make_unique()
adata_sc.obs_names_make_unique()
adata_sc





#import kidney cortex spRNA-seq data of DKD patients 
adata=sc.read_visium('/path_to_file/',genome=None, count_file='filtered_feature_bc_matrix.h5', library_id=None, load_images=True, source_image_path=None) #processed data downloaded under GEO accession: GSE183456 by Lake BB et al Nature, 2023
adata





#preprocess the spRNA-seq data
adata.obs['sample'] = list(adata.uns['spatial'].keys())[0]
adata.var_names_make_unique()

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:,adata.var['total_counts']>100]
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=3000,target_sum=1e4)
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata_sp=adata.copy()
sc.pl.spatial(adata, color=['pct_counts_in_top_500_genes'])





# apply Tangram to map the scRNA-seq data onto spRNA-seq data
tg=ov.space.Tangram(adata_sc,adata_sp,clusters='major_celltype')
tg.train(mode="clusters",num_epochs=500,device="cpu")





#obtain loacation for each cell type
adata_plot=tg.cell2location()
adata_plot.obs.columns
annotation_list=['aTEC_PT_like','aTEC_Mes_like','aTEC_PT_prolf',]

sc.pl.spatial(adata_plot, cmap='afmhot',
                  color=annotation_list,
                  ncols=4, size=1.3,
                  img_key='hires',
                 )





#obtain spatial expression profile for selective genes
annotation_list=['SLC4A4','LRP2','NPHS2']

sc.pl.spatial(adata_plot, cmap='afmhot',
                  # show first 8 cell types
                  color=annotation_list,
                  ncols=6, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.02'
                 )


				 
#import kidney medulla spRNA-seq data of DKD patients 
adata=sc.read_visium('/path_to_file/',genome=None, count_file='filtered_feature_bc_matrix.h5', library_id=None, load_images=True, source_image_path=None) #processed data downloaded under GEO accession: GSE183456 by Lake BB et al Nature, 2023
adata





#preprocess the spRNA-seq data
adata.obs['sample'] = list(adata.uns['spatial'].keys())[0]
adata.var_names_make_unique()

sc.pp.calculate_qc_metrics(adata, inplace=True)
adata = adata[:,adata.var['total_counts']>100]
adata=ov.pp.preprocess(adata,mode='shiftlog|pearson',n_HVGs=3000,target_sum=1e4)
adata.raw = adata
adata = adata[:, adata.var.highly_variable_features]
adata_sp=adata.copy()
sc.pl.spatial(adata, color=['pct_counts_in_top_500_genes'])





# apply Tangram to map the scRNA-seq data onto spRNA-seq data
tg=ov.space.Tangram(adata_sc,adata_sp,clusters='major_celltype')
tg.train(mode="clusters",num_epochs=500,device="cpu")





#obtain loacation for each cell type
adata_plot=tg.cell2location()
adata_plot.obs.columns
annotation_list=['Progenitor_aTEC',]

sc.pl.spatial(adata_plot, cmap='afmhot',
                  color=annotation_list,
                  ncols=4, size=1.3,
                  img_key='hires',
                 )





#obtain spatial expression profile for selective genes
annotation_list=['SLC14A1','UMOD','SOX4','PAX2']

sc.pl.spatial(adata_plot, cmap='afmhot',
                  # show first 8 cell types
                  color=annotation_list,
                  ncols=6, size=1.3,
                  img_key='hires',
                  # limit color scale at 99.2% quantile of cell abundance
                  vmin=0, vmax='p99.02'
                 )

				 
