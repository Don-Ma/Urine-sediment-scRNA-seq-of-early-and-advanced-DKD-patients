import cell2cell as c2c
import scanpy as sc
import pandas as pd

get_ipython().run_line_magic('matplotlib', 'inline')
import warnings
warnings.filterwarnings('ignore')


#Annotation of urinary aTEC using marker genes of adaptive states
#import urine scRNA-seq data of DKD patients generated from this study
atec_adata=ov.read('urine_all_celltype_DKD.h5ad')
atec_atec_adata.var_names_make_unique()
atec_adata.obs_names_make_unique()
atec_adata





# download the Ligand-Receptor Pairs (LR)
lr_pairs = pd.read_csv('https://raw.githubusercontent.com/LewisLabUCSD/Ligand-Receptor-Pairs/master/Human/Human-2020-Jin-LR-pairs.csv') #provided with CellChat by Jin et al. 2021, Nature Communications
lr_pairs = lr_pairs.astype(str)





#prepare metadata for the single cells
meta = rnaseq.obs.copy()





#integrates the scRNA-seq and LR datasets
interactions = c2c.analysis.SingleCellInteractions(rnaseq_data=rnaseq.to_df().T,
                                                   ppi_data=lr_pairs,
                                                   metadata=meta,
                                                   interaction_columns=('ligand_symbol', 'receptor_symbol'),
                                                   communication_score='expression_thresholding',
                                                   expression_threshold=0.1, # values after aggregation
                                                   cci_score='bray_curtis',
                                                   cci_type='undirected',
                                                   aggregation_method='nn_cell_fraction',
                                                   barcode_col='index',
                                                   celltype_col='major_celltype',
                                                   complex_sep='&',
                                                   verbose=False)





# compute communication scores for each LR pair
interactions.compute_pairwise_communication_scores()





# compute CCI scores for each pair of cells
interactions.compute_pairwise_cci_scores()





# perform permutation analysis
cci_pvals = interactions.permute_cell_labels(evaluation='interactions', 
                                             permutations=1000, #1000 or above
                                             fdr_correction=True,
                                             verbose=True)





# visualizations
# generate a metadata for the cell types
meta.major_celltype.unique()
group_meta = pd.DataFrame(columns=['Celltype', 'Group'])
group_meta['Celltype'] = meta.major_celltype.unique().tolist()
group_meta['Group'] = ['Leukocyte', 'Epithelial','Epithelial','Leukocyte','Fibrosis','Mesenchymal','Fibrosis']





# generate colors for the groups in the metadata
colors = c2c.plotting.get_colors_from_labels(labels=group_meta['Group'].unique().tolist(),
                                             cmap='coolwarm'
                                            )





# visualize communication scores for each LR pair and each cell pair
interaction_clustermap = c2c.plotting.clustermap_ccc(interactions,
                                                     metric='jaccard',
                                                     method='complete',
                                                     metadata=group_meta,
                                                     sample_col='Celltype',
                                                     group_col='Group',
                                                     colors=colors,
                                                     row_fontsize=14,
                                                     title='Active ligand-receptor pairs for interacting cells',
                                                     filename=None,
                                                     cell_labels=('SENDER-CELLS', 'RECEIVER-CELLS'),
                                                     **{'figsize' : (10,9),
                                                       }
                                                     )

# Add a legend to know the groups of the sender and receiver cells:
l1 = c2c.plotting.generate_legend(color_dict=colors,
                                  loc='center left',
                                  bbox_to_anchor=(20, -2), # Indicated where to include it
                                  ncol=1, fancybox=True,
                                  shadow=True,
                                  title='Groups',
                                  fontsize=14,
                                 )





# visualized with Circos plot
#communication between aTEC_PT & aTEC_Mes 
sender_cells = [ 'aTEC_PT_like',  'aTEC_PT_prolf']
receiver_cells = ['aTEC_Mes',] 
ligands = ['LAMB3','CD99','SPP1','COL6A1','MDK','CDH1','COL4A2','COL6A2','CEACAM1','COL1A1','TNC','LAMB1','ANGPTL4','THBS2']#refer to the clustermap result
receptors = ['CD99','CD44','ITGAV&ITGB8','SDC2','SDC4','SDC1','ITGAV&ITGB5','ITGA3&ITGB1','CEACAM5','ITGAV&ITGB3','ITGAV&ITGB6']#refer to the clustermap result

c2c.plotting.circos_plot(interaction_space=interactions,
                         sender_cells=sender_cells,
                         receiver_cells=receiver_cells,
                         ligands=ligands,
                         receptors=receptors,
                         excluded_score=0,
                         metadata=group_meta,
                         sample_col='Celltype',
                         group_col='Group',
                         colors=colors,
                         fontsize=15,
                        )





#communication between aTEC_DCT/aTEC_LOH & leuko
sender_cells =  ['aTEC_DCT_like', 'aTEC_LOH_like']
receiver_cells =['MDC','M2_Macro']
ligands = ['LAMB3','CD99','SPP1','COL6A1','MDK','CDH1','COL4A2','COL6A2','CEACAM1','COL1A1','TNC','LAMB1','ANGPTL4','THBS2']#refer to the clustermap result
receptors = ['CD99','CD44','ITGAV&ITGB8','SDC2','SDC4','SDC1','ITGAV&ITGB5','ITGA3&ITGB1','CEACAM5','ITGAV&ITGB3','ITGAV&ITGB6']#refer to the clustermap result

c2c.plotting.circos_plot(interaction_space=interactions,
                         sender_cells=sender_cells,
                         receiver_cells=receiver_cells,
                         ligands=ligands,
                         receptors=receptors,
                         excluded_score=0,
                         metadata=group_meta,
                         sample_col='Celltype',
                         group_col='Group',
                         colors=colors,
                         fontsize=15,
                        )

