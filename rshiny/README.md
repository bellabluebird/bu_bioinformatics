# rnaseq dashboard

available online here: https://bellabluebird.shinyapps.io/rnaseq_dashboard/

interactive shiny app for exploring huntington's disease RNA-seq data from GSE64810. built for genomics research and data exploration. designed specifically for GSE64810 but can be adapted for other RNA-seq datasets.

this dashboard lets you analyze differential gene expression data comparing huntington's disease patients vs healthy controls in brain tissue. you can upload your data files and explore them through multiple analysis modules.

## data files needed

- **GSE64810_series_matrix.txt** - sample metadata from GEO
- **GSE64810_mlhd_DESeq2_diffexp_DESeq2_outlier_trimmed_adjust.txt** - differential expression results  
- **GSE64810_mlhd_DESeq2_norm_counts_adjust.txt** - normalized expression counts
- **pathway_gene_sets.gmt** - gene sets for pathway analysis (from MSigDB)

## analysis information + modules

### sample metadata
- clinical characteristics and quality metrics
- age, RIN scores, sequencing depth analysis
- HD-specific data: CAG repeats, disease duration, vonsattel grades
- statistical comparisons between groups

### volcano plot
- interactive differential expression visualization
- adjustable significance thresholds
- gene selection and filtering
- downloadable results tables

### expression analysis
- qc diagnostic plots
- PCA analysis for sample clustering
- expression heatmaps with hierarchical clustering
- variance analysis and gene rankings

### GSEA analysis
- gene set enrichment analysis using fgsea
- multiple ranking methods for DESeq2 data
- pathway visualization and statistics
- leading edge gene identification
  
### features
- handles large datasets efficiently (up to 500MB uploads)
- interactive plots with plotly
- professional styling with custom CSS
- modular architecture for easy maintenance
- comprehensive error handling and validation
- downloadable results in multiple formats

### dataset info
- **study**: huntington's disease vs control brain samples
- **tissue**: prefrontal cortex (brodmann area 9)
- **samples**: 69 total (20 HD patients + 49 controls)
- **platform**: illumina RNA-sequencing
- **genes**: ~28,000 in final dataset
