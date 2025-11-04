# BulkR
<img src="https://github.com/vgrozd/BulkR/blob/main/BulkR_Tilev2.png" width="100">

## Tools for single-cell and pseudobulk analysis  
BulkR is a collection of tools and synthetic data examples for the analyis of single-cell and bulk sequencing data. 


## Installation 

BulkR is implemented in R and can be installed in R with: 
```
devtools::install_github("vgrozd/BulkR") 
```

## Usage 

For detailed usade, see the BulkR vignette. 

* #### Combine DGE analysis results

With this function, differential gene expression analysis results can be combined into a single object. 
Acceptet inputs are DESeq2 results as a DESeq2 results object or a data.frame and MAST results with standard 
column names. 

```combine_results()``` formats the inputs to the same format (column names) and returns a list of formatted 
data.frames or a single dataframe with concatenated columns. 


#### Extract significant results 

With this funcition, names of significant features can be extracted. 
Alternatively, a data.frame with filtered for results with ```fdr<threshold``` will be returned. 

#### Hierarchical gene exression 

This function returns a BulkedExpression object with pseudo-bulk gene expression across hierarchical or unrelated 
sample factors. 

### Hierarchical differential expression analysis 

This function automates differential gene expression analysis across multiple levels. 

One of the two default methods, "DESeq2"  or "MAST" can be specified. 
Alternatively, a user-defined analysis function can be accepted if it returns an accepted results format. 

### Systematic SVA 

This function performs SVA analysis with increasing number of SVs and stores the data in a BulkedExpression object. 
The number of SVs can be specified. 



