
# DAAMS
Differential abundance analysis for shotgun metagenomics. This repo is inspired by https://github.com/nearinj/Comparison_of_DA_microbiome_methods. Changes are made here to include multiple variables for some methods and also more suitable for shotgun metagenomics data. 


For quick start take a look at [analysis.Rmd](analysis.Rmd)


## Summary

This is a simple wrapper that invokes several methods for differential abundance analysis.

For a single fixed variable the following methods are run:

1. DeSeq
2. Maaslin2
3. EdgeR
4. MetagenomeSeq
5. ANCOM2
6. ANCOM-BC

For two fixed variables the following methods are run:
1. DeSeq
2. Maaslin2
3. ANCOM-BC

For fixed variables and a random variable (like subject in repeated measures) the following methods are run:

1. DeSeq
2. Maaslin2
3. ANCOM-BC

## Input

The main function to use is `wrapper_daa()` and following inputs are required by this function:

1. A comma separated species table is required. See [example](example/merged_species_table.csv)

2. A metadata file, see [example](example/TS_metadata.txt)

3. The name of the column in the metadata file that is the variable of interest. 

4. Set  min_abund: Minimum total abundance of a taxa required to be kept for further analysis

5. Set min_prevelance: Minimum number of samples that a taxa should be present in.

6. Output_dir: Output directory where all the results are saved.


An [analysis.Rmd](analysis.Rmd) file has the complete workflow that can be used for Metaphlan3 results and for pathway abundace tables from Humann3.

## Method details and shortcommings

Most of the methods used here have a lot of flexibility and the functions can take a lot more parameters as input. The table below explains briefly how these methods are used here. Also all output files are found in the output directory specified by the user. 

|Method      |  function Name| Output files|Information | What features are used     | What is not used|
| :---        |    :----:   |          ---: | ---:|---:|---:|
| DESeq2      | Deseq_fun()|deseq.csv|  [UserGuide](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html)     |  Works for two variables | Ramdom variables are added to the model like (Fixed + random). No interaction terms are allowed|
| EdgeR   |edgeR_fun() |edgeR.csv| [UserGuide](https://bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf)  |  Works for one|| variable    | Need to include more than one variable.s
| Maaslin2   |maaslin2_fun()| Maaslin2.csv + Maaslin2 dir| [UserGuide](https://www.bioconductor.org/packages/release/bioc/vignettes/Maaslin2/inst/doc/maaslin2.html)    |  works for two variables + random variables ,CSS normalization is default and no transformation is applied   | Not much missing but you can change the normalization method or add transformation    |
| MetagenomeSeq   |metagenomeSeq_fun()|Metagenomeseq.csv|[UserGuide](https://www.bioconductor.org/packages/devel/bioc/vignettes/metagenomeSeq/inst/doc/metagenomeSeq.pdf)|  works for one variable   |   It can be extended to multiple variables as per the user guide but looks more difficult to include   |  
| ANCOM2   |  AMCOM2_fun() |ANCOM2.csv|[UserGuide](https://github.com/FrederickHuangLin/ANCOM) | works for one variable  + 1 random variable     | structural zero feature is not being used right now but it can be use if the ANCOM2fun is called directly |
| ANCOM-BC  |ANCOMBC_fun() |ANCOMBC.csv + ACNOMBC_sample_fracs.csv|[UserGuide](http://www.bioconductor.org/packages/release/bioc/vignettes/ANCOMBC/inst/doc/ANCOMBC.html)   | works for two variables + 1 random variable      | Bias correction is not done but the sample fractions are saved if anyone wants to do it later.    |


## Required R packages
tidyverse v1.3.1
limma v3.50.0 
knitr v1.37
edgeR v3.36.0
Maaslin2 v1.8.0
metagenomeSeq v1.36.0
rstatix v0.7.0 
phyloseq v1.38.0
DESeq2 v1.34.0
ANCOMBC V1.4.0
DT v0.20


