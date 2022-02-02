
# DAAMS
Differential abundance analysis for shotgun metagenomics


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

Following inputs files required:

1. A comma separated species table is required. See [example](example/merged_species_table.csv)

2. A metadata file, see [example](example/TS_metadata.txt)

3. The name of the column in the metadata file that is the variable of interest. 

4. Set  min_abund: Minimum total abundance of a taxa required to be kept for further analysis

5. Set min_prevelance: Minimum number of samples that a taxa should be present in.

6. Output_dir: Output directory where all the results are saved.


An analysis.Rmd file has the complete workflow that can be used for Metaphlan3 results and for pathway abundace tables from Humann3.

## Good to know

Most of the methods used here have a lot of flexibility and the functions can take a lot more parameters as input. The table bellow explains the limitations and 

|Method      | Description | What features are used     | What is missing|
| :---        |    :----:   |          ---: | ---:|
| DeSeq2      |        |    |     |
| EdgeR   |    |       |     |
| Maaslin2   |    |       |     |
| MetagenomeSeq   |    |       |     |
| ANCOM2   |    |       |     |
| ANCOM-BC  |    |       |     |






