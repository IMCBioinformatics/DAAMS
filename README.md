
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

For fixed variables and a random variable (like subject in repeated measures) the following methods are run:

1. DeSeq
2. Maaslin2

## Input

A species table is required. See example 
