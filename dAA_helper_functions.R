# Requires the following packages
library(tidyverse)
library(limma)
library(knitr)
library(edgeR)
library(Maaslin2)
library(metagenomeSeq)
library(rstatix)
library(phyloseq)
library(DESeq2)
# For ANCOM this file is required
source("~/projects/Differential_analysis_smg/ancom_v2_1.R")


#' Title
#'
#' @param min_abun  This is the minimum total abundance each taxa should have
#' @param min_prev Minimum proportion of samples that should have a taxa to keep it between 0 and 1
#' @param metadata A metadata file with rownames that match the column names of ASV_table
#' @param ASV_table Rows are taxa  and ony species are present
#' @param fixed_vars The variables that need to tested for. Maximum is 2
#' @param rand_vars  Variables like subjects that are measured multiple times can be placed here
#' @param maaslin2_normalization Normalization method for maaslin2 can be one of c("TSS", "CLR", "CSS", "NONE", "TMM") default = CSS
#' @param maaslin2_transform  Tranformation method for maaslin2 can be one of c("LOG", "LOGIT", "AST", "NONE") default = None
#'
#' @return returns a dataframe with results from suited DA methods
#' @export
#'
#' @examples
wrapper_daa <-
  function(min_abun = 0,
           min_prev = 0.1,
           metadata,
           ASV_table,
           fixed_vars,
           rand_vars = NULL,
           maaslin2_normalization="CSS",
           maaslin2_transform="None") {
    ## metadata
    metadata <- tryCatch(
      metadata %>% select(fixed_vars, rand_vars),
      error = function(e)
        stop(paste(
          "fixed_vars or rand_vars not found in provided metadata "
        ))
    )
    groupings <- metadata %>% drop_na()
   
    print(c("Testing for variables:", fixed_vars))
    print(c("Random variables used:",rand_vars))
    
    
    ## ASV table filtering based on min_abundance and min prev
    print(c("ASV_table size before filterting:", dim(ASV_table)))
    ASV_table <-
      ASV_table[which(rowSums(ASV_table) > min_abun),]
    print(c("ASV_table size after min_abun filterting:", dim(ASV_table)))
    ASV_table <-
      ASV_table[which(apply(ASV_table, 1, function(x)
        length(which(x > 0)) > ncol(ASV_table) * min_prev)),]
    print(c("ASV_table size after min_prev filterting:", dim(ASV_table)))
    
    rows_to_keep <-
      intersect(colnames(otu_table), rownames(groupings))
    groupings <- groupings[rows_to_keep, , drop = F]
    ASV_table <- ASV_table[, rows_to_keep]
    idx = -which(rowSums(ASV_table) == 0)
    if (length(idx) > 0)
      ASV_table <- ASV_table[idx, ]
    print(c("ASV_table size after filtering for samples present in metadata:", dim(ASV_table)))
    
    
    
    if (length(fixed_vars) == 1) {
      if(!is.null(rand_vars)){
      fixed_formula_deseq <- paste0("~", fixed_vars,"+",rand_vars)
      fixed_formula_deseq <- as.formula(fixed_formula_deseq)
      } else {
        fixed_formula_deseq <- paste0("~", fixed_vars)
        fixed_formula_deseq <- as.formula(fixed_formula_deseq)
      }
      
      print(c("Formula being used for Deseq",fixed_formula_deseq))
      
      all_res <- left_join(deseq_fun(ASV_table, groupings,fixed_formula_deseq),
                           edgeR_fun(ASV_table, groupings),
                           by = "microbe") %>%
        left_join(., maaslin2_fun(abundance_table = ASV_table,
                                  metadata = groupings,
                                  normalization = maaslin2_normalization,
                                  transform = maaslin2_transform,
                                  fixed_effect = fixed_vars,
                                  random_effects = rand_vars,
                                  output = "test"), by = "microbe") %>%
        left_join(., metagenomeSeq_fun(ASV_table, groupings), by = "microbe") %>%
        left_join(., ANCOM2_fun(ASV_table, groupings), by = "microbe")
      
    }
    else if (length(fixed_vars) == 2) {
      if(!is.null(rand_vars)){
      fixed_formula_deseq <-
        paste0("~", fixed_vars[1], "+", fixed_vars[2],"+",rand_vars)
      }else{
        fixed_formula_deseq <-
          paste0("~", fixed_vars[1], "+", fixed_vars[2])
      }
      
    } else if (length(fixed_vars) > 2) {
      stop(
        paste(
          "More than 2 variables in fixed effects are not allowed. Please reduce the number."
        )
      )
    }

    
    if (length(fixed_vars) + length(rand_vars) > 1) {
      fixed_formula_deseq <- as.formula(fixed_formula_deseq)
      
      all_res <-
        left_join(
          deseq_fun(ASV_table, groupings, fixed_formula_deseq),
          maaslin2_fun(
            abundance_table = ASV_table,
            metadata = groupings,
            normalization = maaslin2_normalization,
            transform = maaslin2_transform,
            fixed_effect = fixed_vars,
            random_effects = rand_vars,
            output = "test"
          ),
          by = c("microbe","Comparison")
        )
    }
    return(all_res)
  }




#Run Deseq2
#' Title
#'
#' @param abundance_table = A table with species as rows and samples as columns. Must be raw values
#' @param metadata = a dataframe with a single variable to test on with only two factors
#'
#' @return
#' @export
#'
#' @examples

deseq_fun <- function(abundance_table, metadata, formula) {
  dds <- DESeqDataSetFromMatrix(countData = abundance_table,
                                        colData = metadata,
                                        design = formula)
  print("In Deseq")

  
  dds_res <- DESeq(dds, sfType = "poscounts")
  
  if (length(resultsNames(dds_res)) > 2) {
    res1 <-
    results(
        dds_res,
        tidy = T,
        format = "DataFrame",
        name = resultsNames(dds_res)[2]
      )
    res1$Comparison <- paste0(strsplit2(resultsNames(dds_res)[2],split = "_")[,1:2],collapse = "")
      
    
    res2 <-
      results(
        dds_res,
        tidy = T,
        format = "DataFrame",
        name = resultsNames(dds_res)[3]
      )
    res2$Comparison <- paste0(strsplit2(resultsNames(dds_res)[3],split = "_")[,1:2],collapse = "")
    res <- bind_rows(res1, res2)
    
  }
  else{
  
    res <- results(dds_res, tidy = T, format = "DataFrame")
    res$Comparison <- resultsNames(dds_res)[2]
  }
  
  
  
  deseq <-
    res %>% dplyr::rename(microbe = row) %>%
    dplyr::rename(deseq_LFC = log2FoldChange) %>%
    dplyr::rename(deseq_padj = padj) %>%
    select(microbe, Comparison , deseq_LFC, deseq_padj)
  
  #  return(deseq)
  return(deseq)
}



edgeR_fun <-
  function(abundance_table,
           metadata,
           adjust.method = "fdr") {
    ### Taken from phyloseq authors at: https://joey711.github.io/phyloseq-extensions/edgeR.html
    phyloseq_to_edgeR = function(physeq, group, method = "RLE", ...) {
      # Enforce orientation.
      if (!taxa_are_rows(physeq)) {
        physeq <- t(physeq)
      }
      x = as(otu_table(physeq), "matrix")
      # Add one to protect against overflow, log(0) issues.
      x = x + 1
      # Check `group` argument
      if (identical(all.equal(length(group), 1), TRUE) &
          nsamples(physeq) > 1) {
        # Assume that group was a sample variable name (must be categorical)
        group = get_variable(physeq, group)
      }
      # Define gene annotations (`genes`) as tax_table
      taxonomy = tax_table(physeq, errorIfNULL = FALSE)
      if (!is.null(taxonomy)) {
        taxonomy = data.frame(as(taxonomy, "matrix"))
      }
      # Now turn into a DGEList
      y = DGEList(
        counts = x,
        group = group,
        genes = taxonomy,
        remove.zeros = TRUE,
        ...
      )
      # Calculate the normalization factors
      z = edgeR::calcNormFactors(y, method = method)
      # Check for division by zero inside `calcNormFactors`
      if (!all(is.finite(z$samples$norm.factors))) {
        stop(
          "Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument"
        )
      }
      # Estimate dispersions
      return(estimateTagwiseDisp(estimateCommonDisp(z)))
    }
    
    OTU <- phyloseq::otu_table(abundance_table, taxa_are_rows = T)
    sampledata <- phyloseq::sample_data(metadata, errorIfNULL = T)
    phylo <- phyloseq::merge_phyloseq(OTU, sampledata)
    
    test <-
      phyloseq_to_edgeR(physeq = phylo, group = colnames(metadata)[1])
    
    et = exactTest(test)
    
    tt = topTags(
      et,
      n = nrow(test$table),
      adjust.method = adjust.method,
      sort.by = "PValue"
    )
    res <- tt@.Data[[1]]
    edgeR_res <-
      res %>%  rownames_to_column("microbe") %>%
      dplyr::rename(edgeR_LFC = logFC) %>%
      dplyr::rename(edgeR_padj = FDR) %>%
      select(microbe, edgeR_LFC, edgeR_padj)
    return(edgeR_res)
  }


maaslin2_fun <-
  function(abundance_table,
           metadata,
           normalization = "CSS",
           transform = "NONE",
           fixed_effect,
           random_effects = NULL,
           output = "temp") {
    abundance_table <-
      data.frame(
        t(abundance_table),
        check.rows = F,
        check.names = F,
        stringsAsFactors = F
      )
    
    fit_data <-
      Maaslin2(
        input_data = abundance_table,
        input_metadata = metadata,
        output = output,
        transform = transform,
        normalization = normalization,
        fixed_effects = fixed_effect,
        random_effects = random_effects,
        standardize = FALSE,
        plot_heatmap = F,
        plot_scatter = F
      )
    
    M2 <-
      fit_data$results %>% dplyr::rename(microbe = feature) %>%
      dplyr::rename(M2_coef = coef) %>% dplyr::rename(M2_padj = qval) %>%
      dplyr::rename(Comparison = name) %>% dplyr::rename(Num_not_zero = N.not.zero)%>%
    select(microbe, M2_coef, M2_padj,Comparison, Num_not_zero)
 
    return(M2)
  }


metagenomeSeq_fun <- function(abundance_table, metadata) {
  data_list <- list()
  data_list[["counts"]] <- abundance_table
  data_list[["taxa"]] <- rownames(abundance_table)
  
  pheno <- AnnotatedDataFrame(metadata)
  pheno
  counts <- AnnotatedDataFrame(abundance_table)
  feature_data <- data.frame("ASV" = rownames(abundance_table),
                             "ASV2" = rownames(abundance_table))
  feature_data <- AnnotatedDataFrame(feature_data)
  rownames(feature_data) <- feature_data@data$ASV
  
  
  test_obj <-
    newMRexperiment(
      counts = data_list$counts,
      phenoData = pheno,
      featureData = feature_data
    )
  
  p <- metagenomeSeq::cumNormStat(test_obj, pFlag = F)
  p
  
  test_obj_norm <- cumNorm(test_obj, p = p)
  
  fromula <-
    as.formula(paste(~ 1, colnames(metadata)[1], sep = " + "))
  pd <- pData(test_obj_norm)
  mod <- model.matrix(fromula, data = pd)
  regres <- fitFeatureModel(test_obj_norm, mod)
  
  res_table <-
    MRfulltable(regres, number = length(rownames(abundance_table))) %>% rownames_to_column("microbe") %>%
    dplyr::rename(MSeq_OR = oddsRatio) %>%
    dplyr::rename(MSeq_LFC = logFC) %>%
    dplyr::rename(MSeq_padj = adjPvalues) %>%
    select(microbe, MSeq_OR, MSeq_LFC, MSeq_padj)
  
  return(res_table)
}



ANCOM2_fun <-
  function(abundance_table,
           metadata,
           rand_formula = NULL,
           p_adj_method = "BH",
           alpha = 0.05,
           out_cut = 0.05,
           zero_cut = 0.90,
           lib_cut = 1000) {
    metadata$Sample <- rownames(metadata)
    prepro <-
      feature_table_pre_process(
        feature_table = abundance_table,
        meta_data = metadata,
        sample_var = 'Sample',
        group_var = NULL,
        out_cut = out_cut,
        zero_cut = zero_cut,
        lib_cut = lib_cut,
        neg_lb = FALSE
      )
    
    feature_table <- prepro$feature_table
    metadata <- prepro$meta_data
    struc_zero <- prepro$structure_zeros
    
    #run ancom
    main_var <- colnames(metadata)[1]
    p_adj_method = p_adj_method
    alpha = alpha
    adj_formula = NULL
    rand_formula = rand_formula
    res <-
      ANCOM(
        feature_table = feature_table,
        meta_data = metadata,
        struc_zero = struc_zero,
        main_var = main_var,
        p_adj_method = p_adj_method,
        alpha = alpha,
        adj_formula = adj_formula,
        rand_formula = rand_formula
      )
    
    ANCOM2_res <-
      res$out %>% dplyr::rename(microbe = taxa_id) %>%
      dplyr::rename(ANCOM2_W = W) %>%
      dplyr::rename(ANCOM2_0.6 = detected_0.6) %>%
      select(microbe, ANCOM2_W, ANCOM2_0.6)
    
    return(ANCOM2_res)
  }


wilcox_AST <- function(abundance_table, metadata) {
  ### WILCOX + AST
  
  # Arc Sine Square Root Transformation
  AST <- function(x) {
    return(sign(x) * asin(sqrt(abs(x))))
  }
  abundance_table_rel <- apply(abundance_table, 2, function(x)
    x / sum(x))
  gathered <-
    apply(abundance_table_rel, 2, AST) %>%
    as.data.frame() %>%
    rownames_to_column("microbe") %>%
    gather(key = Samples, value = abundance,-microbe)
  metadata$Samples <- rownames(metadata)
  gathered <- left_join(gathered, metadata, by = "Samples")
  WRS_res <-
    left_join(
      gathered %>% rstatix::group_by(microbe) %>%
        wilcox_test(abundance ~ Groupings, detailed = T) %>%
        adjust_pvalue(method = "fdr"),
      gathered %>%
        rstatix::group_by(microbe) %>%
        wilcox_effsize(abundance ~ Groupings),
      by = "microbe"
    ) %>%
    select(microbe, effsize, p.adj) %>%
    dplyr::rename(WRS_effsize = effsize) %>%
    dplyr::rename(WRS_padj = p.adj)
  
  return(WRS_res)
}
