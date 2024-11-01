rm(list = ls())

library(ade4)

args            = commandArgs(trailingOnly=TRUE)
experiment_name = args[1]
jDR_method      = args[2]
omics           = args[3]
B               = as.numeric(args[4])

cohorts = c("_female_", "_male_", "_pooled_")
cohort  = names(which(sapply(cohorts, function(x) grepl(pattern = x, x = experiment_name))))
cohort  = gsub(pattern = "_", replacement = "", x = cohort)

if (jDR_method == "JIVEORIGINE"){
  jDR_pattern = "JIVE"
} else {
  jDR_pattern = jDR_method
}

pattern = paste0(experiment_name, "_jdr_", jDR_pattern, "_omics_", omics, "_nbChunk_")
print(pattern)
files   = list.files(path = paste0("results/03_FeatureSelection/", cohort), pattern = pattern, full.names = T)
print(files)

results = c()
idx_ref = grep("_nbChunk_reference.RDS", files)
if (length(idx_ref) != 1){
  print(files)
  stop("Issue with the reference file")
}
for (i in c(1:length(files))[-idx_ref]){
  file           = files[i]
  tmp            = readRDS(file = file)
  idx_error      = which(sapply(tmp, class) == "try-error")
  tmp[idx_error] = NULL 
  results        = append(results, tmp)
  print(i)
  rm(tmp)
}
reference_result = readRDS(file = files[idx_ref])

idx_NA = which(is.na(unlist(sapply(results, function(x) x["best_cor"]))))
if (length(idx_NA) > 0){
  results[idx_NA] = NULL
}

check_nbComp = unique(unlist(sapply(results, function(x) x["best_cor"])))
if (length(check_nbComp) != 1){
  print(table(unlist(sapply(results, function(x) x["best_cor"]))))
  warning("More than one component is best accross bootstrap samples")
}

print(length(results))
if (length(results) > B){
  results[(B+1):(length(results))] = NULL
}
print(length(results))

omics = unlist(strsplit(x = omics, split = "_"))

if (jDR_method == "JIVE"){
  bootstrapped_features = lapply(1:length(results), function(x){
    reference = Reduce("rbind", reference_result[[1]]$features_w)
    cur       = Reduce("rbind", results[[x]]$features_w)
    rotation  = procuste(dfX = reference, dfY = cur, nf = dim(reference)[2])
    return(rotation$rotY)
  }) 
} else if (jDR_method == "RGCCA"){
  bootstrapped_features = lapply(1:length(results), function(x){
    output = c()
    for (omic in omics){
      reference   = reference_result[[1]]$features_w[[omic]]
      cur         = results[[x]]$features_w[[omic]]
      flip_or_not = sign(diag(cor(cur, reference)))
      output      = rbind(output, cur %*% diag(flip_or_not))
    }
    colnames(output) = colnames(reference)
    return(output)
  })
} else if (jDR_method == "JIVEORIGINE"){
  reference = Reduce("rbind", reference_result[[1]]$features_w)
  reference = as.matrix(reference[, reference_result[[1]]$best_cor])
  bootstrapped_features = sapply(1:length(results), function(x){
    cur           = Reduce("rbind", results[[x]]$features_w)
    cur           = as.matrix(cur[, results[[x]]$best_cor])
    flip_or_not   = drop(sign(cor(cur, reference)))
    cur           = cur*flip_or_not
    if (dim(cur)[1] != dim(reference)[1]){
      stop("Different number of vrairable than reference")
    }
    if ((sum(rownames(cur) == rownames(reference))) != dim(reference)[1]){
      stop("Variable not aligned to reference")
    }
    return(cur*flip_or_not)
  })
}

if (jDR_method == "JIVEORIGINE"){
  bootstrapped_features_best_cor           = bootstrapped_features
  rownames(bootstrapped_features_best_cor) = rownames(reference)
} else {
  bootstrapped_features_best_cor           = sapply(1:length(bootstrapped_features), function(x) bootstrapped_features[[x]][, reference_result[[1]]$best_cor])
  rownames(bootstrapped_features_best_cor) = apply(sapply(bootstrapped_features, rownames), 1, unique)
  
}

empirical_p_values          = apply(bootstrapped_features_best_cor, 1, function(x) max(min(sum(x < 0), sum(x >= 0))/max(sum(x < 0), sum(x >= 0)), 1/length(x)))
adjusted_empirical_p_values = p.adjust(p = empirical_p_values, method = "BH")
all_features_mean               = apply(bootstrapped_features_best_cor, 1, mean)
all_features_sd                 = apply(bootstrapped_features_best_cor, 1, sd)

final_results = c()
for (omic in omics){
  freq_all_features           = results[[1]]$features_w[[omic]][ , 1]*0
  all_features_selected       = unlist(sapply(results, function(x) x[["selectedFeatures"]][[omic]]))
  all_features_length         = unique(unlist(sapply(results, function(x) length(x[["selectedFeatures"]][[omic]]))))
  tmp_freq                    = table(all_features_selected)
  idx_tmp                     = match(names(tmp_freq), names(freq_all_features))
  freq_all_features[idx_tmp]  = tmp_freq
  features_selectionOccurence = sort(tmp_freq, decreasing = T)
  features_idx_end            = max(which(features_selectionOccurence == features_selectionOccurence[min(all_features_length)]))
  features_top_10             = names(features_selectionOccurence)[1:features_idx_end]
  features_adj_pval           = adjusted_empirical_p_values[match(rownames(results[[1]]$features_w[[omic]]), names(adjusted_empirical_p_values))]
  robust_features             = which((features_adj_pval <= 0.05) & (freq_all_features >= 0.8*B))
  features_p_values           = empirical_p_values[match(rownames(results[[1]]$features_w[[omic]]), names(empirical_p_values))]
  features_mean               = all_features_mean[match(rownames(results[[1]]$features_w[[omic]]), names(all_features_mean))]
  features_sd                 = all_features_sd[match(rownames(results[[1]]$features_w[[omic]]), names(all_features_sd))]
  print(robust_features)
  final_results[[omic]]       = list(freq_all_features           = freq_all_features,
                                     all_features_selected       = all_features_selected, 
                                     features_selectionOccurence = features_selectionOccurence,
                                     features_top_10             = features_top_10,
                                     features_adj_pval           = features_adj_pval,
                                     robust_features             = robust_features, 
                                     features_p_values           = features_p_values,
                                     features_mean               = features_mean,
                                     features_sd                 = features_sd)
}

saveRDS(final_results, file = paste0("results/04_ComputeIntervals/", 
                                     cohort, 
                                     "/final_results_", 
                                     experiment_name, 
                                     "_jdr_",
                                     jDR_method,
                                     "_omics_",
                                     paste(omics, collapse = "_"),
                                     "_with_extended_info.RDS"))
