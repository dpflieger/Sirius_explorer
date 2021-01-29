library(tidyverse)
library(data.table)

json_test <- "~/Projects/Gaquerel/Elser/Sirius/4_Tissue+Exudates-sirius301120_38/595_Tissue+Exudates-sirius301120_2955/trees/C30H48O16_[M+Na]+.json"
res <- jsonlite::fromJSON(json_test, flatten = TRUE)

fragments_mz <- res$fragments$mz
names(fragments_mz) <- res$fragments$id

fragments_mz["0"]

setDT(res$losses)

names(res$losses)
names(res$fragments)

res$losses[, source_mz := fragments_mz[as.character(source)]]
res$losses[, targets_mz := fragments_mz[as.character(target)]]

res$losses[, delta_mz := source_mz - targets_mz]


test <- fread("/Volumes/4TB/Users/dpflieger/Projects/Gaquerel/Elser/Sirius/4_Tissue+Exudates-sirius301120_38/formula_identifications.tsv")
colnames(test)
