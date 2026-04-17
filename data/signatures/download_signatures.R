library(ExperimentHub); library(rhdf5)
eh <- ExperimentHub()

# query(eh, c("signatureSearchData", "lincs"))
# lincs_path <- eh[['EH3226']]

#In this case the loaded data instance includes moderated Z-scores from DE analyses of 12,328 genes for 
# 8,140 compound treatments across a total of 30 cell lines corresponding to 45,956 expression signatures.
# This data set can be used by all set-based and correlation-based GESS methods provided by the signatureSearch
# package.


### using expression data from cmap
#query(eh, c("signatureSearchData"))
ssd = query(eh, c("signatureSearchData", "cmap_expr")); ssd$description; ssd$sourceurl
cmap_expr_path <- eh[["EH3224"]] 
#rhdf5::h5ls(cmap_rank_path)


print("Download completed!")
