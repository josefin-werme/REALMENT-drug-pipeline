### ----------
# NOTE: this is part of a much larger (and partially outdated) script for processing the results
## For reference, I included the relevant functions mentioned in filter_signatures.R and signatureSearchProcessData.R
## The use of environments was done to optimise memory usage when procressing large amounts of results, but this isnt necessary for these functions alone. If preferred, you can just use regular variables instead. 
### ----------

## Packages
library(data.table)
library(XML); library(xml2) # to extract DrugBank data
library(parallel)
library(signatureSearch)


#### 1. FUNCTIONS ####
# check number of available cores
# parallel::detectCores()

## Base directory
base.dir = "/path/to/your/base/dir/" # TODO: set to your path

# VARIABLES
vars = list()
vars$snp = c("all","z196","z165")
vars$item = c("d29000","d20003-ATC", "d20003-ATC4N","icd10")
vars$method = c("spearman","spearman-top10","fisher-top10","pearson","spearman-top25","pearson-top25","pearson-top10")
vars$tissue = c("thyroid","brain_nucleus_accumbens","brain_cortex","whole_blood")


##### DEFINE ENVIRONMENTS #####
# create enviornments from variable names
create_envs = function(env.ids, parent.env) {
	tmp = sapply(env.ids,   function(x) { x = new.env(parent = parent.env) }) 
	for (i in names(tmp)) { parent.env[[i]] = tmp[[i]] }; rm(tmp) # add to datenv
}

## DEFINE ENVIRONMENTS
#define_envs = function() {
	### DRUG INFO DATA ###
	if (!exists("drug")) { drug = new.env(parent = globalenv()) }
	
	### RESULTS DATA ###
	### Level 1 -- parent environment for all data envs
	if (!exists("datenv")) {
		datenv = new.env(parent = globalenv())
		
		### Level 2 -- snp filtering
		# create environments from variable names in v
		create_envs(env.ids = vars$snp, parent.env = datenv)
		
		# for (snp in vars$snp) {  }
		
		## Level 3 -- phenotype item
		for (snp in vars$snp) {
			create_envs(env.ids = vars$item,  parent.env = datenv[[snp]])
			
			## Level 4 -- analysis method
			for (item in vars$item) {
				create_envs(env.ids = vars$meth,  parent.env = datenv[[snp]][[item]])
					
				## Level 5 -- tissue
				for (method in vars$method) {
					create_envs(env.ids = vars$tissue,  parent.env = datenv[[snp]][[item]][[method]])
				}
			}
		}
	} else {
		print("Warning: Data environment already exists, not creating new one")
	}
#}

#show_structure()
show_structure = function() {
	names = ls(datenv)
	print(paste("data:       ", paste(names, collapse = ", ")))
	
	v1 = names[1]; names = ls(datenv[[v1]])
	print(paste0("..", v1, ":       ", paste(names, collapse = ", ")))
	
	v2 = names[1]; names = ls(datenv[[v1]][[v2]])
	print(paste0("...", v2, ":   ", paste(names, collapse = ", ")))
	
	v3 = names[1]; names = ls(datenv[[v1]][[v2]][[v3]])
	print(paste0("...", v3, ":   ", paste(names, collapse = ", ")))
}

# set temp variable for an env based on v
set_env = function (v) { return(datenv[[v$snp]][[v$item]][[v$method]][[v$tissue]]) }
# ls(data[[v$snp]][[v$item]]$`spearman-top10`)

# show subset of ATC concerning codes or labels in a vector
show_atc = function(levels, names, class = NULL) {
	atc.sub = drug$atc[unlist(apply(drug$atc, 2, function(x) which(x %in% names))),] # match by columns, extract all rows with match
	headers = c(paste0("code",levels), paste0("class",levels))
	atc.sub = unique(atc.sub[, colnames(atc.sub) %in% headers])
	if (!is.null(class)) {
		atc.sub = subset(atc.sub, code1 == class)
	}
	return(atc.sub)
}
#show_atc(levels = c(1,3,5), names = "V09I")
#show_atc(levels = c(1,3,5), names = "V09I")

####### Create drug ID mapping data set for lincs2 and filter out experimental compounds #######
# NOTE: this doesnt need to be a function but its a bit neater like this (and prevents the code from running accidentally)
filter_compounds = function() {
	# compound meta data
	comp = fread(paste0(base.dir, "/your/path/to/compoundinfo_beta.txt"),        # TODO: set to your path
				 data.table = F)[,c("pert_id","cmap_name","compound_aliases")] 
	
	# Excluding drugs that don't have a 'common' ID or alias since I am assuming that these must be novel / experimental 
	# and would therefore not be relevant for evaluation here
	unknown.compounds = comp$cmap_name == comp$pert_id & comp$compound_aliases == "" # if cmap name is BRD format, and no alias exist
	comp = comp[!unknown.compounds,]
	
	### Format IDs
	# reason being that many cmap_name IDs are in BRD format. But I also dont want to just take the alias because they are a bit weirdly formatted sometimes
	cmapBRD = comp$cmap_name == comp$pert_id # get cmap_names that are BRD format
	comp$drug_id = comp$cmap_name # set new ID to cmap
	comp[cmapBRD,]$drug_id = comp[cmapBRD,]$compound_aliases # set cmap brd IDs to alias instead
	comp$drug_id
	
	comp[comp$drug_id == "quetiapine",] # lots of duplicates for many genes (due to other columns that were removed in beginning)
	dim(comp[comp$drug_id == "", ])
	
	### Creating final ID conversion file with only drugs of interest, and no duplicates
	comp_id = comp[,c('pert_id','drug_id')]
	comp_id = unique(comp_id)
	
	dim(comp_id)
	head(comp_id)
	
	write.table(comp_id, "compoundIDs_filtered.txt", row.names = F, quote = F, sep = '\t')
}

####### END #######




##### Extract data from DrugBank #####
### Drug ID synonyms ###
extract_drug_ids = function() {
	db_path = file.path(base.dir, "db_database.xml")   
	if (!file.exists(db_path)) { 
		print_parallel("ERROR: ATC code file path in 'extract_synonyms()' incorrect")
	}
	
	db = xmlParse(db_path)
	dbTop = xmlRoot(db) # get root node
	
	#drug_ids = list()
	drug$ids = list()
	#drug_ids = new.env(parent=globalenv())
	
	drug$ids = mclapply(1:as.numeric(xmlSize(dbTop)), mc.cores=10, function(i) {
		#drug$ids = mclapply(1:3, mc.cores=3, function(i) { # i=1
		name = tolower(xmlToList(dbTop[[i]][["name"]])) # common name
		syn = xmlToList(dbTop[[i]][["synonyms"]])
		syn = tolower(sapply(syn, "[[", "text"))
		names(syn) = NULL
		ids = list(name = name, syn = syn)
		return(ids)
	})
	
	# Format
	# this is a bit cumbersome, but I cannot do drug$ids[[name]] in loop so need to extract the names somehow)
	drug_names = sapply(drug$ids, "[[", "name")
	drug$ids = lapply(drug$ids, "[[", "syn")
	names(drug$ids) = drug_names
	
	## Check if all common drug names are also listed as synonyms
	# com_in_syn = sapply(names(drugbank_ids), function(x) x %in% drugbank_ids[[x]])
	# drugbank_ids[!com_in_syn]
	
	## Add common to synonym (and do unique)
	drug$ids = sapply(names(drug$ids), function(x) unique(c(drug$ids[[x]], x)))
	
	## Melt into data frame
	drug$ids = reshape2::melt(drug$ids)
	colnames(drug$ids) = c("synonym","common")
	
	# clean up memory
	rm(list = c("db", "dbTop", "drug_names"))
	gc()
}


###	 Harmonize drug names using DrugBank data	###
harmonise_ids = function(env, phen = NULL, colname) { # env = data$all[[item]][[method]]; colname = "drug_id"; print_args(harmonise_ids)
	## Extract drug ID synonyms from drugbank (only needs to be done first time function is run)
	if (!exists("ids", envir = drug)) {
		print_parallel("extracting drug IDs for mapping")
		extract_drug_ids()
	}
	
	# Match drug data with synonymous DrugBank IDs to harmonise
	env[[phen]]$global_id = drug$ids$common[match(tolower(env[[phen]][[colname]]), drug$ids$synonym)]
	env[[phen]]$failed = NA
	env[[phen]]$failed[is.na(env[[phen]]$global_id)] = env[[phen]][[colname]][is.na(env[[phen]]$global_id)]
}

harmonise_ids_vec = function(ids) { # env = data$all[[item]][[method]]; colname = "drug_id"; print_args(harmonise_ids)
	## Extract drug ID synonyms from drugbank (only needs to be done first time function is run)
	if (!exists("ids", envir = drug)) {
		print_parallel("extracting drug IDs for mapping")
		extract_drug_ids()
	}
	
	# Match drug data with synonymous DrugBank IDs to harmonise
	global_id = drug$ids$common[match(tolower(ids), drug$ids$synonym)]
	return(global_id)
}


