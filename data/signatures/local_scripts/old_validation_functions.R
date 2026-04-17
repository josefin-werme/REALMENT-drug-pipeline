## Packages
library(data.table)
library(XML); library(xml2) # to extract DrugBank data
library(parallel)
library(signatureSearch)


#### 3. To do list ####
# TODO: check if package MOA differs from that in compoundinfo_beta.txt
# TODO: Analyse automatically 
# Different cell types - 
# TODO: redefine the ATC code phenotypes with the new data - or select some more cardiac phenotypes
# TODO: Perform enrichment analysis with KEGG
# TODO: analyse with fisher
# TODO: analyse diabetes, parkinsons, hypothyroidism
# TODO: another way of using ATC could be to see whether you can predict within class ?
# TODO:  >>>>>> >>>>> MAIN !!!!!! <<<<<< <<<<<< 
# TODO: write function to use the annotated MOAss directly for enrichmnet
# TODO: also instead of averaging across tissues, why not just subset????
#### end ####


#### 1. FUNCTIONS ####


# check number of available cores
# parallel::detectCores()

## Base directory
base.dir = "~/surfdrive/POSTDOC/Genes_to_Drugs/"

# VARIABLES
vars = list()
vars$snp = c("all","z196","z165")
vars$item = c("d29000","d20003-ATC", "d20003-ATC4N","icd10")
vars$method = c("spearman","spearman-top10","fisher-top10","pearson","spearman-top25","pearson-top25","pearson-top10")
vars$tissue = c("thyroid","brain_nucleus_accumbens","brain_cortex","whole_blood")

# v$snp = c("all","z96","z65")
# v$item = c("d29","d23")
# v$meth = c("spear","spear10","fish10")

# vmap = list()
# vmap$short = v
# vmap$full$snp = c("","z196","z165")
# vmap$full$item = c("d29000","d20003-ATC")
# vmap$full$meth = c("spearman","fisher")


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
	comp = fread(paste0(base.dir, "/DrugSignatures/LINCS2020/compoundinfo_beta.txt"), data.table = F)[,c("pert_id","cmap_name","compound_aliases")] 
	
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




##### ~ BASIC PROCESSING FUNCTIONS FOR RESULTS ~ #####

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

### Annnotate results ###
#  Extracts drug IDs, cell types, cell subtypes, and add MOA
annotate_results = function(env, phen, v) {
	if (!exists("cells", envir = drug)) { drug$cells = fread(paste0(base.dir, "/DrugSignatures/LINCS2020/cellinfo_beta.txt"), data.table = F) }
	
	if (v$dataset == "nath") {
		if (!exists("josy.mat")) { load(paste0(base.dir,"DrugSignatures/Nathans_signatures/normal.drug.info.RData")) }
		env[[phen]]$pert_id = josy.info$pert_id[match(env[[phen]]$signature, josy.info$sig_id)]
		env[[phen]]$cell_id = josy.info$cell_iname[match(env[[phen]]$signature, josy.info$sig_id)]
	} else {
		if (!exists("comp", envir=drug)) { 
			drug$comp = fread(paste0(base.dir, "/DrugSignatures/LINCS2020/compoundIDs_filtered.txt"), data.table = F) # d$comp[d$comp$pert_id %in% names(tail(sort(table(d$comp$pert_id)))),]  # TODO: How to deal with duplicate pert_ids ??? exclude them entirely?? (there are not many so perhaps not huge issue for now, but good to deal with soon)
		}
		vars = strsplit(env[[phen]]$signature, "__")
		env[[phen]]$pert_id = sapply(vars, "[[", 1)
		env[[phen]]$cell_id = sapply(vars, "[[", 2)
	}
	env[[phen]]$celltype = drug$cells[match(env[[phen]]$cell_id, drug$cells$cell_iname),]$cell_type
	env[[phen]]$subtype = drug$cells[match(env[[phen]]$cell_id, drug$cells$cell_iname),]$subtype
	env[[phen]]$drug_id = drug$comp[match(env[[phen]]$pert_id, drug$comp$pert_id),]$drug_id
	
	# Annotate MOA
	if (!exists("clue_moa_list")) { data("clue_moa_list") }
	env[[phen]] = addMOA(env[[phen]], "drug_id", clue_moa_list)
}

### Read in and process a single phenotype ###
# annotating drugs, cells, MOAss + harmonising IDs
# process_phenotype = function(item, phen) {
# 	# item = "d29000"; method="fisher-top10"; dataset = "lincs2"; tissue="whole_blood"; phen="depression"
# 	## Read in
# 	data = fread(paste0(paste0("whole_blood_",item,"_"), phen, "_lincs2_spearman.out"), data.table=F)
# 	## Annotate
# 	data = annotate_results(data, T)
# 	## Harmonise drug names
# 	data = harmonise_ids(data, "drug_id")
# 	return(data)
# }

# ### Read in by filename ###
# process_resultsfile = function(env, phen, filename) { # filename = input_files[i]; phen = phenos[i]
# 	env[[phen]] = fread(paste0(base.dir,"validation/signatures/",filename), data.table=F) ## Read in
# 	
# 	# print(phen)
# 	# print(filename)
# 	# print(head(env[[phen]]))
# 	
# 	annotate_results(env, phen, T) ## Annotate
# 	try = try(harmonise_ids(env, phen,"drug_id")) ## Harmonise drug names
# 	if (class(try)=="try-error") { stop("!!!! NOTE: You need to add current condition to 'vars' and re-define the datenv environment !!!!") }
# }
# env[[phen]]$global_id

### Process all phenotypes for a specific item / method / tissue / data set ###
# item = v$item; method=v$method; dataset = "lincs2"; tissue="whole_blood";n.cores=5; 
# env = set_env(v)
# asdf
process_results = function(env, v, phenos = NULL, n.cores = 10, test=F) {
	# print_args(process_results) # print arguments (debug)
	# n.cores = as.numeric(n.cores)
	input_files = get_files(v, phen_sub = phenos)
	#get_files(v, phen_sub)
	if (is.null(phenos)) { phenos = get_phenos(v) }
	
	for(i in 1:length(input_files)) {	#
	#tmp = mclapply.noreturn(1:length(phenos), mc.cores = n.cores, function(i) {
		# i=1
	#	print_parallel(phenos[i])
		phen = phenos[i]
		env[[phen]] = fread(paste0(base.dir,"validation/signatures/",input_files[i]), data.table=F) ## Read in
		annotate_results(env, phen, v) ## Annotate
		try = try(harmonise_ids(env, phen,"drug_id")) ## Harmonise drug names
		if (class(try)=="try-error") { stop("!!!! NOTE: You need to add current condition to 'vars' and re-define the datenv environment !!!!") }
		
	}#)
	#for(phen in phenos) { env[[phen]] = tmp[[phen]] }
}

get_files = function(v, phen_sub = NULL) {
	all_files = dir(paste0(base.dir,"validation/signatures")) # list all results files
	prefix = paste0(v$tissue, "_", v$item,"_") # defining this separately so phenotypes can be easily extracted below
	
	relevant_files = grep(prefix,        # filter by args
						  grep(paste0(v$dataset,"_",v$method, ".out"),  # need to be general here since not all have snp filtering
						  all_files, value = T), value = T)
	#print(relevant_files)
	
	# grep based on snp filtering ("all" was never labelled as it was the first, hence why it is selected by excludion)
	if(v$snp == "all") {
		for (alt.snp in vars$snp[vars$snp!="all"]) {
			relevant_files = relevant_files[-grep(alt.snp, relevant_files)]
		}
	} else {
		relevant_files = grep(paste0("_",v$snp,"_"), relevant_files, value=T)
	}
	#print(relevant_files)
	if (!is.null(phen_sub)) {
		# extract phenotypes from file names
		all_phenos = get_phenos(v)
		# subset filenames to phen_sub
		relevant_files = relevant_files[all_phenos %in% phen_sub]
	}
	#print(relevant_files)
	return(relevant_files)
}

get_phenos = function(v, atc_subset = NULL) {
	relevant_files = get_files(v)
	# get phenos
	phenos = gsub(paste0(v$tissue,"_",v$item,"_"), "",
				gsub(paste0("_",v$dataset,"_",v$method,".out"),"", relevant_files))
	if (v$snp != "all") {
		phenos = gsub(paste0("_",v$snp), "", phenos)
	}
	if (!is.null(atc_subset)) {
		phenos = grep(paste0("^",atc_subset), phenos, value=T) # if e.g. only phenotypes from N should be selected
	}
	return(phenos)
}

##### END #####


#### UTILS ####
#### multi.match 
# can output multiple columns, and will collapse column values separately. Outputs list by default. 
multi.match = function(values, match.df, match.col, out.col, collapse = T) {
	matches = lapply(values, function(x) {
		unique(
			match.df[which(match.df[,match.col] %in% x), out.col]
		) 
	})
	
	
	if (collapse) {
		matches = lapply(matches, function(x) {
			apply(x, 2, paste, collapse="; ")
		})
	}
	
	names(matches) = values
	return(matches)
}
# # get all rows that match per drug
# match.col = "global_id"; out.col = "code3"; match.df = drug$atc; values = dat$global_id; out.col = c("code3","class3")
# values=mania$ids; match.df= drug$atc; match.col= "global_id"; out.col = c("code3","class3"); collapse=T
# ----------
# Functions
print_args = function(func) {
	arg_ids = names(formals(func))
	arg_defs = paste0("'", sapply(arg_ids, function(x) eval(parse(text = x))), "'")
	args = cbind(arg_ids, arg_defs)
	print(paste0(apply(args, 1, paste, collapse=" = "), collapse = "; "))
}

print_parallel <- function(...){
	system(sprintf('echo "\n%s\n"', paste0(..., collapse="")))
}

# Environment related
# return all objects in current environment
list_objects = function(env = globalenv()) {
	# TODO: make work for other envs than global... 
	
	# if (!exists("utils")) utils = new.env(parent = globalenv())
	# utils$objs = as.data.frame(t(t(sort(sapply(ls(globalenv()), function(x) mode(eval(parse(text = x))))))))
	# utils$objs = data.frame(ID = rownames(utils$objs), TYPE = utils$objs$V1)
	# utils$objs
	objs = as.data.frame(t(t(sort(sapply(ls(env), function(x) mode(eval(parse(text = x))))))))
	objs = data.frame(ID = rownames(objs), Type = objs$V1)
	return(objs)
}
# mode(eval(parse(text = objects(data))))

# return all runctions in enviornment
list_functions = function(env = globalenv()) {
	print(subset(list_objects(env), Type == "function")$ID)
}

# return all runctions in enviornment
list_environments = function(env = globalenv()) {
	print(subset(list_objects(env), Type == "environment")$ID)
}

clear_environments = function(env = globalenv()) {
	envs = list_environments(env)
	rm(list=envs)
}

read_atc = function() {
	#if (!exists("atc", envir = drug)) {
		drug$atc = fread(paste0(base.dir,"atc_jw.txt"), data.table = F)
		drug$atc[,c("class1","class2","class3")] = apply(drug$atc[,c("class1","class2","class3")], 2, stringr::str_to_title)
		harmonise_ids(drug, "atc", "class5")
	#}
}

# list environments
# search(globalenv())
# ls(globalenv())
# ?base::search
# 
# ?is.function()
# ?is.function()
# 
# Sys.getenv()
# 
# objects()
# base::ls(all.names = T)

##### END #####



##### AVERAGE ACROSS CELL TYPES #####
# Average results across cell lines
## NOTE: this result HAS to have gone through the processing step first
# average_cells
aggregate_results = function(x, type) {
	# x = "normal soft tissue sample"
	# TODO: correlate estimates across subtypes? (to get some confidence in score)
	print(paste(type, x))
	head(data[data[[type]] == x,])
	
	est = aggregate(est ~ pert_id, data[data[[type]] == x,], median)
	p = aggregate(p ~ pert_id, data[data[[type]] == x,], median)
	if (! all(est$drug_id == p$drug_id)) stop("match ids") # safety check (to avoid having to match)
	
	# annotate
	dat = annotate_results(data.frame(est, p = p$p), full = F)
	return(dat[order(dat$p),])
}

## aggregate_cells
average_results = function(data) {
	celltypes = unique(data$celltype)
	subtypes = unique(subset(data, celltype == "normal")$subtype) # NO CANCER subtypes
	#subtypes = table(subset(data, celltype == "normal")$subtype) # NO CANCER subtypes
	#names(subtypes[subtypes > 5])
	
	ct = lapply(celltypes, function(x) aggregate_results(x, "celltype")); names(ct) = celltypes
	st = lapply(subtypes, function(x) aggregate_results(x, "subtype")); names(st) = subtypes
	return(list(orig = data, celltypes = ct, normal_subtypes = st))
}
##### END #####


###### ATC ENRICHMENT ANALYSIS FUNCTION ######
# allows subsetting by cell type and cell subtype
# lvl=1; lvl.sub = NA; cellType = NA; subType = NA; p.thresh = .05; neg.only=F
fisher = function(qtop=NULL, qbot=NULL, list=NULL) { # vectors of T/F
	#if (!is.null(list)) { qtop = list$qtop; qbot = list$qbot }
	if (!is.null(list)) { qtop = list[[1]]; qbot = list[[2]] }
	fishdat = data.frame(
		"sig" = c(sum(qtop == T), sum(qtop == F)),
		"not_sig" = c(sum(qbot == T), sum(qbot == F)),
		row.names = c("related", "not_related"),
		stringsAsFactors = FALSE
	)
	fish = try(fisher.test(fishdat, alternative = "greater"), silent=T)
	return(data.frame(
		odds_ratio = round(as.numeric(fish$estimate), 2),
		p = signif(fish$p.value, 2)))
}

atc_enrichment = function(env, phen, v, lvl, lvl.sub = NA, cellType = NA, subType = NA, p.thresh = .05, neg.only) {
	if (!exists("atc", envir = drug)) { read_atc() }
	#print_args(atc_enrichment)
	
	## Subset the results to only drugs listed in the ATC file
	tmp.df = env[[phen]][env[[phen]]$global_id %in% drug$atc$global_id,]
	
	if (neg.only) { print("WARNING: neg.only=T"); tmp.df = tmp.df[tmp.df$est < 0,] } # including only negative correlations
	
	## Subset the data to desired celltype / subtype
	if (!is.na(cellType)) {
		tmp.df = subset(tmp.df, celltype == cellType)
	}
	if (!is.na(subType)) {
		tmp.df = subset(tmp.df, subtype == subType)
	}
	
	#lvl.sub = "N"
	if (!is.na(lvl.sub)) {
		atc_codes = unique(subset(drug$atc, code1 == lvl.sub)[[paste0("code",lvl)]])
	} else {
		atc_codes = unique(drug$atc[[paste0("code",lvl)]])
	}
	
	## Iterate over ATC codes to evaluate enrichment
	#for (i in out$code) {
	out = do.call(rbind, mclapply(atc_codes, mc.cores=10, function(i) {
		# i = "N05A"
		# i = "N05AF04"
		atc.sub = drug$atc[drug$atc[[paste0("code",lvl)]] == i,] # TODO: precompute theese
		
		qtop = subset(tmp.df, p < p.thresh)$global_id %in% atc.sub$global_id
		qbot = subset(tmp.df, p > p.thresh)$global_id %in% atc.sub$global_id
		
		fish = fisher(qtop = qtop, qbot = qbot)
		
		#out[out$code == i,]$odds_ratio = round(as.numeric(fish$estimate), 2)
		#out[out$code == i,]$p = signif(fish$p.value, 2)
		data.frame(
			atc = i,
			odds_ratio = fish$odds_ratio,
			p = fish$p
		)
	}))
	out = out[order(out$p),] # order by significance
	return(out)
}




# flag "by.atc" expects F or a lvl1 code
# sign expects: F, "-", "+"
moa_enrichment = function(env, phen, v, cellType = NA, subType = NA, p.thresh = .05, by.atc = F, sign = NA) { # p.thresh = .05
	## Subset data to only drugs with MOA annotation
	tmp.df = env[[phen]][!is.na(env[[phen]]$MOAss),]
	
	if (!is.na(sign)) {
		if (! sign %in% c("-","+")) { stop("Invalid sign argument") }
		# head(tmp.df, 10)[,c("est","global_id")]; sign="+"
		sign1 = as.numeric(paste0(sign,1)) # -1 / 1
		tmp.df = subset(tmp.df, est * sign1 > 0) # subsetting by multiplying the est by -1 or 1
	}
	
	## If ATC code is specified, subset data to drugs belonging to code of interest
	if (by.atc != F) {
		atc.drugs = subset(drug$atc, code1 == by.atc)$global_id
		tmp.df = subset(tmp.df, global_id %in% atc.drugs)
	}
	## Subset the data to desired celltype / subtype
	if (!is.na(cellType)) {
		tmp.df = subset(tmp.df, celltype == cellType)
	}
	if (!is.na(subType)) {
		tmp.df = subset(tmp.df, subtype == subType)
	}
	## Prep output data frame
	out = data.frame(moa = unique(unlist(strsplit(tmp.df$MOAss, "; "))),
					 odds_ratio = NA,
					 p = NA)
	
	## Splitting data set
	top.moa = unlist(strsplit(subset(tmp.df, p < p.thresh)$MOAss, "; "))
	bot.moa = unlist(strsplit(subset(tmp.df, p > p.thresh)$MOAss, "; "))
	
	## Iterate over MOAs to evaluate enrichment
	for (i in out$moa) { # i= top.moa[1]
		qtop = top.moa == i
		qbot = bot.moa == i
		fish = fisher(qtop = qtop, qbot = qbot)
		
		out[out$moa == i,]$odds_ratio = fish$odds_ratio
		out[out$moa == i,]$p = fish$p
	}
	out = out[order(out$p),] # order by significance
	return(out)
}

# This is only convenient if you want to run all ATC classes (i.e. code = NULL)
# .. otherwise you can just as well run the original function with the code arg set
moa_by_atc = function(env, phen, v, code = NA, p.thresh = .05, sign = NA) { # code = "N"; p.thresh = .5; sign = NA
	# if no class specified, run all
	if (is.na(code)) {
		moatc = list()
		codes = unique(drug$atc$code1)
		for (c in codes) {
			moatc[[c]] = moa_enrichment(env, phen, v, by.atc = c, p.thresh = p.thresh)
		}
		names(moatc) = lapply(codes, function(x) paste0(show_atc(1,x), collapse=": "))
	} else {
		moatc = moa_enrichment(env, phen, v, by.atc = code, p.thresh = p.thresh) 
	}
	return(moatc)
}

###### END ######




##### Extract data from DrugBank #####
### Drug ID synonyms ###
extract_drug_ids = function() {
	db_path = file.path(base.dir, "ATC code extraction/Nathans/db_database.xml")
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

### Pathway information ###

# Note: these seem to just be drug related pathways ? 
# pathways = list()
# for (i in 1:as.numeric(xmlSize(dbTop))) {
# 	# common name
# 	drug = xmlToList(dbTop[[i]][["name"]]) 
# 	# get pathway info
# 	path = xmlToList(dbTop[[i]][["pathways"]])
# 	# store
# 	if (!is.null(path)) {
# 		pathways[[drug]] = paste0(sapply(path, "[[", "name"), collapse = " | ")
# 	}
# }

##### END #####





#### FORMAT ATC ENRICHMENT RESULTS ####
map_atc_labels = function(names) {
	names.c1 = names[names %in% paste0(LETTERS, 0)] # extract all names like N0/H0 etc
	names[names %in% names.c1] = gsub("0","",names.c1) # remove the 0
	atc_ids = lapply(names, atc_class); names(atc_ids) = names # map code to label
	return(atc_ids)
}

# Add ATC class to phenotype names ### (only for 23000-ATC)
format_enrich_ids = function(orig_ids) {
	atc_labels = map_atc_labels(orig_ids) # map ATC code to name of corresponding class label
	code_to_label = cbind(names(atc_labels), as.vector(atc_labels)) # combine code and label into single phenotype ID
	new_ids = apply(code_to_label, 1, paste0, collapse=": ")
	return(new_ids)
}
### GET ATC CLASS FROM CODE
# (anly level)
atc_class = function(code, all.levels = F) {
	code.loc = apply(drug$atc[,paste0("code",1:5)], 2, function(x) x %in% code)
	arr.ind = which(code.loc, arr.ind = T)
	level = unique(arr.ind[,2])
	atc.sub = drug$atc[arr.ind[,1],]
	if (all.levels) {
		return(list(unique(atc.sub$class1),
					unique(atc.sub$class2),
					unique(atc.sub$class3),
					unique(atc.sub$class3),
					unique(atc.sub$class4),
					unique(atc.sub$class5)))
	} else {
		return(unique(atc.sub[[paste0("class",level)]]))
	}
}

### Format table for writing
format_atc_out = function(dat, class = NULL) { # dat = sig[[paste0(code,": ",label)]][["4"]]
	if (is.null(class)) { class = c(4,5) }
	
	dat = head(dat, 10)
	dat = cbind(
		#sapply(1:5, function(x) sapply(lapply(dat$atc, atc_class, T), "[[", x)),
		dat,
		sapply(lapply(dat$atc, atc_class, T), "[[", class[1]),
		sapply(lapply(dat$atc, atc_class, T), "[[", class[2]))
	colnames(dat) = c("code","odds-ratio","p","class","subclass")
	return(dat)
}
### Get output path/filename
outname = function(prefix, phen, v) {
	paste0(base.dir,"REALMENT/",prefix,"_",v$outname,"_",phen,".txt")
}
### Write
write_atc_out = function(dat, phen, v) {
	write.table(dat, outname("atc-enrich", phen, v), row.names = F, quote=F, sep='\t')
}

#### END ####





#####################################################
#######		 	UKB MEDICATIONS OVERLAP 		#####
#####################################################
# use this to pre-select drugs for evaluation
# use drugs with higher N (filter in cluster)

## putting this in function to prevent the code from running accidentally

process_ukb_medications = function () {
	atc = drug$atc
	##### UKB medication coding data #####
	meds = fread(paste0(base.dir,"ukb_medications_coding.tsv"), data.table=F)
	colnames(meds)[2] = "drug"
	
	## Add N
	meds_n = fread(paste0(base.dir,"medications_N.csv"), data.table=F) # (this data set was just copied from the data showcase and pasted into excel)
	meds$N = meds_n$Count[match(meds$drug, meds_n$Category)]
	meds = subset(meds, N > 0) # remove 0 count drugs
	meds = meds[order(meds$N, decreasing = T),] # order by N
	head(meds, 20)
	
	## Harmonise drug IDs
	meds$global_id = harmonise_ids_vec(meds$drug)
	head(meds, 20)
	
	# check which ones failed to convert (but might still be relevant due to large N)
	head(meds[is.na(meds$global_id), c("drug","N")], 40)
	
	# removing ' product' (see below)
	meds$drug = gsub(" product","", meds$drug)
	meds$global_id = harmonise_ids_vec(meds$drug)
	
	# renaming levothyroxine (see below)
	grep("thyroxine", meds$drug, value = T) # note: dextrothyroxine sodium is not relevant
	meds$global_id[meds$drug %in% c("levothyroxine sodium","thyroxine sodium", "sodium thyroxine", "thyroxine")] = "levothyroxine"
	
	# adding insulin to the global id so it can still be included for analysis
	meds$global_id[meds$drug == "insulin"] = "insulin"
	
	# drug bank is too specific...
	
	##### END #####
	
	
	
	#######################
	#### ADD ATC CODES ####
	#######################
	# and check N per code
	
	#if (!exists("atc")) {
	#	atc = fread(paste0(base.dir,"atc_jw.txt"), data.table = F)
	#	atc = harmonise_ids(atc, "class5")
	#}
	
	# add N from med list to atc
	# atc.ukb = subset(drug$atc, code3 %in% colnames(atc.extract))
	atc.ukb = drug$atc
	
	# setting insulin to 'insulin human' so it is included in analysis
	subset(atc.ukb, code3 == "A10A")
	atc.ukb$global_id = gsub("insulin human", "insulin", atc.ukb$global_id) 
	
	atc.ukb$ukb_n = meds$N[match(atc.ukb$global_id, meds$global_id)]
	subset(atc.ukb, global_id == "insulin")
	
	# N per class
	atc.n = aggregate(ukb_n ~ class3, atc.ukb, sum)
	atc.n$code3 = atc$code3[match(atc.n$class3, atc$class3)]
	atc.n$code3 = atc$code3[match(atc.n$class3, atc$class3)]
	atc.n$code2 = atc$code2[match(atc.n$class3, atc$class3)]
	atc.n2 = aggregate(ukb_n ~ code2, atc.n, sum)
	
	# ---- check N ----
	# N for stimulants
	subset(atc.ukb, code3 == "N06B")
	
	## CHECK N FOR SSRIs vs OTHER ANTI DEP
	sum(subset(atc.ukb, code4 == "N06AA")$ukb_n, na.rm = T)
	sum(subset(atc.ukb, code4 == "N06AB")$ukb_n, na.rm = T)
	subset(atc.ukb, code4 == "N06AX")
	
	## Antipsychotics
	sum(subset(atc.ukb, code4 == "N05AA")$ukb_n, na.rm = T)
	sum(subset(atc.ukb, code4 == "N05AB")$ukb_n, na.rm = T)
	sum(subset(atc.ukb, code4 == "N05AX")$ukb_n, na.rm = T)
	
	# Anxioulytics
	sum(subset(atc.ukb, code4 == "N05BA")$ukb_n, na.rm = T)
	sum(subset(atc.ukb, code4 == "N05BB")$ukb_n, na.rm = T)
	
	# Anxiolytics
	sum(subset(atc.ukb, code4 == "N03B")$ukb_n, na.rm = T)
	# ---- ---- ----
	
	# check that important drugs in the top drug list are included in these classes
	# e.g. simvastatin
	
	# Extract drugs for all relevant classes (ordered by N)
	atc.list = lapply(unique(atc.ukb$class3), function(x) subset(atc.ukb, class3 == x)[,c("global_id","ukb_n")])
	#atc.list = lapply(unique(atc.drugs$class), function(x) subset(atc.drugs, class == x)[,c("global_id","ukb_n")])
	atc.list = lapply(atc.list, function(x) {
		y = unique(x[!is.na(x[,2]),])
		y = y[order(y[,2], decreasing = T),]
		return(y)
	})
	names(atc.list) = unique(atc.ukb$class3)
	#atc.list = atc.list[-1] # remove drugs this is the drugs that had NA class 
	
	# Map the list of drugs back to the drug names in the ukb medication list
	atc.meds = lapply(atc.list, function(x) { subset(meds, global_id %in% x[,1])[,-5] }) # map the drugs from each class to UKB medication IDs
	names(atc.meds) = atc$code3[match(names(atc.meds), stringr::str_to_title(atc$class3))]
	atc.meds = do.call(rbind, atc.meds)
	atc.meds$atc_code1 = gsub("\\..*","",rownames(atc.meds))
	subset(atc.meds, atc_code1 == "N05B")
	dim(atc.meds)
	## TODO: WHY ARE THERE DUPLICATES???
	# It is from the global ID mapping to multiple ukb IDs
	
	## create more general codes to group similar drugs together
	atc.meds$atc_code2 = substr(atc.meds$atc_code1, 1, 3)
	atc.meds$atc_code3 = substr(atc.meds$atc_code1, 1, 2)
	write.table(atc.meds, paste0(base.dir, "medications_atc.txt"), row.names = F, quote = F, sep = '\t')
	
	# also check level 4 and 5
	atc = drug$atc
	lvl4 = unique(atc[match(atc.meds$global_id, atc$global_id),][,c("code4","global_id")])
	lvl5 = unique(atc[match(atc.meds$global_id, atc$global_id),][,c("code5","global_id")])
	
	lvl4 = cbind(lvl4, meds[match(lvl4$global_id, meds$global_id),][c("drug","N")])
	lvl5 = cbind(lvl5, meds[match(lvl5$global_id, meds$global_id),][c("drug","N")])
	
	l4n = aggregate(N ~ code4, lvl4, sum)
	l4n[order(l4n$N),]
	
	head(atc.meds)
	
	subset(lvl4, code4=="A01AD")
	subset(lvl4, code4=="A10BH")
	subset(lvl4, code4=="C09AA")
	subset(lvl4, code4=="C07FB")
	
	#l5n = aggregate(N ~ code5, lvl5, sum)
	
	lvl4$code1 = substr(lvl4$code4, 1,1)
	subset(lvl4, code1 == "N")
	l4n_N = rbind(l4n[substr(l4n$code4, 1,4) == "N05A",],# antipsychotics
		l4n[substr(l4n$code4, 1,4) == "N05B",], # anxiolytics
		l4n[substr(l4n$code4, 1,4) == "N06A",], # antidepressants
		l4n[substr(l4n$code4, 1,4) == "N03A",]) # anticonvulsants
	
	lvl4_select = subset(l4n_N, N > 200)
	lvl4_N = subset(lvl4, code4 %in% lvl4_select$code4)
	lvl4_N$coding = meds$coding[match(lvl4_N$global_id, meds$global_id)]
	head(lvl4_N)
	write.table(lvl4_N, paste0(base.dir, "medications_atc_lvl4_Neuro.txt"), row.names = F, quote = F, sep = '\t')
	
	
	
	# ATC classes/codes to extract
	atc.extract = data.frame(
		# N NERVOUS SYSTEM
		N03A = "ANTIEPILEPTICS",
		N05A = "ANTIPSYCHOTICS",
		N05B = "ANXIOLYTICS",
		N04B = "DOPAMINERGIC AGENTS",
		N06A = "ANTIDEPRESSANTS",
		#N06C = "PSYCHOLEPTICS AND PSYCHOANALEPTICS IN COMBINATION",
		# N06B = "PSYCHOSTIMULANTS, AGENTS USED FOR ADHD AND NOOTROPICS", # too smal N
		#N06D = "ANTI-DEMENTIA DRUGS", # small N
		
		# C CARDIOVASCULAR SYSTEM
		# C01 CARDIAC THERAPY
		C01B = "ANTIARRHYTHMICS, CLASS I AND III",
		C01D = "VASODILATORS USED IN CARDIAC DISEASES",
		# C07 BETA BLOCKING AGENTS
		C07A = "BETA BLOCKING AGENTS",
		C07B = "BETA BLOCKING AGENTS AND THIAZIDES",
		C07C = "BETA BLOCKING AGENTS AND OTHER DIURETICS",
		C07D = "BETA BLOCKING AGENTS, THIAZIDES AND OTHER DIURETICS",
		C07E = "BETA BLOCKING AGENTS AND VASODILATORS",
		C07F = "BETA BLOCKING AGENTS, OTHER COMBINATIONS",
		# C08 CALCIUM CHANNEL BLOCKERS
		C08C = "SELECTIVE CALCIUM CHANNEL BLOCKERS WITH MAINLY VASCULAR EFFECTS",
		C08D = "SELECTIVE CALCIUM CHANNEL BLOCKERS WITH DIRECT CARDIAC EFFECTS",
		C08G = "CALCIUM CHANNEL BLOCKERS AND DIURETICS",
		# C10 LIPID MODIFYING AGENTS
		C10A = "LIPID MODIFYING AGENTS, PLAIN",
		C10B = " LIPID MODIFYING AGENTS, COMBINATIONS",
		
		# A10 DRUGS USED IN DIABETES
		A10A = "INSULINS AND ANALOGUES",
		A10B = "BLOOD GLUCOSE LOWERING DRUGS, EXCL. INSULINS",
		A10X = "OTHER DRUGS USED IN DIABETES",
		
		# H03 THYROID THERAPY
		H03A = "THYROID PREPARATIONS",
		#H03B = "ANTITHYROID PREPARATIONS", # too small N
		H03C = "IODINE THERAPY",
		
		M01A = "ANTIINFLAMMATORY AND ANTIRHEUMATIC PRODUCTS, NON-STEROIDS"
	)
}


##### END #####








################################
######### RANK RESULTS #########
################################
#d = e; col = "class"

# subset(d[[1]], class == "Antibiotics")
# sort(table(d[[1]]$class))

rank_results = function(d, col) {
	all = sig = data.frame(var = d[[1]][[col]])
	
	for (phen in names(d)) {
		tmp = data.frame(var = d[[phen]][[col]], 
					   rank = nrow(d[[phen]]):1,
						sig = d[[phen]]$p < .05)
		
		# rank all 
		all[[phen]] = tmp$rank[match(all$var, tmp$var)]
		
		# only sig
		tmp$rank[!tmp$sig] = NA # setting non sig to zero and storing ranks as before
		sig[[phen]] = tmp$rank[match(sig$var, tmp$var)]
	}
	
	all = unique(all)
	sig = unique(sig)
	
	all$mean = apply(all[,-1], 1, mean)
	all = all[order(all$mean, decreasing = T),]
	
	sig$mean = apply(sig[,-1], 1, mean, na.rm=T)
	sig = sig[order(sig$mean, decreasing = T),]
	head(sig, 20)
	
	return(list(sig = sig, all = all))
}


#### END ####




#### MAKE NEW LAPPLY FUNCTION ####
# NOTE: Its not really working as I dont know how to add to package environment

# assign to package environment below
mclapply.noreturn = function (X, FUN, ..., mc.preschedule = TRUE, mc.set.seed = TRUE, 
		  mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L), 
		  mc.cleanup = TRUE, mc.allow.recursive = TRUE, affinity.list = NULL) 
{
	cores <- as.integer(mc.cores)
	if ((is.na(cores) || cores < 1L) && is.null(affinity.list)) 
		stop("'mc.cores' must be >= 1")
	.check_ncores(cores)
	if (isChild() && !isTRUE(mc.allow.recursive)) 
		return(lapply(X = X, FUN = FUN, ...))
	if (!is.vector(X) || is.object(X)) 
		X <- as.list(X)
	if (!is.null(affinity.list) && length(affinity.list) < length(X)) 
		stop("affinity.list and X must have the same length")
	if (mc.set.seed) 
		mc.reset.stream()
	if (length(X) < 2) {
		old.aff <- mcaffinity()
		mcaffinity(affinity.list[[1]])
		res <- lapply(X = X, FUN = FUN, ...)
		mcaffinity(old.aff)
		return(res)
	}
	if (length(X) < cores) 
		cores <- length(X)
	if (cores < 2L && is.null(affinity.list)) 
		return(lapply(X = X, FUN = FUN, ...))
	jobs <- list()
	prepareCleanup()
	on.exit(cleanup(mc.cleanup))
	if (!mc.preschedule) {
		FUN <- match.fun(FUN)
		if (length(X) <= cores && is.null(affinity.list)) {
			jobs <- lapply(seq_along(X), function(i) mcparallel(FUN(X[[i]], 
																	...), name = names(X)[i], mc.set.seed = mc.set.seed, 
																silent = mc.silent))
			res <- mccollect(jobs)
			if (length(res) == length(X)) 
				names(res) <- names(X)
			has.errors <- sum(sapply(res, inherits, "try-error"))
		}
		else {
			sx <- seq_along(X)
			res <- vector("list", length(sx))
			names(res) <- names(X)
			fin <- rep(FALSE, length(X))
			if (!is.null(affinity.list)) {
				cores <- max(unlist(x = affinity.list, recursive = TRUE))
				d0 <- logical(cores)
				cpu.map <- lapply(sx, function(i) {
					data <- d0
					data[as.vector(affinity.list[[i]])] <- TRUE
					data
				})
				ava <- do.call(rbind, cpu.map)
			}
			else {
				ava <- matrix(TRUE, nrow = length(X), ncol = cores)
			}
			jobid <- integer(cores)
			for (i in 1:cores) {
				jobid[i] <- match(TRUE, ava[, i])
				ava[jobid[i], ] <- FALSE
			}
			if (anyNA(jobid)) {
				unused <- which(is.na(jobid))
				jobid <- jobid[-unused]
				ava <- ava[, -unused, drop = FALSE]
			}
			jobs <- lapply(jobid, function(i) mcparallel(list(FUN(X[[i]], 
																  ...)), mc.set.seed = mc.set.seed, silent = mc.silent, 
														 mc.affinity = affinity.list[[i]]))
			jobsp <- processID(jobs)
			has.errors <- 0L
			delivered.result <- 0L
			while (!all(fin)) {
				s <- selectChildren(jobs[!is.na(jobsp)], -1)
				if (is.null(s)) 
					break
				if (is.integer(s)) 
					for (ch in s) {
						ji <- match(TRUE, jobsp == ch)
						ci <- jobid[ji]
						r <- readChild(ch)
						if (is.raw(r)) {
							child.res <- unserialize(r)
							if (inherits(child.res, "try-error")) 
								has.errors <- has.errors + 1L
							if (is.list(child.res)) 
								child.res <- child.res[[1]]
							if (!is.null(child.res)) 
								res[[ci]] <- child.res
							delivered.result <- delivered.result + 
								1L
						}
						else {
							fin[ci] <- TRUE
							jobsp[ji] <- jobid[ji] <- NA
							if (any(ava)) {
								nexti <- which.max(ava[, ji])
								if (!is.na(nexti)) {
									jobid[ji] <- nexti
									jobs[[ji]] <- mcparallel(list(FUN(X[[nexti]], 
																	  ...)), mc.set.seed = mc.set.seed, 
															 silent = mc.silent, mc.affinity = affinity.list[[nexti]])
									jobsp[ji] <- processID(jobs[[ji]])
									ava[nexti, ] <- FALSE
								}
							}
						}
					}
			}
			
			s <- length(X) - delivered.result
			if (
				s > 0) 
				warning(sprintf(ngettext(
					s, "%d parallel function call did not deliver a result", 
										 "%d parallel function calls did not deliver results"), 
								
								s), domain = NA)
		}
		if (has.errors) 
			warning(gettextf("%d function calls resulted in an error", 
							 has.errors), domain = NA)
		return(res)
	}
	if (!is.null(affinity.list)) 
		warning("'mc.preschedule' must be false if 'affinity.list' is used")
	sindex <- lapply(seq_len(cores), function(i) seq(i, length(X), 
													 by = cores))
	schedule <- lapply(seq_len(cores), function(i) X[seq(i, length(X), 
														 by = cores)])
	ch <- list()
	res <- vector("list", length(X))
	names(res) <- names(X)
	cp <- rep(0L, cores)
	fin <- rep(FALSE, cores)
	dr <- rep(FALSE, cores)
	inner.do <- function(core) {
		S <- schedule[[core]]
		f <- mcfork()
		if (isTRUE(mc.set.seed)) 
			mc.advance.stream()
		if (inherits(f, "masterProcess")) {
			on.exit(mcexit(1L, structure("fatal error in wrapper code", 
										 class = "try-error")))
			if (isTRUE(mc.set.seed)) 
				mc.set.stream()
			if (isTRUE(mc.silent)) 
				closeStdout(TRUE)
			sendMaster(try(lapply(X = S, FUN = FUN, ...), silent = TRUE))
			mcexit(0L)
		}
		jobs[[core]] <<- ch[[core]] <<- f
		cp[core] <<- processID(f)
		NULL
	}
	job.res <- lapply(seq_len(cores), inner.do)
	ac <- cp[cp > 0]
	has.errors <- integer(0)
	while (!all(fin)) {
		s <- selectChildren(ac[!fin], -1)
		if (is.null(s)) 
			break
		if (is.integer(s)) 
			for (ch in s) {
				a <- readChild(ch)
				if (is.integer(a)) {
					core <- which(cp == a)
					fin[core] <- TRUE
				}
				else if (is.raw(a)) {
					core <- which(cp == attr(a, "pid"))
					job.res[[core]] <- ijr <- unserialize(a)
					if (inherits(ijr, "try-error")) 
						has.errors <- c(has.errors, core)
					dr[core] <- TRUE
				}
				else if (is.null(a)) {
					core <- which(cp == ch)
					fin[core] <- TRUE
				}
			}
	}
	for (i in seq_len(cores)) {
		this <- job.res[[i]]
		if (inherits(this, "try-error")) {
			for (j in sindex[[i]]) res[[j]] <- this
		}
		else if (!is.null(this)) 
			res[sindex[[i]]] <- this
	}
	
	s <- cores - sum(dr)
	if (
		s > 0) 
		warning(sprintf(ngettext(
			s, "scheduled core %s did not deliver a result, all values of the job will be affected", 
								 "scheduled cores %s did not deliver results, all values of the jobs will be affected"), 
						paste(which(dr == FALSE), collapse = ", ")), domain = NA)
	if (length(has.errors)) {
		if (length(has.errors) == cores) 
			warning("all scheduled cores encountered errors in user code")
		else warning(sprintf(ngettext(length(has.errors), "scheduled core %s encountered error in user code, all values of the job will be affected", 
									  "scheduled cores %s encountered errors in user code, all values of the jobs will be affected"), 
							 paste(has.errors, collapse = ", ")), domain = NA)
	}
#	res 	# jw
}

# assign (necessary for internal functions to be accessed)
environment(mclapply.noreturn) = environment(mclapply)


