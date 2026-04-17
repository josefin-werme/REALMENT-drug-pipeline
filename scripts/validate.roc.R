arg = commandArgs(T); input.dir = arg[1]; in.prefix = arg[2]; code = arg[3]; case_iids = arg[4]; all_iids = arg[5]; n.cores = as.numeric(arg[6])
debug=F
# code="N05BB"; tissue="brain_cortex"; input.dir="drugs"; in.prefix=paste0("maf005-z165-",tissue,"-spearman"); case_iids = paste0("../classN_cases/",code,".cases"); all_iids = "../IIDs.txt"; n.cores = 15; iids_to_process = "1-50000"; outname = paste0("ROC/test",".",iids_to_process); # debug=T

# source('/gpfs/work5/0/vusr0748/realment/scripts/validate.roc.R')

# COMPUTES ROC FOR 'PURE' CASES (i.e. only people that take that particular drug)

library(data.table); library(bettermc)
# - for either cases or controls, extract the p-values for all relevant drugs
# - those from cases will be TPs/FNs, and those from controls are FPs/FNs

### Get drug IDs for input files
sigs = read.table("../signatureIDs_mapped.txt", header=F); colnames(sigs) = "id"
sigs$num = 1:nrow(sigs) # numeric id equivalent to that in the first results based on all ~30k signatures
sig.info = fread("../exemplar_signatures_mappedIDs.txt", data.table=F)
classN.sigs = fread("../classN_signatures.txt",data.table=F,header=F)$V1 # NOTE: subsetting to classN sigs since only these were analysed since whole_blood (cortex and nigra have all)
sigs = sigs[match(classN.sigs, sigs$id),]
infile.drug.ids = sig.info[match(sigs$id, sig.info$compound_id),]$drug # drug ids in same order as input files
num.classN.sigs = sigs$num # get numeric signature ids for whole signature data set corresponding to classN signatures
#rm(sigs, sig.info, classN.sigs)

### If subset to cell type
cell.subtype=T
if (cell.subtype) {
	cellinfo = fread("../cellinfo_beta.txt", header=T)
	celltype = "normal blood sample"
        cells = cellinfo[subtype==celltype]$cell_iname
	cells = c(cells, cellinfo[cell_lineage == "central_nervous_system"]$cell_iname)
	sigs$cell = sig.info[match(sigs$id, sig.info$compound_id),]$cell_iname
	cell.subtype.sigs = sigs$cell %in% cells
	#length(cell.subtype.sigs)
}


### Get target drugs for current ATC class
atc = fread("../atc_processed.txt", data.table=F)
#code.level = paste0("code",nchar(code)-1)
#target.drugs = atc[atc[[code.level]]==code, ]$drug
#target.drugs.idx = infile.drug.ids %in% target.drugs
#rm(atc, code.level, target.drugs)

#### OBTAIN PURE CASES AND OVERLAP BETWEEN CLASSES ####
if (F) {
    for (level in 2:4) {
#	level = 4 # current level
	all.codes = read.table("../classN_phenos.txt") # all code levels
	codes = all.codes[sapply(all.codes, nchar) == level+1] # codes for current level

	cases = list() # read in cases 
	for (i in codes) { cases[[i]] = read.table(paste0("../classN_cases/",i,".cases"))$V1 }
	dupl = table(unlist(cases)) # get duplicate ids (i.e. cases for multiple classes/codes)
	dupl = names(dupl[dupl > 1])

	# prop of dupl cases per code 
	t(t(sort(sapply(codes, function(code) round(sum(cases[[code]] %in% dupl)/length(cases[[code]]),2)))))

	# pure cases
	pure.cases = sapply(codes, function(c) cases[[c]][!cases[[c]] %in% dupl])

	# compare N
	n = data.frame(n.orig = sapply(cases, length), n.pure = sapply(pure.cases, length))
	n$percentage = round(n$n.pure / n$n.orig, 2); n

	# remove classes with 0 cases and save
	pure.cases = pure.cases[sapply(pure.cases, function(x) length(x) > 0)]
	sort(sapply(pure.cases, length))
	for (i in names(pure.cases)) { write.table(pure.cases[[i]], paste0("../classN_cases/",i,".pure"), col.names=F, row.names=F, quote=F) }
    }
}
#### END PURE ####



#### COMPUTE AUC / ENRICHMENT ON INDIVIDUAL LEVEL ####
# i.e. consider true positives as any drug that individual is taking (dont split into cases/controls based on codes)
# need $ukb_phenos/classN_meds_processed_onlycases.dat & $sigs/cellinfo_beta.txt
# tissue="brain_cortex"; input.dir="drugs"; in.prefix=paste0("maf005-z165-",tissue,"-spearman")

# pathways
moa = fread("../compoundIDs_mapped.txt", data.table=F)
moa = moa[match(infile.drug.ids, moa$drug),]
moa = data.frame(drug = moa$drug, moa = tolower(moa$moa))
moa = unique(moa)

terms.excl = c('selective ',' receptor agonist',' receptor antagonist',' receptor modulator',' receptor inverse agonist',
  ' reuptake enhancer',' reuptake inhibitor',' release stimulant',' precursor',' transporter','voltage-gated ',' benzodiazepine site',' uptake',' antiepileptic',' aminotransferase','-norepinephrine',' receptor partial agonist',' benzodiazepine site','voltage-gated ',' channel',
  ' blocker',' inhibitor', ' reputake', ' activator',' agonist')

replace_terms = function(moa, x, new.term, replace) { for (i in replace) { moa[[x]] = gsub(i, new.term, moa[[x]]) }; return(moa) }
moa$broad = moa$moa

x = 'broad'
new.term = ""
replace = terms.excl
moa = replace_terms(moa, x, new.term, replace)
moa$broad = gsub('sert','serotonin',moa$broad)

x = 'cat'
moa[[x]] = moa$broad
new.term = 'catecholamine'
replace = c('adrenergic','norepinephrine','dopamine','catechol o methyltransferase')
moa = replace_terms(moa, x, new.term, replace)
unique(moa[[x]][moa$broad %in% replace])

x = 'cat_mono'
moa[[x]] = moa$cat
new.term = 'monoamine'
replace = c('catecholamine','monoamine oxidase','vesicular monoamine','seroronin','dopamine','tricyclic antidepressant')
moa = replace_terms(moa, x, new.term, replace)
unique(moa[[x]][moa$broad %in% replace])

x = 'cat_mono_exin'
new.term = 'exin'
moa[[x]] = moa$cat_mono
replace = c('potassium','sodium','glutamate','gaba')
moa = replace_terms(moa, x, new.term, replace)
unique(moa[[x]][moa$broad %in% replace])


#sort(unique(moa$broad))
#sort(unique(moa$cat))
#sort(unique(moa$mono.cat))
#sort(unique(moa$mono.cat.exin))

# people taking only clozapine
#cloz.only = casefile[clozapine==1 & n.drugs==1]$IID 
#iid = cloz.only[1]
#iid.list = cloz.only
#iid = "3284762" # N06AA

simple_auc <- function(TPR, FPR){
  # inputs already sorted, best scores first
  dFPR <- c(diff(FPR), 0)
  dTPR <- c(diff(TPR), 0)
  sum(TPR * dFPR) + sum(dTPR * dFPR)/2
}

p.thresh = seq(.01,1,.01)

### Set the number of cores based on number of drugs
#n.cores = round(7000 / length(num.classN.sigs)) # 7200 is based on 500 drugs using almost all memory, so took 450*16 and rounded down
#n.cores = min(16, n.cores) # set to a max of 16
#n.cores = max(1, n.cores) # and min of 1

get_aucs = function(iid.list, MOA, level) {
	print(paste0("** Calculating AUC for MOA = ",MOA)) 
	code.level = paste0("code",level)
	auc = unlist(mclapply(iid.list, mc.cores = n.cores, function(iid) {
		infile = paste0(iid,"-",in.prefix,".drugs")
		d = try(fread(file.path(input.dir, infile), data.table=T)[,-2], silent=T)
		if (class(d)[1]!="try-error") {  # TODO: this is just cause not all outputfiles are finished for nigra... FIX THAT 

			# since the first predictions were done with all ~30k signatures, those files will need to be subsetted to class N signatures
			setkey(d, sig)
			if (nrow(d) > length(num.classN.sigs)) {
				d = d[J(num.classN.sigs), nomatch=0L]
			}
			d$drug = infile.drug.ids

			if (cell.subtype) {
				cell.subtype.sigs = d$sig[cell.subtype.sigs]
				d = d[J(cell.subtype.sigs), nomatch=0L]
			}

			case.drugs = gsub("_"," ",colnames(drugs.by.iid)[drugs.by.iid[IID==iid]==1])
			case.atc = unique(atc[atc$drug %in% case.drugs, ][[code.level]])
			target.drugs = unique(atc[atc[[code.level]] %in% case.atc, ]$drug)

			if (MOA!="no_moa") {
				target.moas = unique(moa[moa$drug %in% target.drugs,][[MOA]])
				target.drugs = unique(moa[moa[[MOA]] %in% target.moas,]$drug) # extract all drugs with shared pathways
			}

			# get TPR and FPR
			rates = do.call(rbind, lapply(p.thresh, function(pt) {
				sig.targ = d[p < pt]$drug %in% target.drugs
				nsig.targ = d[p > pt]$drug %in% target.drugs
				data.frame(
					TPR = sum(sig.targ) / (sum(sig.targ) + sum(nsig.targ)),
					FPR = sum(!sig.targ) / (sum(!sig.targ) + sum(!nsig.targ))
				)
			}))
			simple_auc(rates$TPR, rates$FPR)
		}
	}))
	return(auc)
}

# tissue="brain_substantia_nigra"; input.dir="drugs"; in.prefix=paste0("maf005-z165-",tissue,"-spearman"); case_iids = paste0("../classN_cases/",code,".cases"); all_iids = "../IIDs.txt"; n.cores = 15; iids_to_process = "1-50000"; outname = paste0("ROC/test",".",iids_to_process); # debug=Ti
#tissue="whole_blood"; input.dir="drugs"; in.prefix=paste0("maf005-z165-",tissue,"-spearman"); case_iids = paste0("../classN_cases/",code,".cases"); all_iids = "../IIDs.txt"; n.cores = 15; iids_to_process = "1-50000"; outname = paste0("ROC/test",".",iids_to_process); # debug=T
# source('/gpfs/work5/0/vusr0748/realment/scripts/validate.roc.R')
setwd(paste0("../",in.prefix))
debug=F
clndrugs = read.table("../drugs.subset.txt")$V1
cln = read.table("/gpfs/work5/0/vusr0748/realment/scripts/classN_ROC.txt")$V1
cln = cln[cln %in% unique(subset(atc, drug %in% clndrugs)$code4)]
read.table("../classN_ROC.sub.txt")$V1

for (code in cln) {
print(code)

# code = "N06DA"
level = nchar(code)-1 # for not pure, code can be multidrug_level2 multidrug_level3 multidrug_level4

# drug info for classN cases
drugs.by.iid = fread("../classN_meds_processed_onlycases.dat",header=T)
colnames(drugs.by.iid)[colnames(drugs.by.iid)=="cases"] = "n.drugs"

# extract IIDs for pure cases or people with multiple drugs
if (level < 10) {
	iid.list = read.table(paste0("../classN_cases/",code,".pure"), header=F)$V1
} else {
	level = gsub("notpure_code","",code)
	iid.list = drugs.by.iid[n.drugs > 1]$IID
}
iid=iid.list[2]; i="no_moa" # for running interactively

if (cell.subtype) { moas = "no_moa" } else { moas = c("no_moa", "moa","broad","cat","cat_mono","cat_mono_exin") }
mean.aucs = data.frame(moa = moas, median = NA)

if (!debug) {
	for (i in moas) {
		aucs = get_aucs(iid.list, MOA=i, level); mean(aucs)
		if (cell.subtype) { outname = paste0("ROC/",code,".",i,".aucs_sub") } else { outname = paste0("ROC/",code,".",i,".aucs") }
		fwrite(data.table(IID = iid.list, AUC = signif(aucs, 4)), outname, sep="\t")
		mean.aucs$median[mean.aucs$moa==i] = signif(median(aucs, na.rm=T),4)
	}
}
if (cell.subtype) { outname = paste0("ROC/",code,".",i,".median.aucs_sub") } else { outname = paste0("ROC/",code,".",i,".median.aucs") }
write.table(mean.aucs, outname, row.names=F, quote=F, sep='\t')

}
