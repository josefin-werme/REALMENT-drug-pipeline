## input.info can be either:
## - a (vector of) filename(s)
##   - input.dir can be a vector of same length, to specify different directories for each file (use "" or "." for current directory)
## - a data.frame with the requisite columns
## - a list of PhenotypeInfo objects (output from read.info() / create.info())

#' @export
process.input = function(input.info, sample.overlap.file, ref.prefix, phenotypes=NULL, chromosomes="all", input.dir=NULL, defer.loading=F, snp.duplicates=c("drop", "first", "last"), trim.data=F) {
	log.message("Processing input", .indent=0)

	if (is.character(input.info)) { ## treat as (vector of) filenames
		if (!is.null(input.dir)) {
			input.dir[input.dir == ""] = "."
			if (length(input.dir) == 1) input.dir = rep(input.dir, length(input.info))
			else if (length(input.dir) != length(input.info)) input.error("number of elements of input.dir argument is inconsistent with number of file names in input.info")
		} else input.dir = rep(".", length(input.info))
		check.files.exist(input.info)

		pheno.info = list()
		for (i in seq_along(input.info)) {
			log.message("reading phenotype info file ", input.info[i])
			pheno.info = c(pheno.info, validate.info(input.info[i], phenotypes=phenotypes, input.dir=input.dir[i], skip.unknown=(length(input.info) > 1)))
		}
	} else if (is.data.frame(input.info)) {
		log.message("loading phenotype info from input data.frame")
		if (length(input.dir) > 1) input.error("input.dir argument contains multiple elements, with data.frame value for input.info")
		pheno.info = validate.info(input.info, phenotypes=phenotypes, input.dir=input.dir)
	} else if (is.list(input.info)) { ## treat as output from create.info
		if (!all(sapply(input.info, has.type, "RawPhenotype"))) input.error("not all elements in input.info list are RawPhenotype information objects")
		pheno.info = input.info
	} else input.error("argument 'input.info' is of invalid type")

	names(pheno.info) = sapply(pheno.info, function(pi) {pi$name()})
	if (any(duplicated(names(pheno.info)))) input.error("duplicate phenotypes in input")
	if (!is.null(phenotypes)) pheno.info = pheno.info[check.phenotypes(phenotypes, names(pheno.info))]

	sampling.covariance = if (!is.null(sample.overlap.file)) read.sampling.covariance(sample.overlap.file, phenotypes=names(pheno.info), trim=F)

	data = DataSet$new(ref.prefix, snp.duplicates=match.arg(snp.duplicates), trim.data=trim.data, chromosomes=chromosomes)
	for (i in 1:length(pheno.info)) data$add.phenotype(pheno.info[[i]])
	if (!is.null(sampling.covariance)) data$set.sampling.correlation(sampling.covariance)

	if (!defer.loading) {
		data$load()
		if (length(pheno.info) > 1) check.overlap(data)
	}

	return(data)
}


## phenotype can be a named list of phenotype=filename pairs, in which case the filename argument can be omitted
create.info = function(phenotype, filename, ..., input.dir=NULL) {
	if (is.list(phenotype) && !is.null(names(phenotype))) {
		filename = unlist(phenotype)
		phenotype = names(phenotype)
	}

	input.info = list(phenotype=phenotype, filename=filename, ...)
	lengths = sapply(input.info, length)
	if (!all(lengths == 1 | lengths == max(lengths))) input.error("lengths of input are inconsistent")
	return(validate.info(data.frame(input.info), input.dir=input.dir))
}


## input.info can be either a filename, or a data.frame
read.info = function(input.info, phenotypes=NULL, input.dir=NULL) {
	if (length(input.dir) > 1) input.error("input.dir argument contains multiple elements")
	if (input.dir %in% c("", ".")) input.dir = NULL

	if (is.character(input.info) && length(input.info) > 1) input.error("input.info argument contains multiple string elements")
	if (!(is.character(input.info) || is.data.frame(input.info))) input.error("argument input.info has invalid type")

	return(validate.info(input.info, phenotypes=phenotypes, input.dir=input.dir))
}


read.sampling.covariance = function(covar.file, phenotypes=NULL, trim=T) {
	log.message("reading sampling covariance file ", covar.file)
	check.files.exist(covar.file)
	covar = check.covariance(as.matrix(read.table(covar.file, header=T, check.names=F)), allow.NA=T, type="sampling covariance matrix")

	if (!is.null(phenotypes)) {
		exist = phenotypes %in% covar$get.names()
		if (!any(exist)) input.error("none of the input phenotypes are listed in sampling covariance file; please verify that the correct file has been provided")
		if (!all(exist)) {
			throw.warning("%phenotype% ", items.and=phenotypes[!exist], " %is% not listed in sampling covariance file, assuming independent samples")
			covar = covar$expand(phenotypes[!exist], covar.value=0)
		}
		if (trim) covar = covar$subset(phenotypes)
	}

	return(covar)
}


read.plink.data = function(ref.prefix) {return(plink.interface(ref.prefix))}


create.locus = function(..., unit.type=NULL) {
	annot = create.loci(..., unit.type=unit.type)
	if (annot$size() > 1) throw.warning("input contains more than one locus, using the first")
	else if (annot$size() == 0) input.error("input contains no loci")
	return(annot$get.locus(1))
}

create.loci = function(..., unit.type=NULL) {args = list(...)
	if (length(args) == 1 && is.list(args[[1]])) args = args[[1]]
	input = ListInput$new(args, globals$annot.header.index)

	if (!input$has.parameters("locus")) input.error("no 'locus' argument has been provided provided")
	if (!input$has.parameters("snps")) {
		use.param = c("locus", "chromosome", "start", "stop")
		if (!input$has.parameters(use.param[-1])) input.error("either 'snps' or a combination of ", items.and=use.param[-1], " parameters is required")
	} else use.param = c("locus", "snps")

	data = input$get.subset(use.param)
	if (!all(sapply(data[-1], length) == length(data[[1]]))) input.error("input arguments ", items.and=names(data), " have different lengths")

	if (use.param[2] == "snps") return(ListAnnotation$new(data, type=unit.type))
	else return(RangeAnnotation$new(as.data.frame(data), type=unit.type))
}


read.loci = function(loc.file, unit.type=NULL) {
	check.files.exist(loc.file)
	input = TextFileInput$new(loc.file, globals$annot.header.index, input.type="locus")

	if (input$has.parameters("snps")) return(ListAnnotation$new(input$get.subset("locus", "snps"), source=loc.file, type=unit.type))
	else return(RangeAnnotation$new(input$get.subset("locus", "chromosome", "start", "stop"), source=loc.file, type=unit.type))
}


##############################



## validates input data.frame or file and creates corresponding PhenotypeInfo objects
## if a file, will also read global parameters (if any) from file header
## subsets to specified selection of phenotypes if specified (allows partial name matching); disregards unknown phenotypes if skip.unknown=T
## argument global.override is an optional named list of parameters to set to a global value (overrides those specified in file header)
validate.info = function(input, phenotypes=NULL, input.dir=NULL, skip.unknown=F, global.override=NULL) {
	if (!is.null(input.dir)) check.dirs.exist(input.dir)

	## load in and subset data.frame with info
	if (is.data.frame(input)) input.data = DataframeInput$new(input, globals$info.header.index, input.type="info")
	else input.data = TextFileInput$new(input, globals$info.header.index, input.type="phenotype info")
	input.info = input.data$get.subset("phenotype", "filename", optional=names(globals$info.header.index))  ## select parameters and return as data.frame
	input.info[input.info == "." | input.info == ""] = NA   ## treat . and empty as missing value code,


	## process global parameters, add them as constant columns to inupt
	global.params = list.merge(if ("get.globals" %in% names(input.data)) input.data$get.globals() else list(), global.override, drop.null=T)
	if (length(global.params) > 0) {
		params = ListInput$new(global.params, globals$info.header.index, use.mode="last")$get.data()

		existing = names(params) %in% names(input.info)
		if (any(existing)) throw.warning("both global %value% and input %column% %has% been specified for %parameter% ", items.and=names(map)[existing], ", global %value% will be ignored")

		params = params[!existing]
		for (param in c("unit", "parameters", "properties")) {if (param %in% names(params)) params[[param]] = paste(params[[param]], collapse=";")}

		IF.ANY(multiple=sapply(params, length) > 1, THEN=input.error("global %value% specified for %parameter% ", items.and=names(params)[multiple], " %[+each]% %has% multiple elements"))
		input.info = list.merge(input.info, params, drop.null=T)
	}


	## validate phenotype names and input files
	if (any(duplicated(input.info$phenotype))) input.error("duplicate phenotypes in input")
	if (any(grepl("*", input.info$phenotype, fixed=T))) input.error("character '*' is not allowed in phenotype names")
	if (any(grepl("::", input.info$phenotype, fixed=T))) input.error("character string '::' is not allowed in phenotype names")
	if (!is.null(phenotypes)) {
		phenotypes = check.phenotypes(phenotypes, input.info$phenotype, discard.unknown=skip.unknown)
		input.info = input.info[match(phenotypes, input.info$phenotype),]
	}
	if (nrow(input.info) == 0) return(list())
	if (!is.null(input.dir) && input.dir != ".") input.info$filename = paste0(input.dir, "/", input.info$filename)
	check.files.exist(input.info$filename, resolve.chr=T)


	## determine continuous vs binary, and validate associated parameters
	input.types = c("continuous", "binary")
	if ("input.type" %in% names(input.info)) {
		match = pmatch(input.info$input.type, input.types, duplicates.ok=T)
		IF.ANY(invalid=unique(input.info$input.type[!is.na(input.info$input.type) & is.na(match)]), THEN=input.error("unknown input.type %value% ", items.and=invalid))
		input.info$input.type = input.types[match]
	}

	binary.param = c("no.cases", "no.controls", "case.proportion", "prevalence")
	binary.param = binary.param[binary.param %in% names(input.info)]
	for (param in c("sample.size", binary.param)) {if (param %in% names(input.info)) input.info[[param]] = suppressWarnings(as.numeric(input.info[[param]]))}

	if (length(binary.param) > 0) {
		for (param in binary.param) input.info[[param]][!is.na(input.info[[param]]) & input.info[[param]] <= 0] = NA
		is.binary = apply(!is.na(input.info[binary.param]), 1, any)

		if ("input.type" %in% names(input.info)) input.info$input.type[is.na(input.info$input.type)] = input.types[1+is.binary[is.na(input.info$input.type)]]
		else input.info$input.type = input.types[1+is.binary]
	}

	if (!("input.type" %in% names(input.info))) input.error("unable to determine binary/continuous status of any phenotype, please specify")
	IF.ANY(unknown=is.na(input.info$input.type), THEN=input.error("unable to determine binary/continuous status of %phenotype% ", items.and=input.info$phenotype[unknown]))


	## create PhenotypeInfo objects
	pheno.info = list()
	for (i in 1:nrow(input.info)) pheno.info[[input.info$phenotype[i]]] = RawPhenotype$new(input.info[i,])

	return(pheno.info)
}


## helper function to check level of overlap in SNPs between loaded phenotypes in DataSet object
check.overlap = function(data) {
	check.types(data="DataSet")

	pattern.matrix = NULL
	if (data$no.pheno() > 1) {
		ld.reference = data$get.reference(); sum.stats = data$get.sumstats()

		pattern.matrix = add.dimnames(matrix(F, nrow=ld.reference$max.id(), ncol=length(sum.stats)), col.names=names(sum.stats), row.names=1:ld.reference$max.id())
		for (i in 1:length(sum.stats)) pattern.matrix[sum.stats[[i]]$get.ids()$internal.id,i] = T
		if (!data$get.chromosomes()$has.all()) pattern.matrix = pattern.matrix[ld.reference$get.snp.info()$chromosome %in% data$get.chromosomes()$get(),]
		pattern.matrix = pattern.matrix[column.any(pattern.matrix),]

		variants.used = nrow(pattern.matrix)
		overlap.global = sum(column.all(pattern.matrix))
		overlap.proportion = overlap.global/variants.used
		pairwise.overlap = t(pattern.matrix) %*% pattern.matrix
		log.message("data contains ", variants.used, " variants aligned for at least one phenotype, of which ", overlap.global, " variants shared across all phenotypes (", round(100*overlap.proportion, 2), "%)")

		prop.thresh = globals$variant.overlap.warning.proportion
		if (ncol(pattern.matrix) > 2) {
			pairwise.proportion = pairwise.overlap / (outer(diag(pairwise.overlap), diag(pairwise.overlap), FUN="+") - pairwise.overlap)

			proportion.range = range(pairwise.proportion[lower.tri(pairwise.proportion)])
			log.message("pairwise overlap proportions are between ", round(100*proportion.range[1], 2), "% and ", round(100*proportion.range[2],2), "% ", .indent=2)
			if (proportion.range[1] < prop.thresh) throw.warning(ifelse(proportion.range[2] < prop.thresh, "all", "some"), " phenotype pairs have variant overlap proportion of less than ", 100*prop.thresh, "%; analysis may be unstable")
		} else {
			if (overlap.proportion < prop.thresh) throw.warning("variant overlap proportion is less than ", 100*prop.thresh, "%; analysis may be unstable")
		}
		if (max(pairwise.overlap[lower.tri(pairwise.overlap)]) <= globals$variant.overlap.failure.count) input.error("variant overlap between phenotypes is negligible")
	}

	invisible(pattern.matrix)
}





