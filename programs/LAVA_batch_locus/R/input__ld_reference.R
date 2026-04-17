#' @include global.R

## alignment helper object, translates ordered alelle pairs into numeric code indicating type (absolute value) and strand (sign)
## stores allele codes for reference to match other input
AlleleAligner = R6::R6Class("AlleleAligner",
	private = list(
		alleles = c("A", "C", "G", "T", "a", "c", "t", "g"), allele.pairs = list(AA=NA, CC=NA, GG=NA, TT=NA, 	AC=1, CA=-1, TG=1, GT=-1,		AG=2, GA=-2, TC=2, CT=-2,		AT=11, TA=11, CG=12, GC=12),
		pair.index = NULL, ## data.frame with columns allele1, allele2, code,
		ref.alignment = NULL,  ## allele code vector for input, to compare against

		## create allele pair index from allele string codes
		map.alleles = function(allele1.str, allele2.str) {
			if (length(allele1.str) != length(allele2.str)) input.error("allele string vectors do not have the same length")
			allele1.index = match(allele1.str, private$alleles)
			allele2.index = match(allele2.str, private$alleles)

			return(private$pair.index$code[allele1.index + (allele2.index-1)*length(private$alleles)])
		},

		compare.pairs = function(allele.pair1, allele.pair2) {
			res = rep(NA, length(allele.pair1))
			res[allele.pair1 == -allele.pair2] = -1
			res[allele.pair1 == allele.pair2] = 1
			return(res)
		}
	),
	public = list(
		initialize = function(allele1.str, allele2.str) {
			private$pair.index = data.frame(allele1=rep(private$alleles, times=length(private$alleles)), allele2=rep(private$alleles, each=length(private$alleles)), stringsAsFactors=F)

			private$pair.index$code = unlist(private$allele.pairs[paste0(toupper(private$pair.index$allele1), toupper(private$pair.index$allele2))])
			private$pair.index$code[private$pair.index$code > 10] = NA	# setting strand ambiguous alleles to NA

			private$ref.alignment = private$map.alleles(allele1.str, allele2.str)
		},

		## return values: NA denotes mismatch (including invalid allele pair codes), -1 denotes matching but flipped, 1 denotes matching and same direction
		compare = function(allele1.str, allele2.str, index=NULL) {
			ref.pairs = if (!is.null(index)) private$ref.alignment[index] else private$ref.alignment
			if (length(allele1.str) != length(ref.pairs)) input.error("inconsistent input to function compare()")
			return(private$compare.pairs(ref.pairs, private$map.alleles(allele1.str, allele2.str)))
		}
	)
)


## wrapper for internal IDs for LD reference data; IDs are a numeric index based on the order SNPs are stored in the reference data
## NB: indices are kept unique and sorted
SNPindex = R6::R6Class("SNPindex",
	private = list(internal.ids = integer(0), hash.code = NULL),
	public = list(
		initialize = function(internal.ids) {
			if (length(internal.ids) > 0 && is.numeric(internal.ids)) {
				private$internal.ids = sort(unique(internal.ids))
				if (!is.integer(private$internal.ids)) private$internal.ids = suppressWarnings(as.integer(private$internal.ids))
				if (any(is.na(private$internal.ids)) || private$internal.ids[1] < 0) fatal.error("invalid input to SNPindex object")
			}
			private$hash.code = rlang::hash(private$internal.ids)
		},

		id = function() {return(private$hash.code)},
		snps = function() {return(private$internal.ids)},
		size = function() {return(length(private$internal.ids))},

		equals = function(internal.ids) {
			if (!has.type(internal.ids, "SNPindex")) {
				return(self$size() == length(internal.ids) && all(self$snps() == internal.ids))
			} else return(self$id() == internal.ids$id())
		},

		contains = function(internal.ids) {
			if (has.type(internal.ids, "SNPindex")) internal.ids = internal.ids$snps()
			return(length(internal.ids) <= self$size() && all(internal.ids %in% self$snps()))
		},

		intersection = function(internal.ids) {return(SNPindex$new(private$internal.ids[self$extract.index(internal.ids)]))},

		## return new SNPindex with all IDs within specified range (inclusive)
		subset = function(from, to) {return(SNPindex$new(private$internal.ids[private$internal.ids >= from & private$internal.ids <= to]))},

		## returns numeric index of where this SNPindex matches input, to be used to subset (objects aligned with) internal.ids to SNPs in this SNPindex, in the same order as in this SNPindex
		extract.index = function(internal.ids) {
			if (has.type(internal.ids, "SNPindex")) out = match(self$snps(), internal.ids$snps())
			else out = match(self$snps(), internal.ids)
			return(out[!is.na(out)])
		}
	)
)



#########################################



plink.interface = function(prefix, chromosomes="all") {
	if (grepl("{CHR}", prefix, fixed=T) || grepl("{BLOCK}", prefix, fixed=T)) (PlinkMultiData$new(prefix, chromosomes=chromosomes))
	else return(PlinkData$new(prefix))
}



## base class for different LD reference input formats
LDreferenceData = R6::R6Class("LDreferenceData",
	private = list(
		hash.code = NULL,
		aligner = NULL, ## AlleleAligner object
		alignment.index = NULL, ## vector of allele pair code for each SNP
		snp.info = NULL, ## data.frame with columns internal.id, snp.id, chromosome, position, allele1, allele2
		maximum.id = NA,

		## at present, internal IDs match numeric index into snp.info, as codified here
		from.index = function(index) {return(index)},  ## convert snp.info index to internal ID
		to.index = function(internal.ids) {return(if (has.type(internal.ids, "SNPindex")) internal.ids$snps() else internal.ids)},  ## convert internal ID to snp.info index

		create.id = function() {undefined.error("create.id", class(self)[1])},

		set.snps = function(snp.info) {
			if (any(duplicated(snp.info$snp.id))) private$error("data contains duplicate SNP IDs")

			private$snp.info = data.frame(internal.id=1:nrow(snp.info), snp.info)
			private$maximum.id = nrow(snp.info)
			private$aligner = AlleleAligner$new(snp.info$allele1, snp.info$allele2)
			private$hash.code = private$create.id()
		}
	),
	public = list(
		initialize = function() {fatal.error("cannot initialize base LDreferenceData objects")},

		data.id = function() {return(private$hash.code)},
		equals = function(other) {check.types(other="LDreferenceData"); return(self$data.id() == other$data.id())},

		no.snps = function() {return(nrow(private$snp.info))},
		max.id = function() {return(private$maximum.id)},
		get.snp.info = function(select.ids=NULL, columns=NULL, collapse=T) {
			info = if (!is.null(columns)) private$snp.info[columns] else private$snp.info
			if (!is.null(select.ids)) info = info[private$to.index(select.ids),,drop=F]
			return(if (collapse && ncol(info) == 1) info[,1] else info)
		},

		## sum.stats: data.frame with snp.id, allele1 and allele2 columns
		## updates sum.stats data.frame:
		## - adds chromosome column (if add.chromosome=T)
		## - adds direction column to multiply with statistics later (if skip.align=F), with a value of -1, 0 or 1; set to NA if SNP cannot be aligned
		## - removes the snp.id, allele1 and allele2 columns, and adds an internal.id column corresponding to rows of the reference data SNP info
		align.sumstats = function(sum.stats, add.chromosome=F, skip.align=F) {
			sum.stats$index = match(tolower(sum.stats$snp.id), private$snp.info$snp.id)
			matched = !is.na(sum.stats$index)

			if (!skip.align) {
				sum.stats$direction = NA
				if (any(matched)) sum.stats$direction[matched] = private$aligner$compare(sum.stats$allele1[matched], sum.stats$allele2[matched], index=sum.stats$index[matched])
			}

			if (add.chromosome) {
				sum.stats$chromosome = NA
				sum.stats$chromosome[matched] = private$snp.info$chromosome[sum.stats$index[matched]]
			}

			sum.stats = sum.stats[c("index", names(sum.stats)[!(names(sum.stats) %in% c("index", "snp.id", "allele1", "allele2"))])]
			sum.stats$index = private$from.index(sum.stats$index); names(sum.stats)[1] = "internal.id"
			return(add.rownames(sum.stats, NULL))
		},

		get.block = function(subset) {undefined.error("get.block", class(self)[1])},

		get.chromosomes = function(mode=c("loaded", "available")) {return(validate.chromosomes("all"))},
		set.chromosomes = function(chromosomes) {return(self)}  ## return a copy, unless no change
	)
)



## interface to PLINK data file set; loads .bim file on initialization, SNP IDs are stored in lower case
## call load.data to obtain raw genotype data for specified selection of SNPs
PlinkData = R6::R6Class("PlinkData",
	inherit = LDreferenceData,
	private = list(
		data.prefix = NULL,
		indiv.info = NULL,  ## data.frame with columns indiv.id, family.id, sex, phenotype

		error = function(..., prefix=NULL) {input.error(..., .label=list("PLINK data prefix"=ifelse(is.null(prefix), private$data.prefix, prefix)))},
		data.file = function(prefix, ...) {return(paste0(prefix, ".", c(...)))},

		check.files = function(prefixes) {
			files = unlist(lapply(prefixes, private$data.file, c("bed", "bim", "fam")))
			IF.ANY(missing = !file.exists(files), THEN = private$error("missing %file% for genotype data", .plural=sum(missing), .lines=files[missing]))
		},

		create.id = function() {return(paste(rlang::hash(unlist(private$indiv.info[c("indiv.id", "family.id")])), rlang::hash(private$snp.info$snp.id)))},

		get.printer = function() {return(ObjectPrinter$new("LD Reference (PlinkData)")$add.line("data prefix:", private$data.prefix))},

		read.fam = function(prefix, ids.only=F) {
			indiv = data.table::fread(private$data.file(prefix, "fam"), data.table=F, showProgress=F)
			if (ncol(indiv) != 6) private$error("invalid format for .fam file", prefix=prefix)
			if (nrow(indiv) < globals$ref.minimum.N) private$error("sample size for genotype reference data is smaller than minimum of ", globals$ref.minimum.N)
			indiv = add.names(indiv[,c(2,1,5,6)], c("indiv.id", "family.id", "sex", "phenotype"))

			if (!ids.only) {
				indiv$sex = suppressWarnings(as.numeric(indiv$sex))
				indiv$phenotype = suppressWarnings(as.numeric(indiv$phenotype))

				indiv$sex[!(indiv$sex %in% 1:2)] = NA
				if (all(indiv$phenotype[!is.na(indiv$phenotype)] %in% c(-9,0,1,2))) indiv$phenotype[!(indiv$phenotype %in% 1:2)] = NA  ## set missing values for binary phenotype
			} else indiv = indiv[c("indiv.id", "family.id")]

			return(indiv)
		},

		## internal ID represents index into .bed file
		read.bim = function(prefix) {
			snp.info = data.table::fread(private$data.file(prefix, "bim"), data.table=F, showProgress=F)
			if (ncol(snp.info) != 6) private$error("invalid .bim file, incorrect number of columns", prefix=prefix)

			snp.info = add.names(snp.info[,c(2,1,4:6)], c("snp.id", "chromosome", "position", "allele1", "allele2"))
			snp.info$snp.id = tolower(snp.info$snp.id)

			if (is.character(snp.info$chromosome)) {
				snp.info$chromosome[snp.info$chromosome == "x" | snp.info$chromosome == "X"] = 23
				snp.info$chromosome = suppressWarnings(as.numeric(snp.info$chromosome))
				snp.info$chromosome[is.na(snp.info$chromosome)] = -1
			}

			return(snp.info)
		},

		## assumes subset is a non-empty SNPindex object, validated by load.genotypes
		read.bed = function(subset) {
			invisible(getNamespace("snpStats"))
			out = as(.Call("readbed", private$data.file(private$data.prefix, "bed"), private$indiv.info$indiv.id, self$get.snp.info(subset, "snp.id"), NULL, if (!is.null(subset)) subset$snps(), PACKAGE = "snpStats"), "numeric")
			return(out)
		},

		load.genotypes = function(subset=NULL) {
			if (!is.null(subset)) {
				if (!has.type(subset, "SNPindex")) subset = SNPindex$new(private$from.index(which(private$snp.info$snp.id %in% tolower(subset))))
				if (subset$size() == 0) return(list(data=data.frame(1:nrow(private$indiv.info))[,-1], internal.id=SNPindex$new(integer(0))))
				if (tail(subset$snps(), 1) > private$maximum.id) fatal.error("encountered invalid internal IDs", .label=list("PLINK data prefix"=private$data.prefix))
			}

			out = private$read.bed(subset)
    	if (is.null(subset)) subset = SNPindex$new(private$snp.info$internal.id)
    	return(list(data = add.dimnames(out, NULL), internal.ids = subset))
		}
	),
	public = list(
		initialize = function(prefix)	{
			log.message("reading PLINK reference data set ", prefix)
			private$data.prefix = prefix
			private$check.files(prefix)

			private$indiv.info = private$read.fam(prefix)
			private$set.snps(private$read.bim(prefix))
			log.message("individuals: ", self$no.indiv(), .indent=2)
			log.message("variants: ", self$no.snps(), .indent=2)
		},

		abbreviate = function() {return(paste0("PLINK genotype data (", private$data.prefix, ")"))},
		print = function(...) {
			printer = private$get.printer()
			printer$add.line("data size:", self$no.indiv(), "individuals,", self$no.snps(), "variants")
			cat(printer$to.string())
		},

		no.indiv = function() {return(nrow(private$indiv.info))},
		get.indiv.info = function() {return(private$indiv.info)},

		## subset: list of SNP IDs or a SNPindex object; if NULL, returns all SNPs in data
		## for external SNP IDs, silently ignores any IDs not in data
		## ... should contain any settings to be passed along to LDblock object
		get.block = function(subset, ..., load.individuals=F, decompose=T) {
			data = private$load.genotypes(subset)
			return(LDblockRaw$new(self, data$internal.id, data$data, indiv.info=if (load.individuals) private$indiv.info, decompose=decompose, ...))
		}
	)
)




## multi file set class, split by chromosome and/or block (use {CHR} and/or {BLOCK} token in file prefix specifier)
## treated as a single file set, loads SNP info for all available files
PlinkMultiData = R6::R6Class("PlinkMultiData",
	inherit = PlinkData,
	private = list(
		file.index = NULL,  ## data.frame(prefix, chromosome, block) with file prefixes, and corresponding chromosome and block values (if applicable)
		chromosomes = list(available = NULL, loaded = NULL),

		get.printer = function() {
			printer = ObjectPrinter$new("LD Reference (PlinkMultiData)")$add.line("data prefix:", private$data.prefix)$add.line("no. files:", nrow(private$file.index), indent=2)
			if ("chromosome" %in% names(private$file.index)) printer$add.line("chromosomes:", validate.chromosomes(private$file.index$chromosome)$to.string(), indent=2)
			if ("block" %in% names(private$file.index)) printer$add.line("blocks:", integer.ranges(private$file.index$block, to.string=T), indent=2)
			return(printer)
		},

		resolve.files = function(prefix) {
			mode = list(chromosome=grepl("{CHR}", prefix, fixed=T), block=grepl("{BLOCK}", prefix, fixed=T))
			if (any(unlist(mode))) {
				files = Sys.glob(paste0(gsub("{CHR}", "*", gsub("{BLOCK}", "*", prefix, fixed=T), fixed=T), ".*"))

				matched = extract.tokens(paste0(prefix, ".{EXT}"), files, allowed=c("chr", "block", "ext")); names(matched)[1] = "prefix"

				matched$prefix = substring(matched$prefix, 0, nchar(matched$prefix)-nchar(matched$ext)-1)
				matched = change.names(matched[matched$ext %in% c("bed", "bim", "fam"), names(matched) != "ext", drop=F], chr="chromosome")

				if (mode$block) matched$block = integer.value(matched$block)
				if (mode$chromosome) matched$chromosome = validate.chromosomes(matched$chromosome, condense=F, mode="filter")
				matched = matched[!apply(is.na(matched), 1, any) & !duplicated(matched),,drop=F]

				matched = matched[do.call(order, matched[names(matched) %in% c("chromosome", "block")]),,drop=F]
			} else matched = data.frame(prefix = prefix)

			return(add.rownames(matched, NULL))
		},

		load.data = function(prefixes) {
			private$indiv.info = private$read.fam(prefixes[1])
			if (length(prefixes) > 1) {
				for (pref in prefixes[-1]) {
					fam = private$read.fam(pref, ids.only=T)
					if (nrow(fam) != nrow(private$indiv.info) || !all(fam == private$indiv.info[c("indiv.id", "family.id")])) private$error("individual IDs in file '", private$data.file(pref, "fam"), "' are not consistent with other data files")
				}
			}

			bim = lapply(prefixes, private$read.bim)
			index = data.frame(prefix = private$file.index$prefix, no.snps = sapply(bim, nrow))

			private$file.index = cbind(index, subset(private$file.index, select=-prefix))
			private$set.snps(do.call(rbind, bim))
		},

		## assumes subset is a non-empty SNPindex object, validated by load.genotypes
		read.bed = function(subset) {
			invisible(getNamespace("snpStats"))

			id.range = if (!is.null(subset)) range(subset$snps()) else c(1, self$maximum.id)
			offset = as.integer(c(0, head(cumsum(private$file.index$no.snps), -1)))
			from = offset + 1; to = offset + private$file.index$no.snps
			use = which(from <= id.range[2] & to >= id.range[1])

			out = list()
			for (i in use) {
				bed.file = private$data.file(private$file.index$prefix[i], "bed")
				curr.subset = if (!is.null(subset)) subset$subset(from=from[i], to=to[i])
				out[[i]] = as(.Call("readbed", bed.file, private$indiv.info$indiv.id, self$get.snp.info(curr.subset, "snp.id"), NULL, if (!is.null(curr.subset)) curr.subset$snps() - offset[i], PACKAGE = "snpStats"), "numeric")
			}

			return(do.call(cbind, out))
		}
	),
	public = list(
		initialize = function(prefix, chromosomes="all")	{chromosomes = validate.chromosomes(chromosomes)
			log.message("reading PLINK reference multi-data set ", prefix)
			private$data.prefix = prefix
			private$file.index = private$resolve.files(prefix)

			if (nrow(private$file.index) == 0) private$error("no input files found for specified genotype data prefix")
			if ("chromosome" %in% names(private$file.index)) {
				available = validate.chromosomes(private$file.index$chromosome)
				private$file.index = private$file.index[private$file.index$chromosome %in% chromosomes$get(),]
				if (nrow(private$file.index) == 0) private$error("no input files available for any of requested chromosomes")

				missing = chromosomes$difference(available)
				if (!missing$is.empty() && !(chromosomes$has.all() && missing$equals("X"))) throw.warning("missing requested %chromosome% ", missing$to.string())

				private$chromosomes = list(available = available, loaded = validate.chromosomes(private$file.index$chromosome))
				log.message("chromosomes: ", private$chromosomes$loaded$to.string(), .indent=2)
			} else private$chromosomes = list(available = validate.chromosomes("all"), loaded = validate.chromosomes("all"))
			if ("block" %in% names(private$file.index)) log.message("blocks: ", integer.ranges(private$file.index$block, to.string=T), .indent=2)
			private$check.files(private$file.index$prefix)

			private$load.data(private$file.index$prefix)
			log.message("individuals: ", self$no.indiv(), .indent=2)
			log.message("variants: ", self$no.snps(), .indent=2)
		},

		abbreviate = function() {return(paste0("PLINK genotype multi-data (", private$data.prefix, ")"))},

		no.files = function() {return(nrow(private$file.index))},
		get.file.info = function() {return(private$file.index)},

		get.chromosomes = function(mode=c("loaded", "available")) {return(private$chromosomes[[match.arg(mode)]])},
		set.chromosomes = function(chromosomes) { ## returns a copy, unless no change
			chromosomes = private$chromosomes$available$intersection(chromosomes)
			if (!private$chromosomes$loaded$equals(chromosomes)) return(PlinkMultiData$new(private$data.prefix, chromosomes))
			else return(self)
		}
	)
)



