# OBJECTS


# METHODS

setGeneric("revComp", function(seed) standardGeneric("revComp"))
setMethod("revComp", signature("character"),
	function(seed) {
		revSeed <- c()
		for (nt in seed) {
			if (nt == "U") {
				revSeed <- c("A", revSeed)
			}
			else if (nt == "A") {
				revSeed <- c("T", revSeed)
			}
			else if (nt == "C") {
				revSeed <- c("G", revSeed)
			}
			else if (nt == "G") {
				revSeed <- c("C", revSeed)
			}
		}
		revSeed
})

setGeneric("parseMIR", function(file) standardGeneric("parseMIR"))
setMethod("parseMIR", signature("character"),
	function (file) {
		mature <- readLines(file)
		mature <- mature[!grepl(">", mature)]
		seeds <- c()
		for (item in mature) {
			item <- strsplit(item, "")
			item <- item[[1]]
			seed <- item[2:8]
			seed <- revComp(seed)
			seeds <- c(seeds, paste(seed, collapse=""))
			seed <- c(seed[2:length(seed)], "A")
			seeds <- c(seeds, paste(seed, collapse=""))
		}
		seeds
})

setGeneric("anabel", function(patterns) standardGeneric("anabel"))
setMethod("anabel", signature("character"),
	function(patterns) {
		seeds <- c()
		c <- 1
		for (seed in patterns) {
			seed <- strsplit(seed, "")
			seed <- seed[[1]]
			if (c %% 2 == 0) {
				idx <- 1:(length(seed)-1)
				idx <- sample(idx)
				seed <- paste(paste(seed[idx], collapse=""), "A", sep="")
				seeds <- c(seeds, seed)
			}
			else {
				idx <- 1:length(seed)
				idx <- sample(idx)
				seed <- seed[idx]
				seeds <- c(seeds, paste(seed, collapse=""))
			}
		c <- c + 1
		}
		seeds
})

setGeneric("freq", function(patterns, seqs, width=500) standardGeneric("freq"))
setMethod("freq", signature("character", "DNAStringSet"),
	function(patterns, seqs, width=500) {
		cat("\nSelected width\n")
		print(width)
		pb_min <- 0
		pb_max <- length(seqs)
		pb <- txtProgressBar(min=pb_min, max=pb_max, style=3)
		counts <- c()
		for (i in 1:length(seqs)) {
			seq <- strsplit(toString(seqs[i]), "")[[1]]
			len_seq <- length(seq)
			start <- 1
			if (len_seq < width) {
				pb_min <- pb_min + 1
				setTxtProgressBar(pb, pb_min)
				next
			}
			for (j in seq(start, len_seq, by=width)) {
				#subseq <- toString(subseq(seqs[i], start=j, width=width))
				#subseq <- strsplit(subseq, "")[[1]]
				jw <- j+width
				subseq <- seq[j:jw]
				counter <- 0
				for (k in 1:(length(subseq)-7)) {
					l <- k + 6
					counter <- counter + ifelse(paste(subseq[k:l], collapse="") %in% patterns, 1, 0)
				}
				counts <- c(counts, counter)
			}
			pb_min <- pb_min + 1
			setTxtProgressBar(pb, pb_min)
		}
		cat("\n\n")
		counts
})

setGeneric("miRclump", function(patterns, seqs, counts, pv_cutoff=2, width=500) standardGeneric("miRclump"))
setMethod("miRclump", signature("character", "DNAStringSet", "numeric", "numeric", "numeric"),
	function(patterns, seqs, counts, pv_cutoff=2, width=500) {
		start <- 1
		cat("\nCalculating miR-Clumps\n")
		pb_min <- 0
		pb <- txtProgressBar(min=pb_min, max=length(seqs), style=3)
		results <- c()
		rownames <- c()
		colnames <- c("start", "end", "len", "counts", "log_p_value", "higher")
		names <- names(seqs)
		for (i in 1:length(seqs)) {
			seq <- strsplit(toString(seqs[i]), "")[[1]]
			len_seq <- length(seq)
			if (len_seq < width || width <= 0) {
				pb_min <- pb_min + 1
				setTxtProgressBar(pb, pb_min)
				next
			}
			for (j in seq(start, len_seq, by=width)) {
				#subseq <- toString(subseq(seqs[i], start=j, width=width))
				#subseq <- strsplit(subseq, "")[[1]]
				jw <- j+width
				subseq <- seq[j:jw]
				counter <- 0
				for(k in 1:(length(subseq)-7)) {
					l <- k + 6
					counter <- counter + ifelse(paste(subseq[k:l], collapse="") %in% patterns, 1, 0)
				}
				max_p_value <- length(counts[counts >= max(counts)]) / length(counts)
				p_value <- length(counts[counts >= counter]) / length(counts)
				p_value <- ifelse(p_value == 0, max_p_value, p_value)
				higher_than_max <- ifelse(p_value == 0, TRUE, FALSE)
				log_p_value <- -log10(p_value)
				if (log_p_value >= pv_cutoff) {
					results <- rbind(results, c(j, j+width, len_seq, counter, log_p_value, higher_than_max))
					rownames <- c(rownames, names[i])
				}
			}
			pb_min <- pb_min + 1
			setTxtProgressBar(pb, pb_min) 
		}
		if (!is.null(dim(results))) {
			rownames(results) <- rownames
			colnames(results) <- colnames
		}
		else {
			results <- c("No significant p-values were found")
		}
		results
})

setGeneric("Mapping", function(clumps, gff) standardGeneric("Mapping"))
setMethod("Mapping", signature("matrix", "GenomicRanges"),
	function(clumps, gff) {

})

setGeneric("clump2gff", function(clumps) standardGeneric("clump2gff"))
setMethod("clump2gff", signature("matrix"),
	function(clumps) {

})

