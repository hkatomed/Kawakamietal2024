
## Before running the script:
##   Store S288C_reference_sequence_R64-1-1_20110203.fsa in reference_genome.
##     The FASTA file can be found at: 
##     http://sgd-archive.yeastgenome.org/sequence/S288C_reference/genome_releases/S288C_reference_genome_R64-1-1_20110203.tgz
##   Store the bedGraph files (.bdg.gz) in bedGraph/FY23.YE.
##     The bedGraph files can be obtained at GEA (https://www.ddbj.nig.ac.jp/gea/index-e.html).
##   Store 13059_2018_1398_MOESM2_ESM.xlsx in XLSX.
##     The excel file can be obtained at:
##     https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5807854/bin/13059_2018_1398_MOESM2_ESM.xlsx
##   There is Combinedpeaks_cpm_annot_with_TranscriptID.txt in RECLU/Gene_Cluster_Expression
##     Relevant files can be obtained at GEA (https://www.ddbj.nig.ac.jp/gea/index-e.html).
##   There are S288C_ENY32_A04_01.ab1, FY23_ENY32_A03_01.ab1 and FY24_ENY32_E03_05.ab1 in AB1.

## The term YE in file names stands for yeast extract dextrose peptone (YPD)
## FY23 is a MATa strain.

## Load the packages.
library(Biostrings)	# https://bioconductor.org/packages/release/bioc/html/Biostrings.html
library(openxlsx)	# https://cran.r-project.org/web/packages/openxlsx/index.html
library(sangerseqR)	# https://www.bioconductor.org/packages/release/bioc/html/sangerseqR.html

## Define the home directory
HOMEDIR <- "HOMEDIR"

## Load the genome, rename the chromosomes and save it.
setwd(HOMEDIR)
setwd("reference_genome")
sc.genome <- readDNAStringSet(filepath = "S288C_reference_sequence_R64-1-1_20110203.fsa")
newNames <- paste("chr", c(1:16, "M"), sep = "")
names(sc.genome) <- newNames
setwd("../")
dir.create("RData")
setwd("RData")
save(sc.genome, file = "sc_genome.RData")

## Unzip the .gz files.
setwd(HOMEDIR)
setwd("bedGraph/FY23.YE")
command <- "gunzip -fk *.gz"
system(command)

## Define a function that loads a bedGraph and generate a data.frame.
getTSSscores <- function(fileName){
	inDF <- read.table(file = fileName, skip = 1)
	inDF$V2 <- NULL
	colnames(inDF) <- c("chr", "pos", "score")

	outList <- list()
	for(c in 1:16){
		chrNum <- as.character(as.roman(c))
		chrName <- paste("chr", c, sep = "")
		chrLen <- width(sc.genome)[c]
		temp <- subset(inDF, chr == chrNum)
		outDF <- data.frame(chr = chrName, pos = 1:chrLen, score = 0)
		cat("chrName: ", chrName, " chrLen: ", chrLen, "\n", sep = "")
		for(t in 1:nrow(temp)){
			outDF$score[outDF$pos == temp$pos[t]] <- temp$score[t]
		}
		outList[[chrName]] <- outDF
	}
	return(outList)
}

## Define the bedGraph filenames to be loaded
fileNames <- c("FY23.YE.1.plus.bdg", "FY23.YE.1.minus.bdg", 
		"FY23.YE.2.plus.bdg", "FY23.YE.2.minus.bdg")

## Load the bedGraph files and store the data.frames in a list named TSS_FY23_YE.
TSS_FY23_YE <- list()
TSS_FY23_YE$plus1 <- getTSSscores(fileName = fileNames[1])
TSS_FY23_YE$minus1 <- getTSSscores(fileName = fileNames[2])
TSS_FY23_YE$plus2 <- getTSSscores(fileName = fileNames[3])
TSS_FY23_YE$minus2 <- getTSSscores(fileName = fileNames[4])

## Combine the data.frame objects 
TSS_FY23_YE$plus <- TSS_FY23_YE$plus1
for(i in 1:16){
TSS_FY23_YE$plus[[i]]$score <- TSS_FY23_YE$plus1[[i]]$score + TSS_FY23_YE$plus2[[i]]$score
}

TSS_FY23_YE$minus <- TSS_FY23_YE$minus1
for(i in 1:16){
TSS_FY23_YE$minus[[i]]$score <- TSS_FY23_YE$minus1[[i]]$score + TSS_FY23_YE$minus2[[i]]$score
}

## Save the list
setwd(HOMEDIR)
setwd("RData")
save(TSS_FY23_YE, file = "TSS_FY23_YE.RData")


## Load Chereji's +1 nucleosome data
setwd(HOMEDIR)
setwd("XLSX")
Chereji <- read.xlsx(xlsxFile = "13059_2018_1398_MOESM2_ESM.xlsx")

## Rearrange the data.frame and change the strand values from 1 and -1 to + and -
Chereji <- Chereji[, c(2, 8, 5)]
for(i in 1:nrow(Chereji)){
	Strand <- Chereji$Strand[i]
	if(Strand == 1) Chereji$strand[i] <- "+"
	if(Strand == -1) Chereji$strand[i] <- "-"
}
Chereji <- Chereji[, c(1, 2, 4)]
names(Chereji) <- c("chr", "pos", "strand")

## Change the chromosome names from chrI to chr1
for(i in 1:nrow(Chereji)){
	temp <- as.integer(as.roman(strsplit(Chereji$chr[i], split = "chr")[[1]][2]))
	Chereji$chr[i] <- paste("chr", temp, sep = "")
}

## Save the data.frame of Chereji's +1 nucleosomes
setwd(HOMEDIR)
setwd("RData")
save(Chereji, file = "Chereji_plus1.RData")


## Extract CTSS counts around the +1 nucleosomes (-200:200)
nucData <- Chereji
outList <- list()
outList$sense <- matrix(nrow = nrow(nucData), ncol = 401)
outList$antisense <- matrix(nrow = nrow(nucData), ncol = 401)
for(n in 1:nrow(nucData)){
	chrName <- nucData$chr[n]
	pos <- nucData$pos[n]
	strand <- nucData$strand[n]
	chrNum <- as.integer(strsplit(chrName, split = "chr")[[1]][2])
	chrLen <- width(sc.genome)[chrNum]
	if(pos > 200 & pos + 199 < chrLen){
		if(strand == "+"){
			outList$sense[n,] <- TSS_FY23_YE$plus1[[chrName]]$score[(pos - 200):(pos + 200)]
			outList$antisense[n,] <- - TSS_FY23_YE$minus1[[chrName]]$score[(pos - 200):(pos + 200)]
		}
		if(strand == "-"){
			outList$sense[n,] <- - TSS_FY23_YE$minus1[[chrName]]$score[(pos + 200):(pos - 200)]
			outList$antisense[n,] <- TSS_FY23_YE$plus1[[chrName]]$score[(pos + 200):(pos - 200)]
		}
	}

}
TSS_P1_FY23_YE <- list()
TSS_P1_FY23_YE$sample1 <- outList

outList <- list()
outList$sense <- matrix(nrow = nrow(nucData), ncol = 401)
outList$antisense <- matrix(nrow = nrow(nucData), ncol = 401)
for(n in 1:nrow(nucData)){
	chrName <- nucData$chr[n]
	pos <- nucData$pos[n]
	strand <- nucData$strand[n]
	chrNum <- as.integer(strsplit(chrName, split = "chr")[[1]][2])
	chrLen <- width(sc.genome)[chrNum]
	if(pos > 200 & pos + 199 < chrLen){
		if(strand == "+"){
			outList$sense[n,] <- TSS_FY23_YE$plus2[[chrName]]$score[(pos - 200):(pos + 200)]
			outList$antisense[n,] <- - TSS_FY23_YE$minus2[[chrName]]$score[(pos - 200):(pos + 200)]
		}
		if(strand == "-"){
			outList$sense[n,] <- - TSS_FY23_YE$minus2[[chrName]]$score[(pos + 200):(pos - 200)]
			outList$antisense[n,] <- TSS_FY23_YE$plus2[[chrName]]$score[(pos + 200):(pos - 200)]
		}
	}

}
TSS_P1_FY23_YE$sample2 <- outList

outList <- list()
outList$sense <- matrix(nrow = nrow(nucData), ncol = 401)
outList$antisense <- matrix(nrow = nrow(nucData), ncol = 401)
for(n in 1:nrow(nucData)){
	chrName <- nucData$chr[n]
	pos <- nucData$pos[n]
	strand <- nucData$strand[n]
	chrNum <- as.integer(strsplit(chrName, split = "chr")[[1]][2])
	chrLen <- width(sc.genome)[chrNum]
	if(pos > 200 & pos + 199 < chrLen){
		if(strand == "+"){
			outList$sense[n,] <- TSS_FY23_YE$plus[[chrName]]$score[(pos - 200):(pos + 200)]
			outList$antisense[n,] <- - TSS_FY23_YE$minus[[chrName]]$score[(pos - 200):(pos + 200)]
		}
		if(strand == "-"){
			outList$sense[n,] <- - TSS_FY23_YE$minus[[chrName]]$score[(pos + 200):(pos - 200)]
			outList$antisense[n,] <- TSS_FY23_YE$plus[[chrName]]$score[(pos + 200):(pos - 200)]
		}
	}

}
TSS_P1_FY23_YE$combined <- outList

## Remove the genes with no CTSS count around the +1 nucleosome.
sumScores <- list()
sumScores$sample1$sense <- rowSums(TSS_P1_FY23_YE$sample1$sense)
sumScores$sample1$antisense <- rowSums(TSS_P1_FY23_YE$sample1$antisense)
sumScores$sample2$sense <- rowSums(TSS_P1_FY23_YE$sample2$sense)
sumScores$sample2$antisense <- rowSums(TSS_P1_FY23_YE$sample2$antisense)
sumScores$combined$sense <- rowSums(TSS_P1_FY23_YE$combined$sense)
sumScores$combined$antisense <- rowSums(TSS_P1_FY23_YE$combined$antisense)

## Get order of genes, higher CTSS counts on top.
sumScores$sample1$senseOrder <- order(sumScores$sample1$sense, decreasing = TRUE)
sumScores$sample1$antisenseOrder <- order(sumScores$sample1$antisense, decreasing = TRUE)
sumScores$sample2$senseOrder <- order(sumScores$sample2$sense, decreasing = TRUE)
sumScores$sample2$antisenseOrder <- order(sumScores$sample2$antisense, decreasing = TRUE)
sumScores$combined$senseOrder <- order(sumScores$combined$sense, decreasing = TRUE)
sumScores$combined$antisenseOrder <- order(sumScores$combined$antisense, decreasing = TRUE)

groupAllsample1 <- sumScores$sample1$senseOrder[1:length(sumScores$sample1$senseOrder)]
groupAllsample2 <- sumScores$sample2$senseOrder[1:length(sumScores$sample2$senseOrder)]
groupAllcombined <- sumScores$combined$senseOrder[1:length(sumScores$combined$senseOrder)]

## Save the data.
setwd(HOMEDIR)
setwd("RData")
save(TSS_P1_FY23_YE, file = "TSS_P1_FY23_YE.RData")
save(groupAllsample1, file = "groupAllsample1.RData")
save(groupAllsample2, file = "groupAllsample2.RData")
save(groupAllcombined, file = "groupAllcombined.RData")


## Get nucleotide sequences
WIDTH <- 200
nucData <- Chereji
outChar <- character(length = nrow(nucData))
names(outChar) <- nucData$name

for(n in 1:nrow(nucData)){
	chrName <- nucData$chr[n]
	pos <- nucData$pos[n]
	strand <- nucData$strand[n]
	chrNum <- as.integer(strsplit(chrName, split = "chr")[[1]][2])
	chrLen <- width(sc.genome)[chrNum]
	if(pos > WIDTH & pos + (WIDTH-1) < chrLen){
		SEQ <- sc.genome[[chrName]][(pos - WIDTH):(pos + WIDTH)]
		if(strand == "-")	SEQ <- reverseComplement(SEQ)
		outChar[n] <- as.character(SEQ)
	}
}

SEQ_P1 <- outChar

## save the data
setwd(HOMEDIR)
setwd("RData")
save(SEQ_P1, file = "SEQ_P1.RData")


################################################################
## plotting functions used in the paper (all the highest CTSS)

## Define a function to plot TSS distribution
plotTSSgroup <- function(group, ylim = NA, SAMPLE = TSS_P1_FY23_YE$sample1, lwd = 3, ABLINE = TRUE){
	temp <- SAMPLE$sense[get(group), ]
	for(i in 1:nrow(temp)){
		maxPos <- which.max(temp[i,])
		maxIndex <- which(temp[i,] == temp[i, maxPos])
		maxScore <- 1/length(maxIndex)
		temp[i,] <- 0
		temp[i, maxIndex] <- maxScore
	}
	if(is.na(sum(ylim)) == TRUE){
		ymax <- max(colSums(temp))
		ylim = c(0, ymax)
	}
	geneNum <- nrow(temp)
	plot(x = -200:200, y = colSums(temp), main = paste(group, "(n=", geneNum, ")", sep = ""),  
		type = "n", ylim = ylim, 
		ylab = "Occurrence of highest CTSS peak", xlab = "Distance from +1 dyad (bp)")
	if(ABLINE == TRUE){
		abline(v = -73, lty = 2, col = "gray", lwd = lwd)
		abline(v = -63, lty = 2, col = "gray", lwd = lwd)
		abline(v = -53, lty = 2, col = "gray", lwd = lwd)
		abline(v = -43, lty = 2, col = "gray", lwd = lwd)
		abline(v = -33, lty = 2, col = "gray", lwd = lwd)
		abline(v = -23, lty = 2, col = "gray", lwd = lwd)
		abline(v = -13, lty = 2, col = "gray", lwd = lwd)
		abline(v = -3, lty = 2, col = "gray", lwd = lwd)		
	}
	box(lwd = lwd)
	axis(1, labels = FALSE, lwd.ticks = lwd)
	axis(2, labels = FALSE, lwd.ticks = lwd)
	axis(1, at = c(-73, 73), labels = FALSE, lwd.ticks = lwd)

	# > col2rgb("lightblue")[,1]
	#   red green  blue 
	#   173   216   230 

	polygon(x = c(-73, 73, 73, -73), y = c(-10, -10, 200, 200), 
		# col = rgb(173, 216, 230, max = 255, alpha = 0.5))
		col = gray(0, alpha = 0.1), border = NA)
	# lines(x = -200:200, y = colMeans(temp), col = "red", lwd = lwd)
	lines(x = -200:200, y = colSums(temp), col = "red", lwd = lwd)
	# return(colSums(temp))

	dyad.pos <- 201
	inside <- sum(colMeans(temp)[(dyad.pos - 73):(dyad.pos + 73)])
	outside <- sum(colMeans(temp)[1:(dyad.pos - 74)])+sum(colMeans(temp)[(dyad.pos + 74):401])
	total <- sum(colMeans(temp))
	cat("total: ", total, "; inside: ", inside, "; outside: ", outside, "\n", sep = "")
	inside.ratio <- inside / total
	outside.ratio <- outside / total
	cat("inside ratio: ", inside.ratio, "; outside ratio: ", outside.ratio, "\n", sep = "")
}


## Figure 3A
ylim <- c(0, 150)
lwd <- 3
plotTSSgroup(group = "groupAllsample1", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$sample1, lwd = lwd, ABLINE = FALSE)
plotTSSgroup(group = "groupAllsample2", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$sample2, lwd = lwd, ABLINE = FALSE)
plotTSSgroup(group = "groupAllcombined", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = lwd, ABLINE = FALSE)
plotTSSgroup(group = "groupAllsample1", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$sample1, lwd = lwd, ABLINE = TRUE)
plotTSSgroup(group = "groupAllsample2", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$sample2, lwd = lwd, ABLINE = TRUE)
plotTSSgroup(group = "groupAllcombined", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = lwd, ABLINE = TRUE)

plotTSSgroup(group = "groupAllcombined", ylim = ylim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = 4, ABLINE = FALSE)

## Mark the position -44, the peak
abline(v = -44)


################################################################
## Plotting functions used in the paper (all the highest CTSS)

## isA at -8
plotTSSisAgroup <- function(group, ylim = NA, xlim = NA, SAMPLE = TSS_P1_FY23_YE$sample1, lwd = 3, ABLINE = TRUE){
	temp <- SAMPLE$sense[get(group), ]
	hit.number <- 0
	for(i in 1:nrow(temp)){

		SEQ <- SEQ_P1[get(group)][i]
		ORDER <- order(temp[i,], decreasing = TRUE)

		hit <- FALSE
		for(o in 1:length(ORDER)){
			maxPos <- ORDER[o]
			maxIndex <- which(temp[i,] == temp[i, maxPos])
			if(length(maxIndex) == 1){
				minus8 <- substr(SEQ, start = maxIndex - 8, stop = maxIndex - 8)

				if(minus8 == "A"){
					temp[i,] <- 0
					maxScore <- 1
					temp[i, maxIndex] <- maxScore
					hit <- TRUE
					break
				}
			}
		}
		if(hit == FALSE)	temp[i,] <- 0
		if(hit == TRUE)		hit.number <- hit.number + 1
	}

	if(is.na(sum(ylim)) == TRUE){
		ymax <- max(colSums(temp))
		ylim = c(0, ymax)
	}
	if(is.na(sum(xlim)) == TRUE){
		xlim = c(-WIDTH, WIDTH)
	}
	geneNum <- hit.number
	plot(x = -WIDTH:WIDTH, y = colSums(temp), main = paste("isA ", group, " (n=", hit.number, ")", sep = ""),  
		type = "n", ylim = ylim, xlim = xlim, 
		ylab = "Occurrence of highest CTSS peak", xlab = "Distance from +1 dyad (bp)")
	if(ABLINE == TRUE){
		abline(v = -73, lty = 2, col = "gray", lwd = lwd)
		abline(v = -63, lty = 2, col = "gray", lwd = lwd)
		abline(v = -53, lty = 2, col = "gray", lwd = lwd)
		abline(v = -43, lty = 2, col = "gray", lwd = lwd)
		abline(v = -33, lty = 2, col = "gray", lwd = lwd)
		abline(v = -23, lty = 2, col = "gray", lwd = lwd)
		abline(v = -13, lty = 2, col = "gray", lwd = lwd)
		abline(v = -3, lty = 2, col = "gray", lwd = lwd)		
	}
	box(lwd = lwd)
	axis(1, labels = FALSE, lwd.ticks = lwd)
	axis(2, labels = FALSE, lwd.ticks = lwd)
	axis(1, at = c(-73, 73), labels = FALSE, lwd.ticks = lwd)

	# > col2rgb("lightblue")[,1]
	#   red green  blue 
	#   173   216   230 

	polygon(x = c(-73, 73, 73, -73), y = c(-10, -10, 200, 200), 
		# col = rgb(173, 216, 230, max = 255, alpha = 0.5))
		col = gray(0, alpha = 0.1), border = NA)
	lines(x = -WIDTH:WIDTH, y = colSums(temp), col = "red", lwd = lwd)
	# return(colSums(temp))
}


## Figure 3C, lower
lwd <- 3
ylim <- c(0, 150)
xlim <- c(-100, 100)
plotTSSisAgroup(group = "groupAllcombined", ylim = ylim, xlim = xlim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = lwd, ABLINE = FALSE)
plotTSSisAgroup(group = "groupAllcombined", ylim = ylim, xlim = xlim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = lwd, ABLINE = TRUE)
plotTSSisAgroup(group = "groupAllcombined", ylim = ylim, xlim = xlim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = 4, ABLINE = TRUE)


## nonA at -8
plotTSSnonAgroup <- function(group, ylim = NA, xlim = NA, SAMPLE = TSS_P1_FY23_YE$sample1, lwd = 3, ABLINE = TRUE){
	temp <- SAMPLE$sense[get(group), ]
	hit.number <- 0
	for(i in 1:nrow(temp)){

		SEQ <- SEQ_P1[get(group)][i]
		ORDER <- order(temp[i,], decreasing = TRUE)

		hit <- FALSE
		for(o in 1:length(ORDER)){
			maxPos <- ORDER[o]
			maxIndex <- which(temp[i,] == temp[i, maxPos])
			if(length(maxIndex) == 1){
				minus8 <- substr(SEQ, start = maxIndex - 8, stop = maxIndex - 8)

				if(minus8 != "A"){
					temp[i,] <- 0
					maxScore <- 1
					temp[i, maxIndex] <- maxScore
					hit <- TRUE
					break
				}
			}
		}
		if(hit == FALSE)	temp[i,] <- 0
		if(hit == TRUE)		hit.number <- hit.number + 1
	}

	if(is.na(sum(ylim)) == TRUE){
		ymax <- max(colSums(temp))
		ylim = c(0, ymax)
	}
	if(is.na(sum(xlim)) == TRUE){
		xlim = c(-WIDTH, WIDTH)
	}
	geneNum <- hit.number
	plot(x = -WIDTH:WIDTH, y = colSums(temp), main = paste("nonA ", group, " (n=", geneNum, ")", sep = ""),  
		type = "n", ylim = ylim, xlim = xlim, 
		ylab = "Occurrence of highest CTSS peak", xlab = "Distance from +1 dyad (bp)")
	if(ABLINE == TRUE){
		abline(v = -73, lty = 2, col = "gray", lwd = lwd)
		abline(v = -63, lty = 2, col = "gray", lwd = lwd)
		abline(v = -53, lty = 2, col = "gray", lwd = lwd)
		abline(v = -43, lty = 2, col = "gray", lwd = lwd)
		abline(v = -33, lty = 2, col = "gray", lwd = lwd)
		abline(v = -23, lty = 2, col = "gray", lwd = lwd)
		abline(v = -13, lty = 2, col = "gray", lwd = lwd)
		abline(v = -3, lty = 2, col = "gray", lwd = lwd)		
	}
	box(lwd = lwd)
	axis(1, labels = FALSE, lwd.ticks = lwd)
	axis(2, labels = FALSE, lwd.ticks = lwd)
	axis(1, at = c(-73, 73), labels = FALSE, lwd.ticks = lwd)

	# > col2rgb("lightblue")[,1]
	#   red green  blue 
	#   173   216   230 

	polygon(x = c(-73, 73, 73, -73), y = c(-10, -10, 200, 200), 
		# col = rgb(173, 216, 230, max = 255, alpha = 0.5))
		col = gray(0, alpha = 0.1), border = NA)
	lines(x = -WIDTH:WIDTH, y = colSums(temp), col = "red", lwd = lwd)
	# return(colSums(temp))
}

## Figure 3C, upper
plotTSSnonAgroup(group = "groupAllcombined", ylim = ylim, xlim = xlim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = lwd, ABLINE = FALSE)
plotTSSnonAgroup(group = "groupAllcombined", ylim = ylim, xlim = xlim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = lwd, ABLINE = TRUE)
plotTSSnonAgroup(group = "groupAllcombined", ylim = ylim, xlim = xlim, SAMPLE = TSS_P1_FY23_YE$combined, lwd = 4, ABLINE = TRUE)


##########################################################

## Define functions to get TSS sequences
seqTSSnonA8group <- function(group, SAMPLE = TSS_P1_FY23_YE$sample1){
	WIDTH <- 200
	temp <- SAMPLE$sense[get(group), ]
	SEQTSS <- character(length = nrow(temp))
	names(SEQTSS) <- names(SEQ_P1[get(group)])
	for(i in 1:nrow(temp)){

		SEQ <- SEQ_P1[get(group)][i]
		ORDER <- order(temp[i,], decreasing = TRUE)

		for(o in 1:length(ORDER)){
			maxPos <- ORDER[o]
			maxIndex <- which(temp[i,] == temp[i, maxPos])
			if(length(maxIndex) == 1){
				minus8 <- substr(SEQ, start = maxIndex - 8, stop = maxIndex - 8)

				if(minus8 != "A"){
					if(maxIndex > 100 & maxIndex < (WIDTH*2+1-100)){
						SEQTSS[i] <- substr(SEQ, start = (maxIndex-100), stop = (maxIndex + 100))
					}else{
						SEQTSS[i] <- ""
					}
					break
				}
			}
		}
	}
	SEQTSS <- SEQTSS[which(SEQTSS != "")]
	return(DNAStringSet(SEQTSS))
}

seqTSSwithA8group <- function(group, SAMPLE = TSS_P1_FY23_YE$sample1){
	WIDTH <- 200
	temp <- SAMPLE$sense[get(group), ]
	SEQTSS <- character(length = nrow(temp))
	names(SEQTSS) <- names(SEQ_P1[get(group)])
	for(i in 1:nrow(temp)){

		SEQ <- SEQ_P1[get(group)][i]
		ORDER <- order(temp[i,], decreasing = TRUE)

		for(o in 1:length(ORDER)){
			maxPos <- ORDER[o]
			maxIndex <- which(temp[i,] == temp[i, maxPos])
			if(length(maxIndex) == 1){
				minus8 <- substr(SEQ, start = maxIndex - 8, stop = maxIndex - 8)

				if(minus8 == "A"){
					if(maxIndex > 100 & maxIndex < (WIDTH*2+1-100)){
						SEQTSS[i] <- substr(SEQ, start = (maxIndex-100), stop = (maxIndex + 100))
					}else{
						SEQTSS[i] <- ""
					}
					break
				}
			}
		}
	}
	SEQTSS <- SEQTSS[which(SEQTSS != "")]
	return(DNAStringSet(SEQTSS))
}

seqTSSanygroup <- function(group, SAMPLE = TSS_P1_FY23_YE$sample1){
	WIDTH <- 200
	temp <- SAMPLE$sense[get(group), ]
	SEQTSS <- character(length = nrow(temp))
	names(SEQTSS) <- names(SEQ_P1[get(group)])
	for(i in 1:nrow(temp)){

		SEQ <- SEQ_P1[get(group)][i]
		ORDER <- order(temp[i,], decreasing = TRUE)
		maxPos <- ORDER[1]
		maxIndex <- which(temp[i,] == temp[i, maxPos])

		if(length(maxIndex) == 1){
			if(maxIndex > 100 & maxIndex < (WIDTH*2+1-100)){
				SEQTSS[i] <- substr(SEQ, start = (maxIndex-100), stop = (maxIndex + 100))
			}else{
				SEQTSS[i] <- ""
			}
		}
	}
	SEQTSS <- SEQTSS[which(SEQTSS != "")]
	return(DNAStringSet(SEQTSS))
}


## Get sequences around TSS
SEQTSS_nonA8_201 <- seqTSSnonA8group(group = "groupAllcombined", SAMPLE = TSS_P1_FY23_YE$combined)
SEQTSS_withA8_201 <- seqTSSwithA8group(group = "groupAllcombined", SAMPLE = TSS_P1_FY23_YE$combined)
SEQTSS_any_201 <- seqTSSanygroup(group = "groupAllcombined", SAMPLE = TSS_P1_FY23_YE$combined)

SEQTSS_nonA8_25 <- subseq(SEQTSS_nonA8_201, start = 101-12, end = 101+12)
SEQTSS_withA8_25 <- subseq(SEQTSS_withA8_201, start = 101-12, end = 101+12)
SEQTSS_any_25 <- subseq(SEQTSS_any_201, start = 101-12, end = 101+12)

## Save the sequence for WebLogo. Figure 3B
setwd(HOMEDIR)
setwd("RData")
save(SEQTSS_nonA8_201, file = "SEQTSS_nonA8_201.RData")
save(SEQTSS_withA8_201, file = "SEQTSS_withA8_201.RData")
save(SEQTSS_any_201, file = "SEQTSS_any_201.RData")
save(SEQTSS_nonA8_25, file = "SEQTSS_nonA8_25.RData")
save(SEQTSS_withA8_25, file = "SEQTSS_withA8_25.RData")
save(SEQTSS_any_25, file = "SEQTSS_any_25.RData")
setwd("../")
dir.create("FASTA")
setwd("FASTA")
writeXStringSet(SEQTSS_nonA8_25, filepath = "SEQTSS_nonA8_25.fasta")
writeXStringSet(SEQTSS_withA8_25, filepath = "SEQTSS_withA8_25.fasta")
writeXStringSet(SEQTSS_any_25, filepath = "SEQTSS_any_25.fasta")
setwd("../")
dir.create("TXT")
setwd("TXT")
write(as.character(SEQTSS_nonA8_25), file = "SEQTSS_nonA8_25.txt")
write(as.character(SEQTSS_withA8_25), file = "SEQTSS_withA8_25.txt")
write(as.character(SEQTSS_any_25), file = "SEQTSS_any_25.txt")


## Load the matrix of CAGE signals. 
setwd(HOMEDIR)
setwd("RECLU/Gene_Cluster_Expression")
cpmDF <- read.table(file = "Combinedpeaks_cpm_annot_with_TranscriptID.txt", header = TRUE, sep = "\t")

## Correlation test
cor.test(x = cpmDF[["FY23.SC_FY23.SC.1"]], y = cpmDF[["FY23.SC_FY23.SC.2"]])
cor.test(x = cpmDF[["FY23.YE_FY23.YE.1"]], y = cpmDF[["FY23.YE_FY23.YE.2"]])
cor.test(x = cpmDF[["FY24.SC_FY24.SC.1"]], y = cpmDF[["FY24.SC_FY24.SC.2"]])
cor.test(x = cpmDF[["FY24.YE_FY24.YE.1"]], y = cpmDF[["FY24.YE_FY24.YE.2"]])


## Draw sequence charts. For Figure 2C.
setwd(HOMEDIR)
setwd("AB1")
S288C <- readsangerseq(filename = "S288C_ENY32_A04_01.ab1")
FY23 <- readsangerseq(filename = "FY23_ENY32_A03_01.ab1")
FY24 <- readsangerseq(filename = "FY24_ENY32_E03_05.ab1")

WIDTH <- 50
START <- 417
chromatogram(S288C, trim5 = START, trim3 = length(primarySeq(S288C)) - START - WIDTH, 
	showcalls = "primary", showhets = FALSE, width = 50, ylim = 1.8)

START <- 420
chromatogram(FY23, trim5 = START, trim3 = length(primarySeq(FY23)) - START - WIDTH, 
	showcalls = "primary", showhets = FALSE, width = 50, ylim = 1.8)

START <- 419
chromatogram(FY24, trim5 = START, trim3 = length(primarySeq(FY24)) - START - WIDTH, 
	showcalls = "primary", showhets = FALSE, width = 50, ylim = 1.8)




 

