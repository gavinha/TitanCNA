wigToRangedData <- function(wigfile, verbose = TRUE) {
  if (verbose) { message(paste("Slurping:", wigfile)) }
  input <- readLines(wigfile, warn = FALSE)
  breaks <- c(grep("fixedStep", input), length(input) + 1)
  temp <- NULL
  span <- NULL

  for (i in 1:(length(breaks) - 1)) {
    data_range <- (breaks[i] + 1):(breaks[i + 1] - 1)
    track_info <- input[breaks[i]]
    if (verbose) { message(paste("Parsing:", track_info)) }
    tokens <- strsplit(
      sub("fixedStep chrom=(\\S+) start=(\\d+) step=(\\d+) span=(\\d+)",
      "\\1 \\2 \\3 \\4", track_info, perl = TRUE), " ")[[1]]
    span <- as.integer(tokens[4])
    chr <- rep.int(tokens[1], length(data_range))
    pos <- seq(from = as.integer(tokens[2]), by = as.integer(tokens[3]),
      length.out = length(data_range))
    val <- as.numeric(input[data_range])
    temp <- c(temp, list(data.frame(chr, pos, val)))
  }
  if (verbose) { message("Sorting by decreasing chromosome size") }
  lengths <- as.integer(lapply(temp, nrow))
  temp <- temp[order(lengths, decreasing = TRUE)]
  temp = do.call("rbind", temp)
  output <- RangedData(ranges = IRanges(start = temp$pos, width = span),
    space = temp$chr, value = temp$val)
  return(output)
}

wigToArray <- function(wigfile, verbose = TRUE) {
  if (verbose) { message(paste("Slurping:", wigfile)) }
  input <- readLines(wigfile, warn = FALSE)
  breaks <- c(grep("fixedStep", input), length(input) + 1)
  temp <- NULL
  for (i in 1:(length(breaks) - 1)) {
    data_range <- (breaks[i] + 1):(breaks[i + 1] - 1)
    track_info <- input[breaks[i]]
    if(verbose) { message(paste("Parsing:", track_info)) }
    tokens = strsplit(
      sub("fixedStep chrom=(\\S+) start=(\\d+) step=(\\d+) span=(\\d+)",
      "\\1 \\2 \\3 \\4", track_info, perl = TRUE), " ")[[1]]
    val <- as.numeric(input[data_range])
    temp <- c(temp, list(val))
  }
  if (verbose) { message("Sorting by decreasing chromosome size") }
  lengths <- as.integer(lapply(temp, length))
  temp <- temp[order(lengths, decreasing = TRUE)]
  output <- unlist(temp)
  return(output)
}

rangedDataToWig <- function(correctOutput, file, column = "copy", sample = "R",
    verbose = TRUE) {
  dat <- eval(parse(text = paste("correctOutput", column, sep = "$")))
  if (length(dat) == 0) {
    stop(paste(column, "is not a valid column"))
  }
  dat[is.na(dat)] = -1

  cat(paste("track type=wiggle_0 name=\"", sample, "\"", sep = ""),
    file = file, sep = "\n")
  temp <- data.frame(chr = space(correctOutput), dat)
  width <- start(correctOutput)[2] - start(correctOutput)[1]
  chrs <- levels(space(correctOutput))

  for (i in 1:length(chrs)) {
    chr <- chrs[i]
    out <- temp$dat[temp$chr == chr]
    if (verbose) {
      message(paste("Outputting chromosome ", chr,
        " (", length(out), ")", sep = ""))
    }
    cat(paste("fixedStep chrom=", chr, " start=1 step=", width,
      " span=", width, sep = ""), file = file, append = TRUE, sep = "\n")
    cat(out, file = file, sep = "\n", append = TRUE)
  }
}

rangedDataToSeg <- function(correctOutput, file, column = "copy", sample = "R",
    verbose = TRUE) {
  dat <- eval(parse(text = paste("correctOutput", column, sep = "$")))
  if (length(dat) == 0) {
    stop(paste(column, "is not a valid column"))
  }
  dat[is.na(dat)] = -1

  width <- start(correctOutput)[2] - start(correctOutput)[1]
  out <- data.frame(sample = sample, chr = space(correctOutput), start = start(correctOutput) - 1,
    end = start(correctOutput) + width - 1, value = round(dat, digits = 6))

  write.table(format(out, format = "f", trim = TRUE, drop0trailing = TRUE),
    file = file, quote = FALSE, row.names = FALSE, sep = "\t")
}

wigsToRangedData <- function(readfile, gcfile, mapfile, verbose = FALSE) {
  output <- wigToRangedData(readfile, verbose)
  colnames(output) <- c("reads")
  output$reads <- as.integer(output$reads)
  gc <- wigToArray(gcfile, verbose)

  if (nrow(output) != length(gc)) {
    stop(paste("Number of readcount bins (", nrow(output),
      ") differs from GC count bins (", length(gc), ")", sep = ""));
  }
  map = wigToArray(mapfile, verbose)
  if (nrow(output) != length(map)) {
    stop(paste("Number of readcount bins (", nrow(output),
      ") differs from mappability bins (", length(map), ")", sep = ""));
  }
  output$gc <- gc;
  output$map <- map;
  return(output)
}

correctReadcount <- function(x, mappability = 0.9, samplesize = 50000,
    verbose = TRUE) {
  if (length(x$reads) == 0 | length(x$gc) == 0 | length(x$map) == 0) {
    stop("Missing one of required columns: reads, gc, map")
  }

  if(verbose) { message("Applying filter on data...") }
  x$valid <- TRUE
  x$valid[x$reads <= 0 | x$gc < 0] <- FALSE
  x$ideal <- TRUE
  routlier <- 0.01
  range <- quantile(x$reads[x$valid], prob = c(0, 1 - routlier), na.rm = TRUE)
  doutlier <- 0.001
  domain <- quantile(x$gc[x$valid], prob = c(doutlier, 1 - doutlier),
    na.rm = TRUE)
  x$ideal[!x$valid | x$map < mappability | x$reads <= range[1] |
    x$reads > range[2] | x$gc < domain[1] | x$gc > domain[2]] <- FALSE

  if (verbose) { message("Correcting for GC bias...") }
  set <- which(x$ideal)
  select <- sample(set, min(length(set), samplesize))
  rough = loess(x$reads[select] ~ x$gc[select], span = 0.03)
  i <- seq(0, 1, by = 0.001)
  final = loess(predict(rough, i) ~ i, span = 0.3)
  x$cor.gc <- x$reads / predict(final, x$gc)

  if (verbose) { message("Correcting for mappability bias...") }
  coutlier <- 0.01
  range <- quantile(x$cor.gc[which(x$valid)],
    prob = c(0, 1 - coutlier), na.rm = TRUE)
  set <- which(x$cor.gc < range[2])
  select <- sample(set, min(length(set), samplesize))
  final = approxfun(lowess(x$map[select], x$cor.gc[select]))
  x$cor.map <- x$cor.gc / final(x$map)
  x$copy <- x$cor.map
  x$copy[x$copy <= 0] = NA
  x$copy <- log(x$copy, 2)
  return(x)
}

plotBias <- function(correctOutput, points = 10000, ...) {
  par(mfrow = c(2, 2))
  set <- which(correctOutput$ideal)
  select <- sample(set, min(length(set), points))
  plot(correctOutput$gc[select], correctOutput$reads[select],
    col = densCols(correctOutput$gc[select], correctOutput$reads[select]),
     ylab = "Uncorrected Readcount", xlab = "GC content",
    main = "GC Bias in Uncorrected Readcounts", ...)

  coutlier = 0.001
  range <- quantile(correctOutput$cor.gc[correctOutput$ideal],
    prob = c(0, 1 - coutlier), na.rm = TRUE)
  valid <- which(correctOutput$cor.gc >= range[1] &
    correctOutput$cor.gc <= range[2])
  select <- intersect(valid, select)
  plot(correctOutput$gc[select], correctOutput$cor.gc[select],
    col = densCols(correctOutput$gc[select], correctOutput$cor.gc[select]),
    ylab = "Semi-corrected Readcount", xlab = "GC content",
    main = "GC Bias in Corrected Readcounts", ...)

  select <- sample(valid, min(length(valid), points))
  plot(correctOutput$map[select], correctOutput$cor.gc[select],
    col = densCols(correctOutput$map[select], correctOutput$cor.gc[select]),
    ylab = "Semi-corrected Readcount", xlab = "Mappability",
    main = "Mappability Bias in GC-Corrected Readcounts", ...)

  coutlier = 0.01
  range <- quantile(correctOutput$cor.map, prob = c(0, 1 - coutlier),
    na.rm = TRUE)
  valid <- which(correctOutput$cor.map >= range[1] &
    correctOutput$cor.map <= range[2])
  select <- intersect(valid, select)
  plot(correctOutput$map[select], correctOutput$cor.map[select],
    col = densCols(correctOutput$map[select], correctOutput$cor.map[select]),
    ylab = "Corrected Readcount", xlab = "Mappability",
    main = "GC and Mappability Corrected Readcount", ...)
}

plotCorrection <- function(correctOutput, chr = space(correctOutput)[1], ...) {
  if (!(chr %in% levels(space(correctOutput)))) {
    stop(paste("Invalid chromosome, try one of:",
      paste(levels(space(correctOutput)), collapse = " ")))
  }
  par(mfrow = c(3, 1))
  correctOutput <- correctOutput[paste(chr)]
  pos <- start(correctOutput)
  from <- min(pos)
  to <- max(pos)

  copy <- correctOutput$reads / median(correctOutput$reads, na.rm = TRUE)
  top <- quantile(copy, 0.99)
  bot <- quantile(copy, 0.01)

  set <- which(correctOutput$valid & pos >= from & pos <= to &
    copy >= bot & copy <= top)

  y <- copy[set]
  m <- signif(mad(y, na.rm = TRUE), digits = 3)
  r <- c(min(y, na.rm = TRUE), max(y, na.rm = TRUE))
  plot(pos[set], y, xlab = paste("Position on Chromosome", chr),
    ylab = "Estimated Copy",
    main = paste("Uncorrected Readcount, MAD = ", m), ylim = r, ...)
  m <- signif(mad(correctOutput$cor.gc[set], na.rm = TRUE), digits = 3)
  plot(pos[set], correctOutput$cor.gc[set],
    xlab = paste("Position on Chromosome", chr), ylab = "Estimated Copy",
    main = paste("CG-corrected Readcount, MAD = ", m), ylim = r, ...)
  m <- signif(mad(correctOutput$cor.map[set], na.rm = TRUE), digits = 3)
  plot(pos[set], correctOutput$cor.map[set],
    xlab = paste("Position on Chromosome", chr), ylab = "Estimated Copy",
    main = paste("Mappability and GC-corrected Readcount, MAD = ", m),
    ylim = r, ...)
}
