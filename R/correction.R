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

