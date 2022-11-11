#!/pdc/vol/R/3.2.0/bin/Rscript

######################
# plot_pb_stats.R
######################
# Runs with as many arguments as there are chains to analyze
# If there is an argument burnin=, it will be used to discard burnin
# By default, discards the first 100 lines.
# The folder where the chain lie can be given with the argument folder=
# The output file can be specified with output
# The program used (so far, pb or nhpb) can be specified with program=
#
# Rscript ~/bin/comparative/plot_pb_stats.R output=pb_default.pdf folder=pb_default burnin=100 program=pb pb_def1 pb_def2 pb_def3
#
## Init
rm(list=ls())
require(RColorBrewer, quietly=TRUE)
require(signal)
colors <- c(brewer.pal(n=9, "Set1"), brewer.pal(n=8, "Set2"))
pastels <- c(brewer.pal(n=9, "Pastel1"), brewer.pal(n=8, "Pastel2"))
## Parses arguments
args <- commandArgs(TRUE)
## Debug
#setwd("~/projects/crenEuk/pb_default")
#args <- c("burnin=1000", "folder=~/projects/crenEuk/pb_default", "pb_def1", "pb_def2", "pb_def3", "pb_def4", "pb_def5", "pb_def6" )
#args <- c("burnin=10", "folder=~/projects/crenEuk/nhpb_catbp", "program=nhpb", "nhpb_catbp1", "nhpb_catbp2", "nhpb_catbp3", "nhpb_catbp4", "nhpb_catbp5", "nhpb_catbp6" )

## Get burnin
burnin_idx <- grep("^burnin=", args)
burnin <- unlist(strsplit(args[burnin_idx], "="))[2]
## Remove it from args
if (length(burnin_idx) > 0){
  args <- args[-burnin_idx]
} else {
  burnin <- 100
}

## Get program
program_idx <- grep("^program=", args)
program <- unlist(strsplit(args[program_idx], "="))[2]
## Remove it from args
if (length(program_idx) > 0){
  args <- args[-program_idx]
} else {
  program <- "pb"
}
if (!program %in% c("pb", "nhpb"))
  stop("Program should be pb or nhpb")  

## Get folder
folder_idx <- grep("^folder=", args)
folder <- unlist(strsplit(args[folder_idx], "="))[2]
## Remove it from args
if (length(folder_idx) > 0){
  args <- args[-folder_idx]
} else {
  folder <- "."
}

## Get output
output_idx <- grep("^output=", args)
output <- unlist(strsplit(args[output_idx], "="))[2]
## Remove it from args
if (length(output_idx) > 0){
  args <- args[-output_idx]
} else {
  output <- "pb_stats.pdf"
}

## Collect data
n_chains <- length(args)
if (program == "pb"){
  # colnames <- c("x", "time", "topo", "loglik", "length", "alpha", "nmode", "statent", "statalpha", "rrent", "rrmean") # if cat-gtr
  colnames <- c("x", "time", "topo", "loglik", "length", "alpha", "nmode", "statent", "statalpha") # if cat-lg
  variables <- colnames[4:7]
} else if (program == "nhpb"){
  colnames <- c("x", "logsamp", "numberBP", "length", "HPrate")
  variables <- colnames[2:5]
}
n_var <- length(variables)
data <- list()
min <- as.list(rep(Inf, length(colnames)))
max <- as.list(rep(-Inf, length(colnames)))
maxdens <- as.list(rep(-Inf, length(colnames)))
names(min) <- colnames
names(max) <- colnames
names(maxdens) <- colnames
for (i in 1:n_chains){
  if (program == "pb"){
    df <- read.table(paste(folder, "/", args[i], ".trace", sep=""),
                     h=TRUE, comment.char="")
  } else if (program == "nhpb"){
    x <- as.numeric(readLines(paste(folder, "/mon_", args[i], ".",
                                    variables[1], sep="")))
    df <- data.frame(x=1:length(x))
    for (variable in variables){
      df[,variable] <- as.numeric(readLines(paste(folder, "/mon_", args[i], ".",
                                    variable, sep="")))
    }
  }
  df <- df[-(1:burnin),]
  names(df) <- colnames
  for (col in colnames){
    min[col] <- min(min[[col]], min(df[,col]))
    max[col] <- max(max[[col]], max(df[,col]))
    maxdens[col] <- max(maxdens[[col]], max(density(df[,col])$y))
  }
  data[[i]] <- df
}
names(data) <- args
## Savitsky-Golay parameter
sgp <- min(201, round((max$x - min$x)/10)*2-1)
## Plot
jpeg(output, w=8, h=9, res=150, unit="in")
par(mfrow=c(n_var, 2), mar=c(2,4,0,0)+0.1)
for (i in 1:n_var){
  ## Density
  plot(1, 1, type="n", xlim=c(min[[variables[[i]]]], max[[variables[[i]]]]),
       ylim=c(0, maxdens[[variables[[i]]]]),
       ylab=variables[[i]])
  legend("topright", legend=args, col=colors[1:n_chains], fill=colors[1:n_chains],
         cex=0.8)
  for(j in 1:n_chains){
    d <- density(data[[j]][,variables[[i]]])
    lines(d$y ~ d$x, col=colors[j])
  }
  ## xyplot
  plot(1, 1, type="n", xlim=c(min$x, max$x),
       ylim=c(min[[variables[[i]]]], max[[variables[[i]]]]),
       ylab=variables[[i]])
  ## cloud of dots
  for(j in 1:n_chains){
    x <- data[[j]]$x
    y <- data[[j]][,variables[[i]]]
    points(y ~ x, col=pastels[j], cex=0.1)
  }
  # sg lines
  for(j in 1:n_chains){
    x <- data[[j]]$x
    y <- sgolayfilt(data[[j]][,variables[[i]]], p=3, n=sgp)
    ystats <- boxplot(y, plot=FALSE)$stats[,1]
    idx <- y >= ystats[1] & y <= ystats[5] 
    lines(y[idx] ~ x[idx], col=colors[j], lwd=1)
  }
  
}
dev.off()
## x <- data[[j]]$x
## y <- data[[j]][,variables[[i]]]
## sg <- sgolayfilt(y, p=7)
## plot(x, y)
## lines(x, sg)
