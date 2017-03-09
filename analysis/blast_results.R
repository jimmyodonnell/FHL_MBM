# create a table from the results of blastn on local and genbank databases
################################################################################

EXPORT <- FALSE

#
library(data.table) # fread(); data.table()
library(magrittr) # %>%

INTERACTIVE <- FALSE

if(INTERACTIVE){
  analysis_dir <- dirname(file.choose()) # choose this script file
  setwd(analysis_dir)
} else {
  analysis_dir <- getwd()
}
data_dir <- file.path("..", "data")
fig_dir <- file.path("..", "figures")

#-------------------------------------------------------------------------------
# LOAD DATA

# point to file locations (relative to data directory)

# counts table:
counts_table_file <- "7_no_primers.unique.pick.MACSE.precluster.pick.count_table"

# script used to run blast:
blast_script_file <- "blast/local/blast_nested_e-local.sh"

# blast results files
blast_output_files <- c(
  local = "blast/local/blast_20170222_1214/blasted_20170222_1214_e_all.tsv", 
  genbank = "blast/genbank/blast_results_nogit.txt"
)

# blast_out_headers <- list(
#   c("qseqid", "sallseqid", "pident", "length", "mismatch", "gapopen", 
#     "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", 
#     "stitle"),
#   c("qseqid", "sallseqid", "pident", )
# )

blast_out_raw <- list()
for(i in 1:length(blast_output_files)){
  db_name <- names(blast_output_files)[i]
  blast_out_raw[[db_name]] <- fread(
    input = file.path(data_dir, blast_output_files[i]), 
    sep = "\t"
    )
}
length(blast_out_raw)
names(blast_out_raw)
lapply(blast_out_raw, dim)

# counts table:
counts_table <- fread(input = file.path(data_dir, counts_table_file))
dim(counts_table)
counts_table[1:5, 1:5]


# get header information.
blast_script <- readLines(file.path(data_dir, blast_script_file))
target_line <- blast_script[grep(pattern = "output_format=", blast_script)]
target_str <- strsplit(target_line, split = '"')[[1]][2]
blast_out_header <- substr(x = target_str, start = 3, stop = nchar(target_str))
blast_out_header <- strsplit(blast_out_header, split = " ")[[1]]
length(blast_out_header)
#-------------------------------------------------------------------------------

for(i in 1:length(blast_out_raw)){
  names(blast_out_raw[[i]]) <- blast_out_header
}

for(i in 1:length(blast_out_raw)){
  in_blastout <- length(
    counts_table[,Representative_Sequence] %in% blast_out_raw[[i]][,qseqid])
  msg <- paste0(in_blastout, " of ", nrow(counts_table), 
                " sequences from count table in blast output ", i)
  message(msg)
}

#===============================================================================
# How many hits does each query sequence have?
#-------------------------------------------------------------------------------
hit_table <- list()

for(i in 1:length(blast_out_raw)){
  hits_per_qseq <- table(blast_out_raw[[i]][,qseqid])
  
  hit_table[[i]] <- data.frame(
    N_hits = as.numeric(names(table(hits_per_qseq))),
    Freq = as.numeric(table(hits_per_qseq))
  )
}

#-------------------------------------------------------------------------------
# the local results look suspicious
hits_per_qseq_loc <- sort(table(blast_out_raw[["local"]][,qseqid]))
tail(hits_per_qseq_loc, 10)
length(hits_per_qseq_loc)

gt_10 <- which(hits_per_qseq_loc > 10)
pident_trail <- list()
for(i in 1:length(gt_10)){
  pident_trail[[i]] <- sort(blast_out_raw[["local"]][
    qseqid == names(hits_per_qseq_loc)[gt_10[i]] , pident], decreasing = TRUE)
}
trail_len <- sapply(pident_trail, length)

# plot trajectories of percent identity when there are multiple hits
plot_name <- "perc_ident_trajectory"

if(!exists("legend_text")){legend_text <- list()}
legend_text[plot_name] <- {"
Trajectory of percent identity for hits against local database.
Each line represents a given query sequence, and its path denotes the sorted percent identity.
Lines are colored by the total number of hits.
"}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 4) #, width = 8, height = 3
}


library(viridis)
col_lev <- viridis(length(unique(trail_len)), alpha = 0.5)

par_o <- par()
par(bg = "grey50", mar = c(4,4,1,1))
plot(c(1,max(trail_len)), range(unlist(pident_trail)), type = "n", 
     xlab = "hit", ylab = "percent identity", 
     bg = "gray",
     # log = "y",
     las = 1
)
for(i in 1:length(gt_10)){
  points(pident_trail[[i]], 
    col = col_lev[as.numeric(as.factor(trail_len))[i]], 
    lwd = 2, 
    type = "l")
}
dev.off()

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting hits per query
xmax <- max(sapply(hit_table, function(x) max(x$N_hits)))
ymax <- max(sapply(hit_table, function(x) max(x$Freq)))


plot_name <- "hits_per_query"

if(!exists("legend_text")){legend_text <- list()}
legend_text[plot_name] <- {"
Frequency of the number of hits per query sequence. 
Points are colored by the database against which sequences were queried.
"}

if(EXPORT){
  pdf_file    <- file.path(fig_dir, paste(plot_name, ".pdf", sep = ""))
  legend_file <- file.path(fig_dir, paste(plot_name, "_legend.txt", sep = ""))
  writeLines(legend_text[[plot_name]], con = legend_file)
  pdf(file = pdf_file, width = 5, height = 4) #, width = 8, height = 3
}

par(mar = c(4,4,1,1))
plot(
  c(0,xmax), c(0,ymax), 
  xlab = "Number of hits per query sequence", 
  ylab = "Frequency", 
  type = "n", 
  las = 1
  )
points(hit_table[[2]], col = hsv(1, 1, 1))
points(hit_table[[1]])
legend("topright", 
  legend = c("local", "GenBank"), 
  title = "Blast Database", 
  bty = "n", 
  pch = 1, 
  col = c(2,1)
)
if(EXPORT){
  dev.off()
}
#===============================================================================


# Append sequence database column
for(i in 1:length(blast_out_raw)){
  blast_out_raw[[i]][,db := names(blast_out_raw)[i]]
}

# combine into a single table
blast_out <- rbindlist(l = blast_out_raw, use.names = TRUE)


#===============================================================================
# Assemble taxonomic annotation table
#-------------------------------------------------------------------------------

is.duplicated <- function(x) { 
  # this function returns a vector where each element corresponds to
  # whether or not an element of x is found in x more than once.
  x %in% x[duplicated(x)]
}

#-------------------------------------------------------------------------------
# LOCAL DATABASE
#-------------------------------------------------------------------------------
# Does the sequence have a match in the FHL course barcode database at >97% identity? (blastn)
# from all of the results,
blast_out %>%
  # extract hits in local database with 97.0 percent identity
  .[db == "local" & pident >= 97.0, 
    ] %>% 
  .[, # extract the hits with the max pident
    .SD[which(pident == max(pident))], 
    by = qseqid] %>% 
  .[,
    # if there are multiple equally good hits, paste the sallseqid names together
    list(pident, hit_name = paste(sallseqid, collapse = ";"), db), 
    by = qseqid] %>% 
  # eliminate duplicate rows
  unique -> hits_local
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# GENBANK
#-------------------------------------------------------------------------------
# Does the sequence match a sequence in GenBank at >97% identity? (blastn)

blast_out %>%
  # extract hits in genbank database with 97.0 percent identity
  # that did not have a hit in the local database
  .[db == "genbank" & pident >= 97.0 & !(qseqid %in% hits_local$qseqid), 
    ] %>% 
  .[, # extract the hits with the max pident
    .SD[which(pident == max(pident))], 
    by = qseqid] %>% 
  .[,
    # if there are multiple equally good hits, paste the staxids names together
    list(pident, hit_name = paste(unique(staxids), collapse = ";"), db), #unique
    by = qseqid] %>%
  # eliminate duplicate rows
  unique -> hits_genbank

#-------------------------------------------------------------------------------
blast_hits <- rbind(hits_local, hits_genbank)
#-------------------------------------------------------------------------------


#-------------------------------------------------------------------------------
# get taxon names
#-------------------------------------------------------------------------------


hits_genbank$hit_name %>% 
  strsplit

unlist(strsplit(hits_genbank$hit_name[920:940], split = ";"))

tax_ann <- data.table(qseqid = counts_table[,Representative_Sequence])

library(plyr)
join(tax_ann, blast_hits, type = "left")

merge(x = tax_ann, y = hits_local)

dt1 <- data.table(let = letters[1:26])
reordered <- sample(1:20)
dt2 <- data.table(
  let = letters[reordered], 
  num = reordered
)
library(plyr)
join(dt1, dt2, type = "left")
