# create a table from the results of blastn on local and genbank databases
################################################################################

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
  local = "blast/local/blast_20170222_1214/blasted_20170222_1214_e_all.txt", 
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
dim(blast_out_raw[[2]])
blast_out_raw[[2]][ V2 %like% "BMBM"]


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

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plotting hits per query
xmax <- max(sapply(hit_table, function(x) max(x$N_hits)))
ymax <- max(sapply(hit_table, function(x) max(x$Freq)))


plot_name <- "hits_per_query"

EXPORT <- TRUE

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

dim(blast_out)

#===============================================================================
# Assemble taxonomic annotation table
#-------------------------------------------------------------------------------

tax_ann <- data.table(seq_id = counts_table[,Representative_Sequence])

# tax_ann[,fhl_db := seq_id %in% ]

# Does the sequence have a match in the FHL course barcode database at >97% identity? (blastn)
good_hits_per_qseq <- table(blast_out[ db == "local" & pident >= 97.00, qseqid])


tax_ann[, good_hits := as.numeric(good_hits_per_qseq[tax_ann$seq_id])] 

# how many sequences do not have hits in FHL CO1 DB >97
nrow(
  tax_ann[ good_hits = "NA" , ]
  )
tax_ann[ good_hits > 0, ]

# Does the sequence match a sequence in GenBank at >97% identity? (blastn)
temp[ db == "genbank" & pident >= 97.00 ]
#-------------------------------------------------------------------------------
