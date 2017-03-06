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
  "blast/local/blast_20170222_1214/blasted_20170222_1214_e_all.txt", 
  "blast/genbank/blast_results_nogit.txt"
)

# blast_out_headers <- list(
#   c("qseqid", "sallseqid", "pident", "length", "mismatch", "gapopen", 
#     "qstart", "qend", "sstart", "send", "evalue", "bitscore", "staxids", 
#     "stitle"),
#   c("qseqid", "sallseqid", "pident", )
# )

blast_out_raw <- list()
for(i in 1:length(blast_output_files)){
  blast_out_raw[[i]] <- fread(
    input = file.path(data_dir, blast_output_files[i]), 
    sep = "\t"
    )
}

length(blast_out_raw)
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
  print(msg)
}



temp <- rbind(
  cbind(blast_out_raw[[1]], db = rep("local", nrow(blast_out_raw[[1]]))), 
  cbind(blast_out_raw[[2]], db = rep("genbank", nrow(blast_out_raw[[2]])))
)

# Does the sequence have a match in the FHL course barcode database at >97% identity? (blastn)
unique(temp[ db == "local" & pident >= 97.00 , qseqid])

# Does the sequence match a sequence in GenBank at >97% identity? (blastn)
temp[ db == "genbank" & pident >= 97.00 ]
