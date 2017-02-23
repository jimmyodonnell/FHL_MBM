# create a table from the results of blastn on local and genbank databases
################################################################################

#
library(data.table) # fread(); data.table()
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
blast_script_file <- "blast/blast_nested_e-2.sh"

# blast results files
blast_output_files <- c(
  "blast/genbank/correct_fields_nogit.txt", 
  # "blast/genbank/FHL_blast_GenBank_results.nogit.txt",
  "blast/local/blast_20170222_1214/blasted_20170222_1214_e_all.txt"
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
  
}
temp <- cbind(blast_out_raw[[2]], db = rep("local", nrow(blast_out_raw[[2]])))

temp
