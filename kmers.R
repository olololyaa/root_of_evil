library(ape)
library(kmer)
library(data.table)
library(dplyr)
library(tictoc)


#setwd("working_dir_here")

#'Calculate k-mers frequency in AA FASTA file
#'@param file_path string  
#'@param kmer_length numeric
#'@param label string
#'@param k.len numeric
#'@return dataframe kmers as a colnames
#'@examples
#'\dontrun{
#'calculate_k_mers("/mnt/data/biohack/kmers/input/Alistipes_senegalensis.faa", 7)
#'}
calculate_k_mers <- function(file_path, kmer_length, label = NULL, k.len){
  
  message(paste("Working with", file_path))
  all.obj <- ape::read.FASTA(file_path,
                             type = "AA")
  counts <- kmer::kcount(all.obj, k = k.len, residues = "AA")
  # Converting to Dayhoff(6) compressed alphabet for k > 3
  # Classes: AGPST, C, DENQ, FWY, HKR, ILM
  counts_df <- data.frame(counts)
  df_sum <- data.frame(t(colSums(counts_df)))
  base_name <- basename(file_path)
  rownames(df_sum) <- tools::file_path_sans_ext(base_name)
  df_sum$"label" <- label
  return(df_sum)
}


#'Calculate k-mers frequency in AA FASTA files in a folder
#'@param folder_path string  
#'@param label string
#'@param file_ext string
#'@param k.len numeric
#'@returns a merged dataframe: filenames a rownames, kmers as a colnames
calc_kmers_for_proteoms_in_a_folder <- function(folder_path, label, file_ext = "faa", k.len){
  tictoc::tic(label)
  ext_pattern <- paste0("\\.", file_ext, "$")
  filenames <- list.files(path = folder_path,
                          pattern = ext_pattern,
                          full.names = TRUE)
  
  df <- purrr::map_dfr(filenames, ~calculate_k_mers(.x, label = label, k.len = k.len))
  tictoc::toc()
  return(df)
}

# healthy 

healthy_df_7 <- calc_kmers_for_proteoms_in_a_folder(folder_path = "/mnt/data/biohack/tax_to_genomes/converted_genomes/MH",
                                                    label = "health",
                                                    file_ext = "faa",
                                                    k.len = 7)
data.table::fwrite(healthy_df_7, file = "kmer_cols_counts_MH_7.csv", row.names=TRUE)
# health: 170.322 sec elapsed

# k > 11 : ! attempt to make a table with >= 2^31 elements
# k = 11 R terminated

# neutral
neutral_df_7 <- calc_kmers_for_proteoms_in_a_folder(folder_path = "/mnt/data/biohack/tax_to_genomes/converted_genomes/neutrals",
                                                    label = "neutral",
                                                    file_ext = "faa",
                                                    k.len = 7)
data.table::fwrite(neutral_df_7, file = "kmer_cols_counts_neutrals_7.csv", row.names=TRUE)



# disease
disease_df_7 <- calc_kmers_for_proteoms_in_a_folder(folder_path = "/mnt/data/biohack/tax_to_genomes/converted_genomes/MN",
                                                    label = "disease",
                                                    file_ext = "faa",
                                                    k.len = 7)
data.table::fwrite(disease_df_7, file = "kmer_cols_counts_disease_7.csv", row.names=TRUE)


#data.table::fread("datatable.csv")


merged_df <- dplyr::bind_rows(disease_df_7, healthy_df_7, neutral_df_7) 
merged_df_mod <- merged_df %>% mutate(status = ifelse(label =="disease", "disease","neutrals"))
data.table::fwrite(merged_df_mod, file = "kmer7_cols_counts_merged_status.csv", row.names=TRUE)