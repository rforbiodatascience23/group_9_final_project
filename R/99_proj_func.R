make_dir <- function(path){
  if( !dir.exists(path) ){
    dir.create(path)
  }
}

download_dataset <- function(raw_dir, rcall, file_name) {
  
  download.file(rcall$url, str_c(raw_dir, file_name))
  unzip(str_c(raw_dir, file_name), exdir=raw_dir)
  
}

download_dataset_ncbi <- function(raw_dir) {
  
  file_ncbi_name <- "raw_ncbi_data.txt.gz"
  
  if (file.exists(str_c(raw_dir, file_ncbi_name)))
    return()
        
  url_ncbi <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27272/matrix/GSE27272_series_matrix.txt.gz"
  rcall <- httr::GET(url_ncbi)
    
  download_dataset(raw_dir, rcall, file_ncbi_name)
    
  
}

download_data_annotation_ncbi <- function(raw_dir) {
  
  anotation_file_name <- "raw_annotation.bgx.gz"
  file_path <- str_c(raw_dir, anotation_file_name)
  
  url_ncbi <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL6883&format=file&file=GPL6883%5FHumanRef%2D8%5FV3%5F0%5FR0%5F11282963%5FA%2Ebgx%2Egz"
  rcall <- httr::GET(url_ncbi)
  download.file(rcall$url, str_c(raw_dir, anotation_file_name))
  
}

download_dataset_kaggle <- function(raw_dir) {

  file_kaggle_name <- "raw_kaggle_data.zip"
  if (file.exists(str_c(raw_dir, file_kaggle_name)))
    return()
  
  kaggle_base_url <- "https://www.kaggle.com/api/v1"
  
  kaggle_credentials <- '{"username":"silviagoldasova","key":"d562065662daf3359c1c5faad0461172"}'
  user <- fromJSON(kaggle_credentials, flatten = TRUE)
  
  dataset_name <- "rwilliams7653/gse27272-human-placenta-transcriptome?datasetVersionNumber=3"
  
  url_kaggle <- str_c(kaggle_base_url, "/datasets/download/", dataset_name)
  
  rcall <- httr::GET(url_kaggle, httr::authenticate(user$username, user$key, type="basic"))
  
  download_dataset(raw_dir, rcall, file_kaggle_name)
  
}

read_table <- function(raw_dir, file_name) {

  read_data <- read.table(str_c(raw_dir, file_name, ".txt"),sep="\t", header=TRUE)
  return(read_data)
  
}

read_unstructured_ncbi_table <- function(file_path) {
  
  table <- read.table(file_path, header=FALSE, fill = TRUE, col.names = paste0("col_",seq_len(184)))
  return(table) 
    
}

read_bgx_file <- function(raw_dir, file_name) {
  
  outputFilePath <- gunzip(str_c(raw_dir, file_name), remove=FALSE)
  return (readBGX(outputFilePath))
  
}

quantile_normalisation <- function(df){
  df_rank <- apply(df,2,rank,ties.method="min")
  df_sorted <- data.frame(apply(df, 2, sort))
  df_mean <- apply(df_sorted, 1, mean)
  
  index_to_mean <- function(my_index, my_mean){
    return(my_mean[my_index])
  }
  
  df_final <- apply(df_rank, 2, index_to_mean, my_mean=df_mean)
  rownames(df_final) <- rownames(df)
  return(df_final)
}
