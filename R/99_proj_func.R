make_dir <- function(path){
  if( !dir.exists(path) ){
    dir.create(path)
  }
}

download_dataset <- function(raw_dir) {

  kaggle_base_url <- "https://www.kaggle.com/api/v1"
  
  kaggle_credentials <- '{"username":"silviagoldasova","key":"d562065662daf3359c1c5faad0461172"}'
  user <- fromJSON(kaggle_credentials, flatten = TRUE)
  
  dataset_name <- "rwilliams7653/gse27272-human-placenta-transcriptome?datasetVersionNumber=3"
  url <- str_c(kaggle_base_url, "/datasets/download/", dataset_name)
  
  rcall <- httr::GET(url, httr::authenticate(user$username, user$key, type="basic"))
  download.file(rcall$url, str_c(raw_dir, "raw_data.zip"))
  
  unzip(str_c(raw_dir, "raw_data.zip"), exdir=raw_dir)
}

read_tibble <- function(raw_dir, file_name) {

  read_data <- as_tibble(read.table(str_c(raw_dir, file_name, ".txt"),sep="\t", header=TRUE))
  return(read_data)
  
}