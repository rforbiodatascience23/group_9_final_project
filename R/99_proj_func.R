make_dir <- function(path){
  if( !dir.exists(path) ){
    dir.create(path)
  }
}

download_dataset_ncbi <- function(raw_dir) {
  
  file_ncbi_name <- "raw_ncbi_data.txt.gz"
  
  if (file.exists(str_c(raw_dir, file_ncbi_name)))
    return()
  
  url_ncbi <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE27nnn/GSE27272/matrix/GSE27272_series_matrix.txt.gz"
  
  download.file(url_ncbi, str_c(raw_dir, file_ncbi_name))
  unzip(str_c(raw_dir, file_ncbi_name), exdir=raw_dir)
  
  
}

download_data_annotation_ncbi <- function(raw_dir) {
  
  anotation_file_name <- "raw_annotation.bgx.gz"
  file_path <- str_c(raw_dir, anotation_file_name)
  
  url_ncbi <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GPL6883&format=file&file=GPL6883%5FHumanRef%2D8%5FV3%5F0%5FR0%5F11282963%5FA%2Ebgx%2Egz"
  download.file(url_ncbi, str_c(raw_dir, anotation_file_name))
  
}

read_unstructured_ncbi_table <- function(file_path) {
  
  table <- read.table(file_path, header=FALSE, fill = TRUE, col.names = str_c("col_",seq_len(184)))
  return(table) 
    
}

read_bgx_file <- function(raw_dir, file_name) {
  
  if (file.exists(str_c(raw_dir, "raw_annotation.bgx")))
    return (readBGX(str_c(raw_dir, "raw_annotation.bgx")))
  
  outputFilePath <- gunzip(str_c(raw_dir, file_name), remove=FALSE)
  return (readBGX(outputFilePath))
  
}

get_gene_expr_by_source <- function(gene_expr_df, sample_source_table, source_abr) {
  
  relevant_sample_ids <- sample_source_table |> 
    filter(sample_source == source_abr) |> 
    select(sample_id) |>
    unlist() |> 
    unname()
  
  gene_exp_from_source <- gene_expr_df |> 
    select(Gene, all_of(relevant_sample_ids))
  
  return(gene_exp_from_source)
}

transpose_gene_expr <- function(gene_expr_df) {
  gene_expr_df |> 
    distinct(Gene, .keep_all = TRUE) |>
    pivot_longer(cols = -Gene, names_to = "gene", values_to = "value") |>
    pivot_wider(names_from = "Gene", values_from = "value") |>
    dplyr::rename("geo_accession" = gene)
}

mean_expr <- function(gene_expr_df, term){
  gene_expr_df |>
    filter(smoking_status == "term") |> 
    select(-1:16) |> 
    colMeans()
}

smoking_expr_fc <- function(gene_expr_df){
  mean_expr_smoker = gene_expr_df |>
    filter(smoking_status == "smoker") |> 
    select(-(1:16)) |> 
    colMeans()
  mean_expr_nonsmoker = gene_expr_df |>
    filter(smoking_status == "non-smoker") |> 
    select(-(1:16)) |> 
    colMeans()
  tibble(mean_expr_nonsmoker, 
         mean_expr_smoker) |>  
    mutate(Gene = names(mean_expr_nonsmoker)) |> 
    relocate(Gene) |> 
    mutate(log2fc = log2(mean_expr_smoker) - log2(mean_expr_nonsmoker))
}

create_models <- function(gene_expr_df) {
  gene_expr_df1 <- gene_expr_df |>
    dplyr::select(-names(gene_expr_df)[1],
                  -names(gene_expr_df)[4:16]) |>
    dplyr::mutate(is_smoker = case_when(str_detect(smoking_status, 
                                            "non") == 1 ~ 0,
                                 str_detect(smoking_status, 
                                            "non") == 0 ~ 1),
           .before = 1) |>
    dplyr::select(-smoking_status)
  gene_expr_df1 |>
    tidyr::pivot_longer(cols = colnames(gene_expr_df1[3:length(colnames(gene_expr_df1))]),
                 names_to = "gene",
                 values_to = "expr_level") |>
    dplyr::mutate(log_2_expr_level = log2(expr_level)) |>
    dplyr::select(-expr_level) |>
    drop_na() |>
    dplyr::group_by(gene) |>
    tidyr::nest() |>
    dplyr::ungroup() |>
    # updating data frame. First we group rows by gene, and then we add a variable 
    # containing model fitted to expression level against early_metastasis variable
    dplyr::group_by(gene) |>
    dplyr::mutate(model_object = map(.x = data,
                              .f = ~lm(formula = log_2_expr_level ~ is_smoker,
                                       data = .x))) |>
    # using tidy function to turn model summary into a tibble instead of a list  
    dplyr::mutate(model_object_tidy = map(.x = model_object,
                                   .f = ~ tidy(.x,
                                               conf.int = TRUE,
                                               conf.level = 0.95))) |>
    # unnesting the tidy models, filtering to only get "slope"-coefficients for models and selecting desired variables
    tidyr::unnest(model_object_tidy) |>
    dplyr::filter(term == "is_smoker") |>
    dplyr::select(gene, 
                  p.value, 
                  estimate, 
                  conf.low, 
                  conf.high) |>
    dplyr::ungroup() |>
    # creating variables "q-value" which is the p-value adjusted for multiple comparisons 
    # and "is_significant" which indicates whether observations are statistically 
    # significant when using a level of significance of 0.05
    dplyr::mutate(is_significant = case_when(
      p.value < 0.05 ~ "yes",
      0.05 < p.value ~ "no"
    ))
}
