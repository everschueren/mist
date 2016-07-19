# some qc to make sure config file exists and is well formatted
readMistConfig <- function(config_file) {
    if (!file.exists(config_file)) {
        stop(cat(sprintf("MIST PIPELINE ERROR: YAML file %s does not exist\n", config_file)), call. = F)
    }
    
#     x = readLines(config_file)
#     x = x[x != ""]  #remove \n\n cases (blank Lines)
#     x = gsub(":  ", ": ", gsub(":", ": ", x))  # make sure there is a space between ':' and any character
#     x = gsub("\t", "  ", x)
#     config = paste(x, collapse = "\n")
    config <- tryCatch(yaml::yaml.load_file(config_file), error = function(e) {
        print("MIST PIPELINE ERROR: YAML config file could not be loaded.")
        stop()
    })
    return(config)
}

runMistPipeline <- function(config, controls) {
    ## create an outputdir if it doesn't exist
    if (is.null(config$files$output_dir) || config$files$output_dir == "") 
        config$files$output_dir = sprintf("%s/processed/", getwd())
    
    if(!dir.exists(config$files$output_dir)) dir.create(config$files$output_dir, showWarnings = T)
    
    ## main switches between parts of the pipeline
    if (config$preprocess$enabled) {
        cat(">> PREPROCESSING FILES\n")
        matrix_file = preprocess.main(data_file = config$files$data, keys_file = config$files$keys, output_file = paste(config$files$output_dir, 
            "preprocessed.txt", sep = "/"), filter_data = config$preprocess$filter_contaminants, contaminants_file = config$preprocess$contaminants_file, 
            rm_co = config$preprocess$remove_carryover, collapse_file = config$files$collapse, exclusions_file = config$files$specificity_exclusions, 
            remove_file = config$files$remove, id_colname = config$preprocess$id_colname, prey_colname = config$preprocess$prey_colname, 
            pepcount_colname = config$preprocess$pepcount_colname, mw_colname = config$preprocess$mw_colname)
    }
    if (config$qc$enabled) {
        cat(">> QUALITY CONTROL\n")
        if (!config$preprocess$enabled) {
            ## use previous data matrix instead of the one from pre-processing call
            matrix_file = config$qc$matrix_file
        }
        qc.main(matrix_file = matrix_file, font_scale = config$qc$cluster_font_scale, cluster = config$qc$cluster, 
            ip_dists = config$qc$ip_distributions)
    }
    if (config$mist$enabled) {
        cat(">> MIST\n")
        if (!config$preprocess$enabled) {
            ## use previous data matrix instead of the one from pre-processing call
            matrix_file = config$mist$matrix_file
        }
        mist.results = mist.main(matrix_file = matrix_file, weights = config$mist$weights, w_R = config$mist$reproducibility, 
            w_A = config$mist$abundance, w_S = config$mist$specificity, training_file = config$mist$training_file, 
            training_steps = config$mist$training_steps, controls)
        output_file = gsub(".txt", "_MIST.txt", matrix_file)
        write.table(mist.results, output_file, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")
    }
    
    cat(">> SCORING FINISHED!\n")
}
 
