## parameter grid

mist.train.getParamGrid = function(steps = 0.1) {
    param_weight_detail = seq(from = 0, to = 1, by = 0.1)
    param_grid = expand.grid(param_weight_detail, param_weight_detail, param_weight_detail)
    param_grid = param_grid[apply(param_grid[, 1:3], 1, sum) == 1, ]
    colnames(param_grid) = c("Reproducibility", "Abundance", "Specificity")
    param_grid = data.frame(param_grid, ID = paste(param_grid$Reproducibility, param_grid$Abundance, param_grid$Specificity, 
        sep = "|"), stringsAsFactors = F)
    param_grid
}

mist.train.getPredictionRates = function(prediction, thresholds = seq(from = 0, to = 0.99, by = 0.01)) {
    pos = nrow(prediction[prediction$label == 1, ])
    neg = nrow(prediction[prediction$label == 0, ])
    pos_false = c()
    pos_true = c()
    neg_true = c()
    neg_false = c()
    ID = unique(prediction$ID)
    for (t in thresholds) {
        pos_predicted = prediction[prediction$predicted_scores >= t, ]
        neg_predicted = prediction[prediction$predicted_scores < t, ]
        pf = nrow(pos_predicted[pos_predicted$label == 0, ])
        pos_false = c(pos_false, pf)
        pos_true = c(pos_true, nrow(pos_predicted) - pf)
        nf = nrow(neg_predicted[neg_predicted$label == 1, ])
        neg_false = c(neg_false, nf)
        neg_true = c(neg_true, nrow(neg_predicted) - nf)
    }
    res = data.frame(R_A_S = ID, threshold = thresholds, tpr = pos_true/pos, fpr = pos_false/neg, specificity = neg_true/neg, 
        precision = pos_true/(pos_true + pos_false), fdr = pos_false/(pos_false + pos_true), acc = (pos_true + 
            neg_true)/(pos + neg), f1 = (2 * pos_true)/((2 * pos_true) + pos_false + neg_false), total_pos = pos_true + 
            pos_false)
    res <- arrange(res,tpr,fpr)
    res$AUC <- trapz(res$fpr, res$tpr)
    return(res)
}

# plotROC = function(prediction_rates){ 
#   dim_predictions = round(sqrt(length(unique(prediction_rates$ID))))
#   prediction_rates = prediction_rates[with(prediction_rates,order(ID,tpr,fpr)),] prediction_auc_summ =
# unique(prediction_rates[,c('ID','auc')]) p = ggplot(data=prediction_rates, aes(x=fpr,y=tpr)) p + geom_line()
# + facet_wrap(~ID, nrow=dim_predictions, ncol=dim_predictions) + geom_text(aes(x = 0.6, y = 0.1, label =
# round(auc,digits=3)), parse=T) + theme(axis.text.x=element_text(angle=-90)) }

## generic grid search
mist.train.gridSearch = function(metrics, param_grid, predictionFun = mist.train.getPredictionRates) {
    prediction_rates = c()
    for (i in 1:nrow(param_grid)) {
        param_set = param_grid[i, ]
        # print(sprintf("%s/%s : R %s A %s S %s", i, nrow(param_grid), param_set$Reproducibility, param_set$Abundance, 
        #     param_set$Specificity))
        predicted_scores = (metrics$Reproducibility * param_set$Reproducibility) + (metrics$Abundance * param_set$Abundance) + 
            (metrics$Specificity * param_set$Specificity)
        prediction = data.frame(metrics, param_set, predicted_scores = predicted_scores)
        prediction_rates = rbind(prediction_rates, predictionFun(prediction))
    }
    prediction_rates
}

mist.train.getConfMat <- function(metrics, parameters, confusion_matrix) {
  cat(paste0("\t\tTest Parameters: R ", parameters$R, " A ", parameters$A,
               " S ",parameters$S, "\t Threshold: ", parameters$Thres, "\n")) 
  predicted_scores <- (metrics$Reproducibility * parameters$R) + (metrics$Abundance * parameters$A) + 
    (metrics$Specificity * parameters$S)
  prediction <- data.frame(metrics, parameters, predicted_scores = predicted_scores)
  pos_predicted <- filter(prediction, predicted_scores > parameters$Thres)
  neg_predicted <- filter(prediction, predicted_scores < parameters$Thres)
  
  FP <- nrow(filter(pos_predicted, label == 0))
  TP <- nrow(pos_predicted) - FP
  FN <- nrow(filter(neg_predicted, label == 1))
  TN <- nrow(neg_predicted) - FN
  
  confusion_matrix$FP <- confusion_matrix$FP + FP 
  confusion_matrix$TP <- confusion_matrix$TP + TP 
  confusion_matrix$FN <- confusion_matrix$FN + FN 
  confusion_matrix$TN <- confusion_matrix$TN + TN
  
  return (confusion_matrix)
}

mist.train.kfolds.label <- function(set, k_folds, label) {
  set <- sample_n(set, nrow(set))
  set$group <- cut(seq(1, nrow(set)), breaks = k_folds, labels = F)
  set$label <- rep(label , nrow(set))
  return (set)
} 

mist.train.label = function(metrics, training_set, controls, k_folds) {
  set.seed(10) # for development purposes only
  metrics <- data.table(metrics, key = "Bait")
  metrics <- metrics[!J(controls),]
  setkey(metrics, Bait, Prey)
  write.table(metrics, file = "debug_metrics2.txt", sep = "\t", quote = F)
  positive_set <- mist.train.kfolds.label(filter(metrics[J(training_set), ], !is.na(Abundance)), 
                                          k_folds, 1)
  negative_set <- mist.train.kfolds.label(metrics[!J(training_set), ], 
                                          k_folds, 0)
  cat(paste0("\tTraining Parameters:\n\t\tpos n: ", nrow(positive_set),
             "\n\t\tk: ", k_folds, "\n"))
  
  return (rbind(positive_set, negative_set))
}

mist.train.main = function(metrics, training_set, output_file, training_steps = 0.1, controls, k_folds = 5) {
    # write.table(metrics, file = "debug_metrics1.txt", sep = "\t", quote = F)
    metrics <- mist.train.label(metrics, training_set, controls, k_folds)
    write.table(metrics, file = paste0(output_file, "k_fold_groups.txt"),
                                       sep = "\t", quote = F, eol = "\n", row.names = F)
    param_grid = mist.train.getParamGrid(training_steps)
    confusion_matrix <- c()
    confusion_matrix$TP <- 0
    confusion_matrix$TN <- 0
    confusion_matrix$FP <- 0
    confusion_matrix$FN <- 0
    weights <- c(0,0,0)
    for (i in 1:k_folds) {
      cat(paste0("\tFold: ", i, "/", k_folds, "\n"))
      prediction_rates <- mist.train.gridSearch(filter(metrics, group != i), param_grid)
      write.table(prediction_rates, file = paste0(output_file, i, ".txt"),
                  sep = "\t", quote = F, eol = "\n", row.names = F)
      
      prediction_rates <- arrange(prediction_rates, desc(AUC))
      top_prediction_rates <- arrange(filter(prediction_rates, R_A_S == prediction_rates$R_A_S[1]), desc(f1))
      prediction_weights <- c(as.numeric(unlist(strsplit(as.character(prediction_rates$R_A_S[1]), "\\|"))))
      prediction_params <- data.frame(R = prediction_weights[1], A = prediction_weights[2], S = prediction_weights[3], Thres =  top_prediction_rates$threshold[1],
                                      stringsAsFactors = F)
      confusion_matrix <- mist.train.getConfMat(filter(metrics, group == i), prediction_params, confusion_matrix)
      weights <- weights + (1/k_folds)*prediction_weights
    }
    cat(paste0("\t", k_folds, "-fold cross validation statistics: "))
    cat(sprintf("\n\t\tsensitivity: %s \n\t\tspecificity: %s \n\t\tprecision: %s \n\t\taccuracy: %s \n\t\tf1: %s \n",
                confusion_matrix$TP/(confusion_matrix$TP + confusion_matrix$FN),
                confusion_matrix$TN/(confusion_matrix$TN + confusion_matrix$FP),
                confusion_matrix$TP/(confusion_matrix$TP + confusion_matrix$FP),
                (confusion_matrix$TP + confusion_matrix$TN)/(confusion_matrix$TP + confusion_matrix$FP + confusion_matrix$FN + confusion_matrix$TN),
                (2*confusion_matrix$TP)/(2*confusion_matrix$TP + confusion_matrix$FP + confusion_matrix$FN)))

    weights <- data.frame(R = weights[1], A = weights[2], S = weights[3], stringsAsFactors = F)
    mist_scores = metrics$Reproducibility*weights$R + metrics$Abundance *weights$A + 
      metrics$Specificity *weights$S
    ID <- paste0(weights$R, "|", weights$A, "|", weights$S)
    prediction <- data.frame(metrics, ID, predicted_scores = mist_scores, stringsAsFactors = F)
    prediction_rates = mist.train.getPredictionRates(prediction)
    write.table(prediction_rates, file = paste0(output_file, "prediction_stats.txt"), sep = "\t", row.names = F)
    cat(sprintf("\tWEIGHTS BASED ON TRAINING SET:\n\t  REPRODUCIBILITY: %s\n\t  ABUNDANCE: %s\n\t  SPECIFICITY: %s\n", 
                weights$R, weights$A, weights$S))
    cat("\tAUC: ", prediction_rates$AUC[1], "\tThreshold:", arrange(prediction_rates, desc(f1))$threshold[1], "\n")
    
    return (weights)
}

# metrics = read.delim('tests/entero/processed/preprocessed_MAT_MIST_METRICS.txt',stringsAsFactors=F)
# training_set = read.delim('tests/entero/input/PV_goldset.txt',stringsAsFactors=F, header=F)
# colnames(training_set) = c('Bait','Prey') mist.train.main(metrics, training_set) 