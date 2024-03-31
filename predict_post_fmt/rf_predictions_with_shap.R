#!/usr/bin/Rscript

##-------------------------------##
## Predict Single Taxa with SHAP ##
##-------------------------------##

# This code is designed to predict a single taxa using SHAP (SHapley Additive exPlanations) values.

# Before running the code, ensure that the following variables are properly set:
# - `output_path`: Path to the directory where output files will be saved.
# - `feature_table_for_ml_path`: Path to the feature table.
# - `single_predicted_taxa`: The name of the taxa to be predicted.

## Output Files

# The function generates the following output files:
  
# - `{output_path}/{predicted_taxa}_pairs_metrics.txt`: Metrics including RMSE, R-squared, and Spearman correlation.
# - `{output_path}/{predicted_taxa}_pairs_predictions.txt`: Model predictions.
# - `{output_path}/{predicted_taxa}_pairs_shap.txt`: SHAP values (if calculated).

## Usage
# predict_single_taxa_with_shap(data = feature_table_for_ml, predicted_taxa = single_predicted_taxa)

### Set up paths and libraries
output_path <- "output_path"
feature_table_for_ml_path <- "feature_table_for_ml_path"
single_predicted_taxa <-args[[1]]

### 

library(fastshap)
library(tidymodels)
library(readr)
library(stringr)

print("load library")

args = commandArgs(trailingOnly=TRUE)

print("load args")

predict_single_taxa_with_shap <- function(data, predicted_taxa){
  
  # orgnize data
  predicted_values <- str_c("recipient_post_fmt_", predicted_taxa)
  
  taxa_names <- names(data)
  taxa_names <- taxa_names[!(taxa_names %in% c("subject_post_timepoint_id", "post_timepoint", "subject_id"))]
  taxa_names <- str_remove_all(taxa_names, pattern = c( "donor_pre_fmt_"))%>%
    str_remove_all(., pattern = c("recipient_pre_fmt_"))%>%
    str_remove_all(., pattern = "recipient_post_fmt_")%>%
    unique()
  
  features_used <- map(.x = taxa_names, .f = ~c(str_c("donor_pre_fmt_", .x), str_c("recipient_pre_fmt_", .x)))%>%
    unlist()
  
  columns_to_keep <- c(features_used, predicted_values)
  formula <- as.formula(str_c(predicted_values, "~."))
  
  filtred_data <- select(data , all_of(columns_to_keep), subject_post_timepoint_id, subject_id,post_timepoint ) 
  data_ids <- select(filtred_data, subject_post_timepoint_id, subject_id,post_timepoint)
  
  data_cv <- group_vfold_cv(filtred_data, group = subject_id , v= 10)
  
  print("Data ready...")
  
  #create model
  recipe <- recipe(formula, data = filtred_data)%>%
    update_role(subject_post_timepoint_id,
                subject_id,
                post_timepoint,
                new_role = "ID")
  
  cores <- parallel::detectCores()
  rf_model <- 
    rand_forest() %>%
    #set_args(mtry = tune(), trees = tune()) %>%
    set_engine("ranger" ,  num.threads = cores, importance = "impurity") %>%
    set_mode("regression")
  
  model_workflow <- workflow() %>%
    add_recipe(recipe) %>%
    add_model(rf_model)
  
  #Test the model using 10 fold-cv
  model_results<- fit_resamples(model_workflow, data_cv, 
                                control = control_resamples(save_pred = TRUE),
                                metrics = metric_set(rmse, rsq))
  
  #collect predictions 
  model_predictions <-collect_predictions(model_results)%>%
    rename(predicted_values = .pred,
           true_values = !!sym(predicted_values),
           row = .row)%>%
    arrange(row)%>%
    bind_cols(data_ids)%>%
    select(subject_post_timepoint_id, subject_id,post_timepoint, predicted_values, true_values)%>%
    mutate(predicted_taxa = predicted_taxa)
  
  #collect metric
  spearman_correlation <-cor.test(x= model_predictions[["predicted_values"]], y= model_predictions[["true_values"]], method = c( "spearman"))
  spearman_row <- tibble(metric = "spearman",
                         p_value = spearman_correlation$p.value,
                         mean = spearman_correlation$estimate)
  
  model_metrics <- collect_metrics(model_results)%>%
    rename(metric = .metric)%>%
    select(metric, mean, std_err)%>%
    bind_rows(spearman_row)%>%
    mutate(predicted_taxa = predicted_taxa)
  
  #SHAP values
  
  if(spearman_correlation$estimate>=0.2){
  
    print("spearman_correlation>=0.2, Calculate SHAP...")
  
    model_fit <- fit(model_workflow, data = filtred_data)
    X <-filtred_data%>%
      select(-c(all_of(predicted_values)))%>%
      as.data.frame()
    
    pfun <- function(object, newdata) {predict(object, data = newdata)$predictions}
    
    predict_function <-  function(model, newdata) {predict(model, newdata) %>% pluck(.,1)}
    
    shap_values <- fastshap::explain(model_fit, X = X, pred_wrapper = predict_function)
    
    shap_values_tidy <- shap_values%>%
      as_tibble()%>%
      select(-c(subject_post_timepoint_id, subject_id,post_timepoint))%>%
      bind_cols(data_ids)
      # pivot_longer(-c(subject_post_timepoint_id, subject_id,post_timepoint), names_to = "feature", values_to = "shap")%>%
      # mutate(predicted_taxa = predicted_taxa)%>%
      # filter(shap>0)
   
    print("Save SHAP files...")
    path_shap <-  str_c("output_path", predicted_taxa, "_pairs_shap.txt")
    write_csv(shap_values_tidy, path_shap)
    
    }else{
      print("spearman_correlation<0.2, Skip SHAP calculation")
    }
  
  print("Save metrics and predictions...")
  
  path_metrics <- str_c("output_path", predicted_taxa, "_pairs_metrics.txt")
  path_predictions <- str_c("output_path", predicted_taxa, "_pairs_predictions.txt")
  
  write_csv(model_metrics, path_metrics)
  write_csv(model_predictions, path_predictions)
  
  # return(list(model_predictions = model_predictions,
  #             model_metrics = model_metrics,
  #             shap_values = shap_values_tidy))
  
}

print("load function")

feature_table_for_ml <- readRDS(feature_table_for_ml_path)

print(str_c("Working on: ", single_predicted_taxa))

print("load data")

predict_single_taxa_with_shap(data =feature_table_for_ml_path, predicted_taxa = single_predicted_taxa)

print("Finished!")