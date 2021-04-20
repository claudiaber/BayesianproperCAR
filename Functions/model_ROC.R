modelROC <- function(flt_predictions, flt_real_labels, str_title){
  
  if(base::missing(str_title)){str_title <- "Binary Performance"}
  
  # Check input data
  # checkmate::assert(
  #   checkmate::checkString(str_title,
  #                          min.chars = 1),
  #   checkmate::checkNumeric(flt_predictions, 
  #                           lower = 0, upper = 1, 
  #                           any.missing = FALSE),
  #   checkmate::checkNumeric(flt_real_labels, 
  #                           lower = 0, upper = 1, 
  #                           any.missing = FALSE),
  #   checkmate::check_set_equal(x = flt_real_labels, 
  #                              y = c(0, 1)), 
  #   combine = "and"
  # )
  
  x <- NULL
  y <- NULL
  
  flt_real_labels <- base::as.integer(flt_real_labels)
  
  obj_prediction <- ROCR::prediction(predictions = flt_predictions, 
                                     labels = flt_real_labels, 
                                     label.ordering = c(0,1))
  
  # ROC, AUC
  obj_auc <- ROCR::performance(prediction.obj = obj_prediction, measure = "auc")
  obj_roc <- ROCR::performance(prediction.obj = obj_prediction, measure = "tpr", x.measure = "fpr")
  flt_auc <- base::as.numeric(base::format(obj_auc@y.values[[1]], digits = 3))
  dtf_tmp <- base::data.frame(x = obj_roc@x.values[[1]], y = obj_roc@y.values[[1]])
  
  plot_roc <- ggplot2::ggplot(data = dtf_tmp, ggplot2::aes(x = x, y = y, group = 1)) + 
    ggplot2::geom_line(color = "blue", size = 1.5) + 
    ggplot2::scale_x_continuous(expand = c(0, 0), limits = c(0, 1), minor_breaks = base::seq(0 , 1, 0.05)) +
    ggplot2::scale_y_continuous(expand = c(0, 0), limits = c(0, 1), minor_breaks = base::seq(0 , 1, 0.05)) + 
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
    ggplot2::labs(title = str_title, subtitle = "ROC") +
    ggplot2::xlab("False positive rate") + 
    ggplot2::ylab("True positive rate") + 
    ggplot2::geom_text(x = 0.7, y = 0.3, label = base::paste0("AUC = ", flt_auc), size = 7)
  
  lst_out <- list(
    AUC = flt_auc, 
    PLOT_ROC = plot_roc
  )
  
  return(lst_out)
}
