tableSummaryTex <- function(summarytab,round=2,rownames, label,caption, filename)
{
  DT <- data.frame(lapply(as.data.frame(summarytab),FUN = function(x){
    if(is.numeric(x)){round(x, round)}
    else{x}}),
             row.names=rownames)
  Tabletex <- knitr::kable(DT,
                                     col.names = stringr::str_replace_all(
                                       stringr::str_replace_all(colnames(summarytab),
                                                       pattern = stringr::fixed("_"),
                                                       replacement = "\\_"),
                                       pattern = stringr::fixed("%"),
                                       replacement= "\\%"),
                                     format = "latex", 
                                     caption = paste0("\\label{",label,"}
                                                      ", caption),
                                     escape=FALSE)
  fileTable <-  file(filename)
  writeLines(text = Tabletex, fileTable)
  close(fileTable)
  return(NULL)
}
  