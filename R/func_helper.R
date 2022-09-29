#' if.null
#'
#' This function is wrapper for ifelse(is.null(.),.,.)
if.null=function(vec,val){ifelse(is.null(vec),val,vec)}
