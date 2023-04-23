#### Array ops ####
#
# Only execute once to prepare the cpp code.
#
# library(einsum)
# func_list=list(
#   array_mult_left=einsum_generator('ixk,xj -> ijk',FALSE),
#   array_mult_right=einsum_generator('xik,jx -> jik',FALSE),
#   array_transp=einsum_generator('ijk -> jik',FALSE),
#   array_collapse_left=einsum_generator('ixk,x -> ik',FALSE),
#   array_collapse_right=einsum_generator('xik,x -> ik',FALSE)
# )
#
# for(f_name in names(func_list)){
#   write(paste('#include <Rcpp.h>',
#               'using namespace Rcpp;',
#               '',
#               '// [[Rcpp::export]]',
#               stringr::str_replace_all(func_list[[f_name]],
#                                        'einsum_impl_func',
#                                        f_name),
#               sep='\n'),
#         paste0('src/',f_name,'.cpp'))
# }


# #' @useDynLib kDGLM
# #' @importFrom Rcpp sourceCpp
# NULL


#' array_mult_left
#'
#' @param A A 3-D array with shapes nxmxk.
#' @param B A matrix with shapes mxl
#'
#' @export
array_mult_left <- function(A, B) {
  apply(A, 3, function(x) {
    x %*% B
  }, simplify = FALSE) %>% simplify2array(except = NULL)
}

#' array_mult_right
#'
#' @param A A 3-D array with shapes nxmxk.
#' @param B A matrix with shapes lxn
#'
#' @export
array_mult_right <- function(A, B) {
  apply(A, 3, function(x) {
    B %*% x
  }, simplify = FALSE) %>% simplify2array(except = NULL)
}

#' array_transp
#'
#' @param A A 3-D array
#'
#' @export
array_transp <- function(A) {
  apply(A, 3, function(x) {
    t(x)
  }, simplify = FALSE) %>% simplify2array(except = NULL)
}

#' array_collapse_left
#'
#' @param A A 3-D array with shapes nxmxk.
#' @param B A matrix with shapes mx1
#'
#' @export
array_collapse_left <- function(A, B) {
  array_mult_left(A, B)[, 1, ]
}

#' array_collapse_right
#'
#' @param A A 3-D array with shapes nxmxk.
#' @param B A matrix with shapes 1xn
#'
#' @export
array_collapse_right <- function(A, B) {
  array_mult_right(A, B)[1, , ]
}
