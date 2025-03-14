# DeeDeeExperiment methods -----------------------------------------------------

#' @export
setGeneric("dea", function(x, ...) standardGeneric("dea"))

#' @export
setGeneric("dea<-", function(x, value) standardGeneric("dea<-"))

#' @export
setGeneric("add_dea", function(x, dea, ...) standardGeneric("add_dea"))

#' @export
setGeneric("remove_dea", function(x, dea_name, ...) standardGeneric("remove_dea"))

#' @export
setGeneric("get_dea_df", function(x, dea_name, ...) standardGeneric("get_dea_df"))

#' @export
setGeneric("get_dea_list", function(x, ...) standardGeneric("get_dea_list"))



#' @export
setGeneric("fea", function(x, ...) standardGeneric("fea"))

#' @export
setGeneric("fea<-", function(x, value) standardGeneric("fea<-"))

#' @export
setGeneric("add_fea", function(x, dea, ...) standardGeneric("add_fea"))

#' @export
setGeneric("remove_fea", function(x, fea_name, ...) standardGeneric("remove_fea"))

