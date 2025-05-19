# DeeDeeExperiment methods -----------------------------------------------------

#' @export
setGeneric("dea_info", function(x, ...) standardGeneric("dea_info"))

#' @export
setGeneric("dea_info<-", function(x, value) standardGeneric("dea_info<-"))

#' @export
setGeneric("dea_names", function(x, value) standardGeneric("dea_names"))

#' @export
setGeneric("dea_rename", function(x, old_name, new_name) standardGeneric("dea_rename"))

#' @export
setGeneric("add_dea", function(x, dea, ...) standardGeneric("add_dea"))

#' @export
setGeneric("remove_dea", function(x, dea_name, ...) standardGeneric("remove_dea"))

#' @export
setGeneric("dea", function(x, dea_name = NULL, ...) standardGeneric("dea"))

#' @export
setGeneric("get_dea_list", function(x, ...) standardGeneric("get_dea_list"))



#' @export
setGeneric("fea_info", function(x, ...) standardGeneric("fea_info"))

#' @export
setGeneric("fea_info<-", function(x, value) standardGeneric("fea_info<-"))

#' @export
setGeneric("fea_names", function(x, value) standardGeneric("fea_names"))

#' @export
setGeneric("fea_rename", function(x, old_name, new_name) standardGeneric("fea_rename"))

#' @export
setGeneric("add_fea", function(x, fea, ...) standardGeneric("add_fea"))

#' @export
setGeneric("remove_fea", function(x, fea_name, ...) standardGeneric("remove_fea"))

#' @export
setGeneric("fea", function(x, fea_name = NULL, ...) standardGeneric("fea"))

#' @export
setGeneric("get_fea_list", function(x, dea_name = NULL, ...) standardGeneric("get_fea_list"))


