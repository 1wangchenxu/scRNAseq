me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "R413"))

args <- commandArgs(TRUE)

file <- 'config.yaml'
obj_file <- 'obj.Rda'
outdir <- './'

if (anyNA(c(file, obj_file, outdir))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, " <parameter.yaml> <obj_file> <outdir>\n", call. = FALSE, immediate. = TRUE)
  quit()
}

library(monocle)
library(igraph)
library(dplyr)
library(Seurat)
library(dplyr)

source(tools::file_path_as_absolute(file.path(lib.dir, "Monocle_lib.R")), chdir = TRUE)

file <- tools::file_path_as_absolute(file)

handlers <- list(
  "bool#no" = function(x) {if (x %in% c("false", "FALSE")) FALSE else x}, 
  "bool#yes" = function(x) {if (x %in% c("true", "TRUE")) TRUE else x}
)
parameter <- yaml::yaml.load_file(file, handlers = handlers)
if ("monocle2" %in% names(parameter)) {
  parameter <- parameter$monocle2	
}

parameter$analyse <- parameter$analyse %||% c("run_monocle", "diff_state", "diff_pseudotime", "diff_branch")
parameter$root_state <- parameter$root_state %||% "1"
print(parameter)

obj_file <- tools::file_path_as_absolute(obj_file)


message(">>>>> Initial monocle object")

mnc_obj <- MonocleObject(parameter = parameter, obj_file = obj_file)
gc(verbose = FALSE)

if (!is.null(parameter$feature_use)) {
  features <- readLines(parameter$feature_use)
  features <- intersect(rownames(mnc_obj), features)
  mnc_obj <- setOrderingFilter(mnc_obj, features)
}

dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

if (is.null(parameter$mnc_obj) || !file.exists(parameter$mnc_obj)) {
  message("==> Run monocle <==")
  mnc_obj <- RunMonocle(mnc_obj, parameter)
}

message("==> Rename States <==") 
mnc_obj <- RenameState(mnc_obj, parameter)

message("==> Add misc data <==")
obj <- Load(obj_file)
misc.names <- names(obj@misc)
experimentData(mnc_obj)@other <- obj@misc[setdiff(misc.names, c("counts", "pdata", "fdata"))]
experimentData(mnc_obj)@other$assay <- DefaultAssay(obj)

message("==> Save mnc_obj: Trajectory.obj.rds")
if (!"Groups" %in% colnames(pData(mnc_obj))) {
  pData(mnc_obj)$Groups <- pData(mnc_obj)$Samples
  experimentData(mnc_obj)@other$color.group <- experimentData(mnc_obj)@other$color.sample
}

StatTrajectory(mnc_obj)

saveRDS(mnc_obj, file = "Trajectory.obj.rds")

