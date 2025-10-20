me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

env <- tools::file_path_as_absolute(file.path(lib.dir, "R413"))

args <- commandArgs(TRUE)

mnc_obj_file <- 'Trajectory.obj.rds'
outdir <- './'

if (anyNA(c(mnc_obj_file, outdir))) {
  warning("\n  Usage: ", env, "/bin/Rscript ", me, "<Trajectory.obj.rds> <outdir>\n", call. = FALSE, immediate. = TRUE)
  quit()
}

library(monocle)
library(igraph)
library(Seurat)
library(dplyr)

source(tools::file_path_as_absolute(file.path(lib.dir, "Monocle_lib2.R")), chdir = TRUE)

mnc_obj_file <- tools::file_path_as_absolute(mnc_obj_file)

dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

message("==> Read monocle object: ", mnc_obj_file)
mnc_obj <- readRDS(mnc_obj_file)

if (!"Groups" %in% colnames(pData(mnc_obj))) {
  pData(mnc_obj)$Groups <- pData(mnc_obj)$Samples
  experimentData(mnc_obj)@other$color.group <- experimentData(mnc_obj)@other$color.sample
}

message("==> stat Trajectory <==")
StatTrajectory(mnc_obj)
df_list <- GetTrajectoryData(mnc_obj, is.return = TRUE) # for online
df_list$cells <- df_list$cells %>% mutate(Groups = pData(mnc_obj)$Groups, .after = Clusters)
WriteTable(df_list$cells, "Trajectory.cell.data.xls")
WriteTable(df_list$bone, "Trajectory.bone.data.xls")

assay <- experimentData(mnc_obj)@other$assay %||% "RNA"
StatCluster(pData(mnc_obj), group.by = "Samples", stat.what = "State", outpref = "State.stat", assay = assay, color.gb = experimentData(mnc_obj)@other$color.sample)
StatCluster(pData(mnc_obj), group.by = "Clusters", stat.what = "State", outpref = "State.stat", assay = assay, color.gb = experimentData(mnc_obj)@other$color.cluster)
if ("Groups" %in% colnames(pData(mnc_obj))) {
  StatCluster(pData(mnc_obj), group.by = "Groups", stat.what = "State", outpref = "State.stat", assay = assay, color.gb = experimentData(mnc_obj)@other$color.group)
}
CDS_avg(mnc_obj)

message("==> Done <==")

