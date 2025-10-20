
me <- normalizePath(sub("--file=", "", grep("--file=", commandArgs(), value = TRUE)))
lib.dir <- normalizePath(dirname(me))

args <- commandArgs(TRUE)

file <- 'config.yaml'
mnc_obj_file <- 'Trajectory.obj.rds'
outdir <- './'

if (anyNA(c(file, mnc_obj_file, outdir))) {
  warning("\n  Usage: Rscript ", me, " <parameter.yaml> <Trajectory.obj.rds> <outdir>\n", call. = FALSE, immediate. = TRUE)
  quit()
}

source(tools::file_path_as_absolute(file.path(lib.dir, "Eula.monocle2.R")), chdir = TRUE)

file <- tools::file_path_as_absolute(file)
mnc_obj_file <- tools::file_path_as_absolute(mnc_obj_file)
cores <- 1

handlers <- list(
  "bool#no" = function(x) {if (x %in% c("false", "FALSE")) FALSE else x},
  "bool#yes" = function(x) {if (x %in% c("true", "TRUE")) TRUE else x}
)
parameter <- yaml::yaml.load_file(file, handlers = handlers)
if ("monocle2" %in% names(parameter)) {
  parameter <- parameter$monocle2
}

dir.create(outdir, FALSE, TRUE)
outdir <- tools::file_path_as_absolute(outdir)
setwd(outdir)

message("==> Read monocle object: ", mnc_obj_file)
mnc_obj <- readRDS(mnc_obj_file)

### Diff
message( ">>> Differ <<<" )
if ( "diff_state" %in% parameter$analyse && nlevels(pData(mnc_obj)$State) > 1 ) {
  message( "==> state diff <== " )
  setwd(outdir)
  dir.create(path = "Diff_State", recursive = T, showWarnings = F)
  setwd("Diff_State")

  parameter$diff_state$thres <- parameter$diff_state$thres %||% 1e-7
  parameter$diff_state$use.q <- parameter$diff_state$use.q %||% TRUE
  
  diff_state_res <- Diff.State(mnc_obj, cores = cores, thres = parameter$diff_state$thres, use.q = parameter$diff_state$use.q)
  GetGenesInPseudotime(mnc_obj[rownames(diff_state_res)[1], ], cores = cores)
  setwd(outdir)
}

if ( "diff_pseudotime" %in% parameter$analyse ) {
  message( "==> pseudotime diff <== " )
  setwd(outdir)
  dir.create(path = "Diff_Pseudotime", recursive = T, showWarnings = F)
  setwd("Diff_Pseudotime")

  parameter$diff_pseudotime$thres <- parameter$diff_pseudotime$thres %||% 1e-7
  parameter$diff_pseudotime$use.q <- parameter$diff_pseudotime$use.q %||% TRUE

  diff_Pseudotime_sig <- Diff.Pseudotime(
    mnc_obj, 
    cores = cores, 
    thres = parameter$diff_pseudotime$thres, 
    use.q = parameter$diff_pseudotime$use.q
  )
  GetGenesInPseudotime(mnc_obj[rownames(diff_Pseudotime_sig)[1], ], cores = cores)
  setwd(outdir)
}

if ( "diff_branch" %in% parameter$analyse && nlevels(pData(mnc_obj)$State) > 1) {
  message( "==> branch diff <== " )
  setwd(outdir)
  dir.create(path = "Diff_Branch", recursive = T, showWarnings = F)
  setwd("Diff_Branch")

  parameter$diff_branch$thres <- parameter$diff_branch$thres %||% 1e-7
  parameter$diff_branch$use.q <- parameter$diff_branch$use.q %||% TRUE

  total_BEAM_sig <- Diff.Branch(
    mnc_obj, 
    cores = cores, 
    thres = parameter$diff_branch$thres, 
    use.q = parameter$diff_branch$use.q
  )
  if (length(total_BEAM_sig) > 0) {
    branch_name <- strsplit(as.character(total_BEAM_sig$branch[1]), " -vs- ")[[1]]
    GetGenesBranchedPseudotime(mnc_obj[as.character(total_BEAM_sig$GeneID[1]), ], branch_point = total_BEAM_sig$branch_node[1], branch_labels = branch_name, cores = cores)
  }
  setwd(outdir)
}

message( ">>> __ALL DONE__ <<<" )

