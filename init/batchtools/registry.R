getRegistry <- function(regpath, make.default = FALSE, writeable = TRUE) {
  if (!file.exists(regpath)) {
    makeExperimentRegistry(
      file.dir = regpath,
      source = "init.R",
      seed = 1,
      make.default = make.default
    )
  } else {
    loadRegistry(file.dir = regpath, writeable = writeable, make.default = make.default)
  }
}