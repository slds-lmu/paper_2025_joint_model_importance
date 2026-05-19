library(ggplot2)
library(data.table)
library(purrr)
library(batchtools)
library(mlr3)
library(mlr3learners)
library(mlr3data)
library(mlr3tuning)
library(mlr3pipelines)
library(mlr3fairness)
library(paradox)
library(checkmate)
library(xgboost)

for (loading in c("assets", "batchtools")) {
  list.files(file.path("init", loading), "\\.r$", ignore.case = TRUE, full.names = TRUE) |>
    lapply(source, verbose = FALSE) |>
    invisible()
}