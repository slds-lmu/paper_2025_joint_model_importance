load("data/design.RData")
load("data/results_vic.RData")

task.keys = names(vic) # german credit, compas, bike sharing, synthetic
learner.keys = c("tree", "glmnet", "xgb", "nnet", "svm")

list.tasks <- list(
  gc = task.gc,
  cs = task.cs,
  bs = task.bs,
  st = task.st
)


list.learners.regr <- list(
  xgb = learner.xgb.regr,
  tree = learner.tree.regr,
  nnet = learner.nnet.regr,
  glmnet = learner.regr.glmnet,
  svm = learner.svm.regr
)

list.learners.classif <- list(
  xgb = learner.xgb.classif,
  tree = learner.tree.classif,
  nnet = learner.nnet.classif,
  glmnet = learner.classif.glmnet,
  svm = learner.svm.classif
)