load("data/design_all_but_TreeFARMS.RData")     # `design`, `pre_design`
load("data/results_vic_all_but_TreeFARMS.RData")  # `vic`, `vic_normalized`
load("data/results_modelperformances.RData")     # `res_dt`: task, learner, model.no, test.score, score

task.keys = names(vic)

# Auto-detect learner.keys from the VIC column names so the set stays in
# sync with the current model generation. (The previous hardcoded list
# c("tree", "glmnet", "xgb", "nnet", "svm") missed e.g. "svm.linear" /
# "svm.radial" / "gosdt" in newer runs, which silently dropped models
# from the Rashomon set.) Only learners that have VIC columns are
# considered downstream -- a learner present in res_dt but not in vic
# would break get_RS()'s column-index bookkeeping.
learner.keys = sort(unique(unlist(lapply(vic, function(v) {
  cnames <- colnames(v)
  if (length(cnames) <= 1) character(0)
  else unique(sub(".*_(.*?)_.*", "\\1", cnames[-1]))
}))))

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