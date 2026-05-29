#### fix bs task: FeatureImp does not handle logical features, convert to factor ####
fix_bs_task_for_featureimp <- function(task) {
  task_data <- as.data.frame(task$data())
  task_id <- task$id
  task_target <- task$target_names
  for (i in seq_along(task_data)) {
    if (is.logical(task_data[[i]])) task_data[[i]] <- as.factor(task_data[[i]])
  }
  as_task_regr(task_data, target = task_target, id = task_id)
}

#### wrap bs model with a factor->logical pipeline ####
# bs models were trained on logical holiday / working_day columns. Once we
# convert those to factor (see fix_bs_task_for_featureimp), the saved model
# cannot predict directly — the convert_types pipeop converts them back so
# the underlying learner sees what it was trained on.
fix_bs_model_for_predict <- function(model, instance) {
  holiday_special <- ppl(
    "convert_types", "factor", "logical",
    selector_name(c("holiday", "working_day")),
    id = "holiday.special"
  )
  invisible(holiday_special$train(instance))
  xstate <- model$state
  gr <- holiday_special$clone(deep = TRUE) %>>% model$clone(deep = TRUE)
  lr <- as_learner(gr$clone(deep = TRUE))
  lr$state <- xstate
  lr$state$train_task <- instance$clone(deep = TRUE)$filter(0)
  lr$model <- gr$state
  lr$model[[gr$ids()[[2]]]] <- xstate
  lr
}

#### build nested performance list from the flat res_dt format ####
# res_dt (from data/results_modelperformances.RData) has columns
#   task, learner, model.no, test.score, score (= score type label)
# get_RS() expects a nested list performance[[task]][[learner]] = numeric
# vector of test scores ordered by model.no, parallel to the model columns
# in vic[[task]].
#
# Important: every task must carry an entry for every learner (an empty
# numeric vector is fine) so that sapply(performance, sapply, min) is a
# rectangular matrix and not a ragged list. Without this, apply(., 2, min)
# in get_RS() fails with "dim(X) must have a positive length".
build_performance_from_res_dt = function(res_dt, tasks = NULL, learners = NULL) {
  if (is.null(tasks))    tasks    = sort(unique(res_dt$task))
  if (is.null(learners)) learners = sort(unique(res_dt$learner))
  performance = list()
  for (t in tasks) {
    performance[[t]] = list()
    sub_t = res_dt[task == t]
    for (l in learners) {
      sub_tl = sub_t[learner == l]
      if (nrow(sub_tl) == 0) {
        performance[[t]][[l]] = numeric(0)
      } else {
        data.table::setorder(sub_tl, model.no)
        performance[[t]][[l]] = sub_tl$test.score
      }
    }
  }
  performance
}

#### get Rashomon set ####
# Returns a list `model_index_RS` with one entry per task. Each entry is an
# integer vector of column positions inside vic[[task]]'s model columns
# (1-indexed, post the leading feature column). Equivalently these are the
# row positions inside design[rn == task], so callers can use them to
# subset either object.
#
# A model is in the Rashomon set when its test score is strictly below
#   best_performance(task) * (1 + epsilon).
#
# We walk the VIC columns directly and look the matching score up in
# `performance` by (learner, model.no) parsed from the column name. This
# is robust against the cases that previously broke the cumsum-offset
# version:
#   - learner present in `performance` but not in vic[[task]] -> skipped
#   - tasks where vic has fewer learners than the global learner.keys
#   - performance vector longer than the corresponding vic block (surplus
#     model.no's simply aren't iterated over)
get_RS = function(epsilon, performance, vic) {
  # min(numeric(0)) is Inf with a warning -- harmless when a learner had
  # no models for some task, but the warnings are noisy.
  best_performance = suppressWarnings(apply(sapply(performance, sapply, min), 2, min))
  task.keys = names(vic)

  model_index_RS = list()
  for (task.key in task.keys) {
    threshold  = best_performance[task.key] * (1 + epsilon)
    vic_cnames = colnames(vic[[task.key]])[-1]   # drop feature column
    if (length(vic_cnames) == 0) {
      model_index_RS[[task.key]] = integer(0)
      next
    }
    vic_learner  = sub(".*_(.*?)_.*", "\\1", vic_cnames)
    vic_model_no = as.integer(sub(".*_m([0-9]+)$", "\\1", vic_cnames))

    rs = integer(0)
    for (i in seq_along(vic_cnames)) {
      l   = vic_learner[i]
      mno = vic_model_no[i]
      perf_vec = performance[[task.key]][[l]]
      if (is.null(perf_vec) || length(perf_vec) == 0) next
      if (is.na(mno) || mno > length(perf_vec)) next
      score = perf_vec[mno]
      if (!is.na(score) && score < threshold) rs = c(rs, i)
    }
    model_index_RS[[task.key]] = rs
  }
  model_index_RS
}
