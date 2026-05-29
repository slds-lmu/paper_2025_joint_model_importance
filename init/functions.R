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
get_RS = function(epsilon, performance, vic){
  # epsilon: Prozentsatz, den die Performance der Modelle im RS von der 
  #          Performance des "besten" Modells abweichen darf. (z.B. 0.05 für 5%)
  # performance: Vektor mit den Performance Maßen der Modelle.
  # vic: Liste mit Eintragungen für jeden Datensatz, die jeweils einen data 
  #      frame enthalten, der in der ersten Spalte die features enthält und in 
  #      den übrigen Spalten pro Modell die passenden feature importance Werte.
  #      Die Spalten müssen das Format pfi_[modell]_m[modellnummer] haben, z.B.:
  #      pfi_tree_m1
  best_performance = apply(sapply(performance, sapply, min),2,min)
  task.keys = names(vic)
  learner.keys = unique(as.vector(sapply(performance, 
                                         function(x) names(x))))
  perf_index_RS = list()
  vic_RS = list()
  vic_normalized_RS = list()
  vic_learner_table = list()
  for(task.key in task.keys){
    vic_names = colnames(vic[[task.key]])
    perf_index_RS[[task.key]] = lapply(performance[[task.key]], function(x) which(x < best_performance[task.key]*(1+epsilon)))
    vic_learner = sub(".*_(.*?)_.*", "\\1", vic_names[-1])
    vic_learner_table[[task.key]] = list()
    vic_learner_table[[task.key]]$all = table(vic_learner)[unique(vic_learner)] # also in pre_design so this is a bit useless
    cum_sum_learner = cumsum(vic_learner_table[[task.key]]$all)
    index = 1 # feature column in vic
    for(learner.key in learner.keys){
      if(!is_empty(perf_index_RS[[task.key]][[learner.key]])){
        index_adopt = cum_sum_learner[which(learner.key == learner.keys)-1]
        if(is_empty(index_adopt)){
          model_index = perf_index_RS[[task.key]][[learner.key]]
        } else {
          model_index = perf_index_RS[[task.key]][[learner.key]]+index_adopt
        }
        index = c(index, model_index+1)
      }
    }
    vic_RS[[task.key]] = vic[[task.key]][,index]
    vic_normalized_RS[[task.key]] = vic_normalized[[task.key]][,index]
    perf_index_RS[[task.key]] = index
    
    vic_RS_names = colnames(vic_RS[[task.key]])
    vic_RS_learner = sub(".*_(.*?)_.*", "\\1", vic_RS_names[-1])
    vic_learner_table[[task.key]]$RS = table(vic_RS_learner)[unique(vic_RS_learner)]
  }
  rm(vic_names, vic_learner, cum_sum_learner, index)
  model_index_RS = lapply(perf_index_RS, function(x) x[-1]-1)
}