#### get performance and SVM kernel specification ####
get_performance_and_SVMkernel = function(file_location){
  # file_location: file location within the repository for a ".rds" file 
  #                containing information of the AutoML process. Needs a list
  #                "torun.samples" containing a list per learner.
  #                Each list entry must be a dataframe with at least two 
  #                columns "taskname" and "score".
  run_models = readRDS(file_location)
  learner.keys = names(run_models$torun.samples)
  task.keys = unique(as.vector(sapply(run_models$torun.samples, 
                                      function(x) unique(x$taskname))))
  performance = list()
  kernel = list()
  for(task.key in task.keys){
    performance[[task.key]] = list()
    kernel[[task.key]] = c()
    for(learner.key in learner.keys){
      dt = run_models$torun.samples[[learner.key]]
      performance[[task.key]][[learner.key]] = dt[taskname == task.key]$score
      if(learner.key == "svm"){
        kernel[[task.key]] = c(kernel[[task.key]],paste0("svm_",dt[taskname == task.key]$svm.kernel))
      }
    }
  }
  return(list(performance = performance, kernel = kernel))
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