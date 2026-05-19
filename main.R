source("init.R")

#### settings #################################################################

# run_models_merged location
run_models_merged_location = "data/run_models_merged.rds"

# Rashomon set epsilon
RS_epsilon = 0.01
distance_metrics = c("euclidean", "manhattan")

## Batchtools
regr = getRegistry(regpath = "/media/external/ewaldf/paper_2025_joint_model_importance", make.default = TRUE)
# writeable = TRUE (default) only in one window !!!!

# Define Cluster-Configurations
regr$cluster.functions = makeClusterFunctionsSocket(ncpus = 30)

#### prerequisites ############################################################

# performance and SVM kernel specification
perfor_and_kernel = get_performance_and_SVMkernel(run_models_merged_location)
performance = perfor_and_kernel$performance
kernel = perfor_and_kernel$kernel
rm(perfor_and_kernel)

# get model indices of Rashomon set per task
RS = get_RS(RS_epsilon, performance, vic)
RS_size = data.table(length = unlist(lapply(RS, FUN = length)))
RS_size$task = names(RS)
RS_size$cum_len = cumsum(RS_size$length)


#### calculate distance metrics ###############################################
# Define Problem
addProblem("fromlist", fun = function(data, job, taskname) {
  task = list.tasks[[taskname]]
  # model = readRDS("/media/external/rashomon/datafiles/st/tree/samplemodel_tree_st_0001.rds")
  
  # Fix logical features (for FeatureImp)
  if(taskname == "bs"){
    task_data = as.data.frame(task$data())
    task_id = task$id
    task_target = task$target_names
    for(i in seq_along(task_data)) {
      if (is.logical(task_data[[i]])) task_data[[i]] = as.factor(task_data[[i]])
    }
    task = as_task_regr(task_data, target = task_target, id = task_id)
  }
  
  # Return of the validation split
  generateCanonicalDataSplits(task, ratio = 2 / 3, seed = 1)$validation
})

# Define Algorithm for VIC calculation based on PFI
addAlgorithm("calculate_preds", fun = function(data, instance, job, learnername, model.no) {
  # browser()
  name = sprintf("/media/external/rashomon/datafiles/%s/%s/samplemodel_%s_%s_%04d.rds",
                 job$pars$prob.pars$taskname, learnername, learnername, job$pars$prob.pars$taskname, model.no)
  model = readRDS(name)
  
  # Fix models in case of task bs (logical features)
  if(job$pars$prob.pars$taskname == "bs"){
    # fix model
    holiday.special = ppl("convert_types", "factor", "logical", selector_name(c("holiday", "working_day")), id = "holiday.special")
    invisible(holiday.special$train(instance)) # list.tasks$bs))
    xstate = model$state
    gr = holiday.special$clone(deep = TRUE) %>>% model$clone(deep = TRUE)
    lr = as_learner(gr$clone(deep = TRUE))
    lr$state = xstate
    lr$state$train_task = instance$clone(deep = TRUE)$filter(0)
    lr$model = gr$state
    lr$model[[gr$ids()[[2]]]] = xstate
    model = lr
    rm(gr, lr, xstate, holiday.special)
  }
  
  # Calculate predictions
  model$predict(instance)$response
})

pred_design = design[0]
for(i in names(RS)){#
  subset = design[rn == i]
  pred_design = rbind(pred_design, subset[RS[[i]],])
}

addExperiments(
  prob.designs = list(fromlist = data.table(taskname = pred_design$rn)),
  algo.designs = list(calculate_preds = pred_design[,-"rn"]),
  repls = 1,
  combine = "bind"
)

test = testJob(31)
length(test)

submitJobs(findErrors())
submitJobs()
waitForJobs()

# extract results
res = ijoin(
  getJobPars(),
  reduceResultsDataTable(fun = function(x) list(res = x))
)
res = unwrap(res, sep = ".")
rres = res$result.res
preds = list()
for(i in 1:length(RS_size$length)){
  if(i == 1){
    preds[[RS_size$task[i]]] = data.frame(t(sapply(rres[1:RS_size$cum_len[i]], 
                                                   function(x) x[1:max(lengths(rres[1:RS_size$cum_len[i]]))])))
  } else {
    preds[[RS_size$task[i]]] = data.frame(t(sapply(rres[(RS_size$cum_len[i-1]+1):RS_size$cum_len[i]], 
                                                   function(x) x[1:max(lengths(rres[(RS_size$cum_len[i-1]+1):RS_size$cum_len[i]]))])))
  }
}

test = preds$bs
rownames(test) = paste0("m", RS$bs)
dist_test = dist(test, method = distance_metrics[1])
plot(hclust(dist_test), hang = -1, main = "Cluster Dendrogram bs") # to see a dendrogram of clustered variables

# hc <- hclust(dist_test)
# memb <- cutree(hc, k = 10)
# cent <- NULL
# for(k in 1:10){
#   cent <- rbind(cent, colMeans(test[memb == k, , drop = FALSE]))
# }
# hc1 <- hclust(dist(cent, method = distance_metrics[1]), members = table(memb))
# par(mfrow = c(1, 2))
# plot(hc, hang = -1, main = "Original Tree")
# plot(hc1, hang = -1, main = "Re-start from 10 clusters")
# par(mfrow = c(1, 1))


#### Classical (Metric) Multidimensional Scaling ##############################

cmd_test = cmdscale(dist_test, k = 2, eig = TRUE)
x = cmd_test$points[,1]
y = cmd_test$points[,2]
plot(x, y, type = "n", xlab = "", ylab = "", asp = 1, axes = FALSE,
     main = "cmdscale(bs)")
text(x, y, rownames(test), cex = 0.6)



#### Partitioning (Clustering) Around Medoids #################################

pam_test = cluster::pam(dist_test, k = 2, metric = distance_metrics[1])
plot(pam_test)
