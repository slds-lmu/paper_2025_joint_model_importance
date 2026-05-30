# Diagnostic: why are some Rashomon sets so small?
#
# Compares res_dt and vic per (task, learner): model counts and best scores.
# If vic is missing the best models for a task, the RS will shrink even though
# many models in res_dt sit within epsilon of the global best.

source("init.R")
library(data.table)

RS_epsilon <- 0.01

cat("\n=== Per-task model counts (res_dt vs vic) ===\n\n")
res_counts <- as.matrix(table(res_dt$task, res_dt$learner))
vic_counts <- matrix(0L, nrow = length(vic), ncol = ncol(res_counts),
                     dimnames = list(names(vic), colnames(res_counts)))
for (t in names(vic)) {
  cnames <- colnames(vic[[t]])
  if (length(cnames) <= 1) next
  learners <- sub(".*_(.*?)_.*", "\\1", cnames[-1])
  tbl <- table(learners)
  for (l in names(tbl)) if (l %in% colnames(vic_counts)) vic_counts[t, l] <- tbl[[l]]
}
res_aligned <- res_counts[rownames(vic_counts), colnames(vic_counts), drop = FALSE]
cat("res_dt counts:\n");  print(res_aligned)
cat("\nvic counts:\n");   print(vic_counts)
cat("\nDiff (res_dt - vic) -- positive means models in res_dt but no VIC:\n")
print(res_aligned - vic_counts)

cat("\n=== Per-task best score: full res_dt vs vic-only ===\n\n")
performance <- build_performance_from_res_dt(res_dt,
                                              tasks    = names(vic),
                                              learners = learner.keys)
for (t in names(vic)) {
  best_res <- suppressWarnings(min(res_dt[task == t, test.score]))
  vic_cnames <- colnames(vic[[t]])[-1]
  if (length(vic_cnames) == 0) {
    cat(sprintf("  %-8s  res=%.5f  vic=NA\n", t, best_res))
    next
  }
  vic_learner  <- sub(".*_(.*?)_.*", "\\1", vic_cnames)
  vic_model_no <- as.integer(sub(".*_m([0-9]+)$", "\\1", vic_cnames))
  scores <- mapply(function(l, m) {
    p <- performance[[t]][[l]]
    if (is.null(p) || length(p) == 0 || is.na(m) || m > length(p)) return(NA_real_)
    p[m]
  }, vic_learner, vic_model_no)
  best_vic <- suppressWarnings(min(scores, na.rm = TRUE))
  threshold <- best_res * (1 + RS_epsilon)
  rs_size <- sum(scores < threshold, na.rm = TRUE)
  cat(sprintf("  %-8s  res_best=%.5f  vic_best=%.5f  threshold(res)=%.5f  RS=%d\n",
              t, best_res, best_vic, threshold, rs_size))
}

cat("\n=== eps=1%% RS sizes per (task, learner), counted from res_dt only ===\n")
cat("(this is what RS would be if every res_dt model had VIC)\n\n")
res_best <- res_dt[, .(best = min(test.score)), by = task]
res_dt_x <- merge(res_dt, res_best, by = "task")
rs_full <- res_dt_x[test.score < best * (1 + RS_epsilon),
                    .N, by = .(task, learner)]
print(dcast(rs_full, task ~ learner, value.var = "N", fill = 0))

cat("\nDone.\n")
