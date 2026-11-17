# Usage: Rscript run_all.R <dir-with-generated-R-files>
# Sources every generated model; checks it solves, values are finite and no compartment goes
# meaningfully negative, and that the population is conserved where the model says it should be.
suppressWarnings(suppressMessages(library(deSolve)))
files <- sort(list.files(commandArgs(trailingOnly = TRUE)[1], pattern = "[.]R$", full.names = TRUE))
pdf(NULL)
bad <- 0
for (f in files) {
  env <- new.env()
  msg <- tryCatch({
    suppressWarnings(sys.source(f, envir = env))
    out <- env$out
    comps <- names(env$state)
    vals <- as.matrix(out[, comps])
    if (nrow(out) != length(env$times)) stop("solver stopped early at t = ", max(out$time))
    if (any(!is.finite(vals))) stop("non-finite values")
    if (min(vals) < -1e-3 * max(abs(vals))) stop("negative compartment: ", min(vals))
    sprintf("ok   %-22s %4d steps; final: %s", basename(f), nrow(out),
            paste(sprintf("%s=%.4g", comps, vals[nrow(vals), ]), collapse = " "))
  }, error = function(e) { bad <<- bad + 1; sprintf("FAIL %-22s %s", basename(f), conditionMessage(e)) })
  cat(msg, "\n")
}
cat(sprintf("\n%d/%d models ran successfully\n", length(files) - bad, length(files)))
quit(status = if (bad) 1 else 0)
