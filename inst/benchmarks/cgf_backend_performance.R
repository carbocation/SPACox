library(SPACox)

args = commandArgs(trailingOnly=TRUE)
n = if(length(args) >= 1) as.integer(args[1]) else 500000
length.out = if(length(args) >= 2) as.integer(args[2]) else 501
threads = if(length(args) >= 3) as.integer(args[3]) else NULL

set.seed(2026)
n.zero = min(1000, floor(n/10))
residuals = c(1-rexp(n-n.zero), rep(0, n.zero))
empirical.cgf = getFromNamespace("SPACox_empirical_CGF", "SPACox")

time.backend = function(backend, threads=NULL) {
  gc(FALSE)
  elapsed = system.time({
    value = empirical.cgf(
      residuals,
      c(-100, 100),
      length.out,
      backend=backend,
      threads=threads
    )
  })[["elapsed"]]
  list(value=value, elapsed=unname(elapsed))
}

reference = time.backend("R")
native = time.backend("rust", threads)
difference = abs(reference$value$cumul-native$value$cumul)

cat("residuals:", n, "\n")
cat("zero residuals:", n.zero, "\n")
cat("CGF grid points:", length.out, "\n")
cat("Rust threads:", if(is.null(threads)) "auto" else threads, "\n")
cat("R elapsed:", reference$elapsed, "seconds\n")
cat("Rust elapsed:", native$elapsed, "seconds\n")
cat("speedup:", reference$elapsed/native$elapsed, "x\n")
cat("maximum absolute difference:", max(difference), "\n")
cat(
  "maximum scaled difference:",
  max(difference/pmax(1, abs(reference$value$cumul))),
  "\n"
)
