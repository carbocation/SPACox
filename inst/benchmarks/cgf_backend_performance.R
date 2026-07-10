library(SPACox)

args = commandArgs(trailingOnly=TRUE)
n = if(length(args) >= 1) as.integer(args[1]) else 500000
length.out = if(length(args) >= 2) as.integer(args[2]) else 501
threads = if(length(args) >= 3) as.integer(args[3]) else NULL
run.reference = if(length(args) >= 4) as.logical(args[4]) else TRUE

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

reference = if(run.reference) time.backend("R") else NULL
scalar = time.backend("rust-scalar", threads)
simd = time.backend("rust", threads)
scalar.difference = if(run.reference) {
  abs(reference$value$cumul-scalar$value$cumul)
} else {
  NULL
}
simd.difference = abs(scalar$value$cumul-simd$value$cumul)

cat("residuals:", n, "\n")
cat("zero residuals:", n.zero, "\n")
cat("CGF grid points:", length.out, "\n")
cat("Rust threads:", if(is.null(threads)) "auto" else threads, "\n")
if(run.reference)
  cat("R elapsed:", reference$elapsed, "seconds\n")
cat("Scalar Rust elapsed:", scalar$elapsed, "seconds\n")
cat("Selected Rust kernel:", simd$value$kernel, "\n")
cat("Selected Rust elapsed:", simd$elapsed, "seconds\n")
if(run.reference)
  cat("Scalar Rust speedup over R:", reference$elapsed/scalar$elapsed, "x\n")
cat("Selected-kernel speedup over scalar Rust:", scalar$elapsed/simd$elapsed, "x\n")
if(run.reference) {
  cat("maximum scalar-vs-R absolute difference:", max(scalar.difference), "\n")
  cat(
    "maximum scalar-vs-R scaled difference:",
    max(scalar.difference/pmax(1, abs(reference$value$cumul))),
    "\n"
  )
}
cat("maximum SIMD-vs-scalar absolute difference:", max(simd.difference), "\n")
cat(
  "maximum SIMD-vs-scalar scaled difference:",
  max(simd.difference/pmax(1, abs(scalar$value$cumul))),
  "\n"
)
