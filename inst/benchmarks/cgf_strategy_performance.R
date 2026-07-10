library(SPACox)

args = commandArgs(trailingOnly=TRUE)
n = if(length(args) >= 1) as.integer(args[1]) else 500000L
length.out = if(length(args) >= 2) as.integer(args[2]) else 10000L
threads = if(length(args) >= 3) as.integer(args[3]) else NULL
direct.points = if(length(args) >= 4) as.integer(args[4]) else 300L

set.seed(2026)
n.zero = min(1000L, floor(n/10))
residuals = c(1-rexp(n-n.zero), numeric(n.zero))
lazy.cgf = getFromNamespace("SPACox_lazy_CGF", "SPACox")
eager.cgf = getFromNamespace("SPACox_empirical_CGF", "SPACox")

lazy.setup.elapsed = system.time({
  lazy = lazy.cgf(
    residuals,
    c(-100, 100),
    length.out,
    backend="rust",
    threads=threads
  )
})[["elapsed"]]

# Two requested values per call approximate a binary grouped genotype during
# root finding while ensuring each batch is distinct and therefore uncached.
t.values = seq(-0.03, 0.03, length.out=direct.points)
t.batches = split(t.values, ceiling(seq_along(t.values)/2))
direct.elapsed = system.time({
  for(t in t.batches)
    invisible(lazy$K_1_emp(t))
})[["elapsed"]]

eager.elapsed = system.time({
  eager = eager.cgf(
    residuals,
    c(-100, 100),
    length.out,
    backend="rust",
    threads=threads
  )
})[["elapsed"]]

lazy.elapsed = lazy.setup.elapsed + direct.elapsed
cat("residuals:", n, "\n")
cat("CGF grid points:", length.out, "\n")
cat("simulated exact points:", direct.points, "\n")
cat("Rust threads:", if(is.null(threads)) "auto" else threads, "\n")
cat("lazy setup elapsed:", lazy.setup.elapsed, "seconds\n")
cat("lazy direct-use elapsed:", direct.elapsed, "seconds\n")
cat("lazy total elapsed:", lazy.elapsed, "seconds\n")
cat("lazy final mode:", lazy$state$mode, "\n")
cat("lazy direct points:", lazy$state$direct_points, "\n")
cat("lazy grid builds:", lazy$state$grid_builds, "\n")
cat("eager elapsed:", eager.elapsed, "seconds\n")
cat("strategy speedup:", eager.elapsed/lazy.elapsed, "x\n")
