library(survival)

if(file.exists("R/Library.r")) {
  source("R/Library.r")
} else {
  library(SPACox)
}

assert_true = function(x, msg) {
  if(!isTRUE(x)) stop(msg, call. = FALSE)
}

get_spacox_fun = function(name) {
  if(exists(name, mode="function", inherits=TRUE)) {
    get(name, mode="function", inherits=TRUE)
  } else {
    getFromNamespace(name, "SPACox")
  }
}

old_empirical_cgf = function(mresid, range, length.out) {
  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)
  cumul = matrix(NA_real_, length.out, 4)

  for(i in seq_along(idx1)){
    t = idx1[i]
    e_resid = exp(mresid*t)
    M0 = mean(e_resid)
    M1 = mean(mresid*e_resid)
    M2 = mean(mresid^2*e_resid)
    K0 = log(M0)
    K1 = M1/M0
    K2 = (M0*M2-M1^2)/M0^2
    cumul[i,] = c(t, K0, K1, K2)
  }

  cumul
}

SPACox_empirical_CGF = get_spacox_fun("SPACox_empirical_CGF")
SPACox_lazy_CGF = get_spacox_fun("SPACox_lazy_CGF")
SPACox_lazy_CGF_prepare_workload = get_spacox_fun("SPACox_lazy_CGF_prepare_workload")
SPACox_prepare_analysis_CGF = get_spacox_fun("SPACox_prepare_analysis_CGF")
SPACox_empirical_CGF_R = get_spacox_fun("SPACox_empirical_CGF_R")
SPACox_group_values = get_spacox_fun("SPACox_group_values")
GetProb_SPA = get_spacox_fun("GetProb_SPA")
GetProb_SPA_grouped = get_spacox_fun("GetProb_SPA_grouped")

resid.test = c(0, 0, -0.75, 0, 0.2, 0, 1.1, -0.4, 0)
range.test = c(-3, 3)
length.out.test = 101
new.cgf = SPACox_empirical_CGF(resid.test, range.test, length.out.test, backend="R")
old.cgf = old_empirical_cgf(resid.test, range.test, length.out.test)

assert_true(
  isTRUE(all.equal(new.cgf$cumul, old.cgf, tolerance=1e-12, check.attributes=FALSE)),
  "zero-compressed empirical CGF should match the old scalar formula"
)
assert_true(new.cgf$n_total == length(resid.test), "CGF total count should match residual length")
assert_true(new.cgf$n_zero == sum(resid.test == 0), "CGF zero count should match exact zero residual count")

overflow.resid = c(-8, 0, 0.5, 1)
overflow.cgf = SPACox_empirical_CGF(overflow.resid, c(-100, 100), 101, backend="R")

assert_true(
  all(is.finite(overflow.cgf$cumul)),
  "empirical CGF should remain finite when unshifted exponentials would overflow"
)
assert_true(
  isTRUE(all.equal(overflow.cgf$cumul[1,2], 800-log(length(overflow.resid)), tolerance=1e-12)),
  "log-space empirical CGF should preserve the dominant tail term"
)
assert_true(
  isTRUE(all.equal(overflow.cgf$cumul[1,3], -8, tolerance=1e-12)),
  "log-space empirical CGF derivative should approach the extreme residual"
)
assert_true(
  all(overflow.cgf$cumul[,4] >= 0),
  "empirical CGF second derivative should be non-negative"
)

rust.available = exists(
  "C_spacox_empirical_cgf",
  envir=environment(SPACox_empirical_CGF),
  inherits=TRUE
)

if(rust.available) {
  rust.small.scalar = SPACox_empirical_CGF(
    resid.test,
    range.test,
    length.out.test,
    backend="rust-scalar",
    threads=1
  )
  rust.small = SPACox_empirical_CGF(
    resid.test,
    range.test,
    length.out.test,
    backend="rust",
    threads=1
  )
  assert_true(
    isTRUE(all.equal(rust.small.scalar$cumul, new.cgf$cumul, tolerance=1e-12, check.attributes=FALSE)),
    "Scalar Rust CGF should match the R reference implementation"
  )
  assert_true(
    isTRUE(all.equal(rust.small$cumul, new.cgf$cumul, tolerance=1e-12, check.attributes=FALSE)),
    "SIMD Rust CGF should match the R reference implementation"
  )
  assert_true(
    rust.small$kernel %in% c("avx2", "scalar"),
    "Rust CGF should report the selected native kernel"
  )

  rust.overflow = SPACox_empirical_CGF(
    overflow.resid,
    c(-100, 100),
    101,
    backend="rust",
    threads=1
  )
  assert_true(
    isTRUE(all.equal(rust.overflow$cumul, overflow.cgf$cumul, tolerance=1e-12, check.attributes=FALSE)),
    "SIMD Rust CGF should match the R reference implementation in overflow tails"
  )

  set.seed(22)
  rust.resid = c(rnorm(10000), rep(0, 1000))
  rust.reference = SPACox_empirical_CGF(
    rust.resid,
    c(-100, 100),
    513,
    backend="R"
  )
  rust.single = SPACox_empirical_CGF(
    rust.resid,
    c(-100, 100),
    513,
    backend="rust",
    threads=1
  )
  rust.parallel = SPACox_empirical_CGF(
    rust.resid,
    c(-100, 100),
    513,
    backend="rust",
    threads=2
  )
  rust.scalar = SPACox_empirical_CGF(
    rust.resid,
    c(-100, 100),
    513,
    backend="rust-scalar",
    threads=1
  )
  assert_true(
    isTRUE(all.equal(rust.single$cumul, rust.reference$cumul, tolerance=1e-10, check.attributes=FALSE)),
    "Rust CGF should match the R reference for random residuals"
  )
  assert_true(
    isTRUE(all.equal(rust.scalar$cumul, rust.reference$cumul, tolerance=1e-10, check.attributes=FALSE)),
    "Scalar Rust CGF should match the R reference for random residuals"
  )
  assert_true(
    isTRUE(all.equal(rust.single$cumul, rust.scalar$cumul, tolerance=1e-10, check.attributes=FALSE)),
    "SIMD Rust CGF should match the scalar native kernel"
  )
  assert_true(
    identical(rust.single$cumul, rust.parallel$cumul),
    "Rust CGF should be deterministic across thread counts"
  )
}

set.seed(23)
lazy.resid = c(rnorm(1000), rep(0, 20))
lazy.cgf = SPACox_lazy_CGF(
  lazy.resid,
  c(-100, 100),
  10000,
  backend="R"
)
lazy.t = c(-3, -0.25, 0, 0.5, 4)
lazy.reference = SPACox_empirical_CGF_R(
  lazy.resid[lazy.resid != 0],
  sum(lazy.resid == 0),
  length(lazy.resid),
  lazy.t
)
lazy.values = cbind(
  lazy.cgf$K_org_emp(lazy.t),
  lazy.cgf$K_1_emp(lazy.t),
  lazy.cgf$K_2_emp(lazy.t)
)

assert_true(
  isTRUE(all.equal(lazy.values, lazy.reference, tolerance=1e-12, check.attributes=FALSE)),
  "lazy direct CGF evaluation should match the exact reference implementation"
)
assert_true(lazy.cgf$state$mode == "direct", "a sparse lazy CGF request should use direct evaluation")
assert_true(lazy.cgf$state$direct_points == length(lazy.t), "reusing a cached direct request should not add work")
assert_true(lazy.cgf$state$grid_builds == 0, "a sparse lazy CGF request should not build the grid")

dense.t = seq(-1, 1, length.out=257)
dense.values = lazy.cgf$K_1_emp(dense.t)
assert_true(length(dense.values) == length(dense.t), "dense lazy CGF requests should preserve their length")
assert_true(lazy.cgf$state$mode == "hybrid", "a dense request should add a grid without changing sparse direct evaluation")
assert_true(lazy.cgf$state$grid_builds == 1, "the lazy CGF grid should be built once")
invisible(lazy.cgf$K_2_emp(dense.t))
assert_true(lazy.cgf$state$grid_builds == 1, "the lazy CGF grid should be reused")

workload.cgf = SPACox_lazy_CGF(lazy.resid, c(-100, 100), 1000, backend="R")
SPACox_lazy_CGF_prepare_workload(list(cgf_state=workload.cgf$state), 2)
assert_true(workload.cgf$state$sparse_mode == "grid", "a large workload should select grid evaluation")
assert_true(workload.cgf$state$grid_builds == 0, "a selected grid should remain deferred until SPA is needed")
invisible(workload.cgf$K_1_emp(c(-0.1, 0.1)))
assert_true(workload.cgf$state$mode == "grid", "the first SPA request should build a selected grid")
assert_true(workload.cgf$state$grid_builds == 1, "workload fallback should build one grid")

continued.workload.cgf = SPACox_lazy_CGF(lazy.resid, c(-100, 100), 1000, backend="R")
invisible(continued.workload.cgf$K_1_emp(c(-0.1, 0.1)))
SPACox_lazy_CGF_prepare_workload(
  list(cgf_state=continued.workload.cgf$state),
  1000000
)
assert_true(
  continued.workload.cgf$state$sparse_mode == "grid",
  "a dense workload should select the grid even after prior direct evaluation"
)
assert_true(
  continued.workload.cgf$state$workload_selected,
  "dense workload selection should be recorded after prior direct evaluation"
)

upgraded.workload.cgf = SPACox_lazy_CGF(lazy.resid, c(-100, 100), 10000, backend="R")
SPACox_lazy_CGF_prepare_workload(list(cgf_state=upgraded.workload.cgf$state), 1)
assert_true(
  upgraded.workload.cgf$state$sparse_mode == "direct",
  "a small initial matrix workload should retain direct evaluation"
)
SPACox_lazy_CGF_prepare_workload(list(cgf_state=upgraded.workload.cgf$state), 1000000)
assert_true(
  upgraded.workload.cgf$state$sparse_mode == "grid",
  "a later dense workload should upgrade a prior small matrix workload to the grid"
)

reversed.range.cgf = SPACox_lazy_CGF(
  lazy.resid,
  c(100, -100),
  1000,
  backend="R"
)
forward.range.cgf = SPACox_lazy_CGF(
  lazy.resid,
  c(-100, 100),
  1000,
  backend="R"
)
range.test.points = c(-1, 0, 1)
assert_true(
  isTRUE(all.equal(
    reversed.range.cgf$K_org_emp(range.test.points),
    forward.range.cgf$K_org_emp(range.test.points),
    tolerance=0,
    check.attributes=FALSE
  )),
  "reversed symmetric ranges should retain the pre-Rust CGF behavior"
)

set.seed(21)
n.subjects = 7
n.intervals = 4
dat = data.frame(
  id = rep(seq_len(n.subjects), each=n.intervals),
  start = rep(0:(n.intervals-1), n.subjects),
  stop = rep(1:n.intervals, n.subjects),
  event = 0,
  x = rnorm(n.subjects*n.intervals)
)
dat$event[seq(n.intervals, n.subjects*n.intervals, by=n.intervals)] = c(1, 0, 1, 0, 1, 0, 1)

g.subject = c(2, 2, 1, 0, 2, 1, 0)
geno = matrix(g.subject, ncol=1)
rownames(geno) = as.character(seq_len(n.subjects))
colnames(geno) = "snp1"

obj.null = SPACox_Null_Model(
  Surv(start, stop, event) ~ x,
  data=dat,
  pIDs=as.character(dat$id),
  gIDs=rownames(geno),
  length.out=10000
)

assert_true(obj.null$cgf_strategy == "lazy", "SPACox null models should use lazy CGF evaluation by default")
assert_true(obj.null$cgf_state$mode == "deferred", "fitting a lazy null model should not calculate the CGF")
assert_true(!is.environment(obj.null$cgf_state), "the public null model should not expose mutable CGF state")

obj.null.alias = obj.null
prepared.alias = SPACox_prepare_analysis_CGF(obj.null.alias, 1000)
assert_true(is.environment(prepared.alias$cgf_state), "an analysis should receive private mutable CGF state")
assert_true(prepared.alias$cgf_state$sparse_mode == "grid", "a private dense analysis should select the grid")
assert_true(obj.null$cgf_state$mode == "deferred", "preparing an alias should not mutate the original null model")
assert_true(obj.null.alias$cgf_state$mode == "deferred", "preparing an alias should not mutate the copied null model")
invisible(prepared.alias$K_1_emp(c(-0.1, 0.1)))
assert_true(prepared.alias$cgf_state$grid_builds == 1, "the private analysis should cache its grid")
assert_true(obj.null$cgf_state$grid_builds == 0, "private grid caching should not leak to the original null model")
reused.alias = SPACox_prepare_analysis_CGF(prepared.alias, 1)
assert_true(
  identical(reused.alias$cgf_state, prepared.alias$cgf_state),
  "chunked analyses should reuse the same private evaluator"
)

normal.only = SPACox(
  obj.null,
  geno,
  Cutoff=Inf,
  min.maf=0,
  missing.cutoff=1
)
assert_true(nrow(normal.only) == ncol(geno), "normal-only SPACox should return one row per variant")
assert_true(obj.null$cgf_state$mode == "deferred", "normal-only variants should not calculate the CGF")

spa.matrix = SPACox(
  obj.null.alias,
  geno,
  Cutoff=0,
  min.maf=0,
  missing.cutoff=1
)
spa.single = SPACox.one.SNP(
  g.subject,
  obj.null,
  Cutoff=0,
  min.maf=0,
  missing.cutoff=1
)
assert_true(
  isTRUE(all.equal(
    as.numeric(spa.matrix[1,]),
    spa.single,
    tolerance=1e-10,
    check.attributes=FALSE
  )),
  "matrix and one-SNP analyses should agree with private lazy evaluators"
)
assert_true(obj.null$cgf_state$mode == "deferred", "one-SNP analysis should leave the original model immutable")
assert_true(obj.null.alias$cgf_state$mode == "deferred", "matrix analysis should leave an aliased model immutable")

assert_true(
  is.null(obj.null$obj.coxph$y),
  "SPACox null models should not retain the Cox response matrix by default"
)

obj.null.with.y = SPACox_Null_Model(
  Surv(start, stop, event) ~ x,
  data=dat,
  pIDs=as.character(dat$id),
  gIDs=rownames(geno),
  length.out=20,
  cgf.strategy="eager",
  y=TRUE
)

assert_true(
  !is.null(obj.null.with.y$obj.coxph$y),
  "SPACox null models should retain the Cox response matrix when y=TRUE"
)
assert_true(obj.null.with.y$cgf_strategy == "eager", "eager CGF evaluation should remain available")
assert_true(is.null(obj.null.with.y$cgf_state), "eager CGF evaluation should not create lazy state")

unit.weights = rep(1, nrow(dat))
positional.arguments = list(
  Surv(start, stop, event) ~ x,
  dat,
  as.character(dat$id),
  rownames(geno),
  c(-100, 100),
  20,
  unit.weights
)
positional.arguments$cgf.backend = "R"
positional.arguments$cgf.strategy = "eager"
obj.null.positional.weights = do.call(SPACox_Null_Model, positional.arguments)

named.arguments = positional.arguments[1:6]
names(named.arguments) = c("formula", "data", "pIDs", "gIDs", "range", "length.out")
named.arguments$weights = unit.weights
named.arguments$cgf.backend = "R"
named.arguments$cgf.strategy = "eager"
obj.null.named.weights = do.call(SPACox_Null_Model, named.arguments)
assert_true(
  isTRUE(all.equal(
    obj.null.positional.weights$resid,
    obj.null.named.weights$resid,
    tolerance=0,
    check.attributes=FALSE
  )),
  "positional coxph arguments should continue to pass through the null-model wrapper"
)

g.row = g.subject[obj.null$row_to_genotype]
MAF = mean(g.subject, na.rm=TRUE)/2
S.grouped = sum(g.subject * obj.null$resid_sum_by_genotype)
S.row = sum(g.row * obj.null$resid)
S.var.grouped = obj.null$var.resid * sum(obj.null$row_count_by_genotype * (g.subject - 2*MAF)^2)
S.var.row = obj.null$var.resid * sum((g.row - 2*MAF)^2)

assert_true(isTRUE(all.equal(S.grouped, S.row, tolerance=1e-12)), "grouped score should match row-expanded score")
assert_true(isTRUE(all.equal(S.var.grouped, S.var.row, tolerance=1e-12)), "grouped variance should match row-expanded variance")
assert_true(
  isTRUE(all.equal(as.vector(obj.null$tX_by_genotype %*% g.subject),
                   as.vector(obj.null$tX %*% g.row),
                   tolerance=1e-12)),
  "grouped tX cache should match row-expanded tX product"
)

G1.row = g.row - 2*MAF
G1norm.row = G1.row/sqrt(S.var.row)
z1 = S.row/sqrt(S.var.row)
N1set = which(g.row != 0)
N0 = length(g.row) - length(N1set)
G1N1 = G1norm.row[N1set]
G1N0 = -2*MAF/sqrt(S.var.row)

G1norm.subject = (g.subject - 2*MAF)/sqrt(S.var.grouped)
G1.grouped = SPACox_group_values(G1norm.subject, obj.null$row_count_by_genotype)

p.direct.upper = GetProb_SPA_grouped(obj.null, G1.grouped$values, G1.grouped$counts, abs(z1), lower.tail=FALSE)
p.direct.lower = GetProb_SPA_grouped(obj.null, G1.grouped$values, G1.grouped$counts, -abs(z1), lower.tail=TRUE)
assert_true(obj.null$cgf_state$mode == "deferred", "direct helper evaluation should not mutate the public null model")
assert_true(obj.null$cgf_state$grid_builds == 0, "direct helper evaluation should not cache state on the public null model")

default.grid.cgf = SPACox_empirical_CGF(obj.null$resid, c(-100, 100), 10000, backend="R")
obj.null.default.grid = obj.null
obj.null.default.grid$K_org_emp = default.grid.cgf$K_org_emp
obj.null.default.grid$K_1_emp = default.grid.cgf$K_1_emp
obj.null.default.grid$K_2_emp = default.grid.cgf$K_2_emp
p.default.grid.upper = GetProb_SPA_grouped(
  obj.null.default.grid,
  G1.grouped$values,
  G1.grouped$counts,
  abs(z1),
  lower.tail=FALSE
)
p.default.grid.lower = GetProb_SPA_grouped(
  obj.null.default.grid,
  G1.grouped$values,
  G1.grouped$counts,
  -abs(z1),
  lower.tail=TRUE
)

p.old.upper = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, abs(z1), lower.tail=FALSE)
p.old.lower = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, -abs(z1), lower.tail=TRUE)
p.new.upper = GetProb_SPA_grouped(obj.null, G1.grouped$values, G1.grouped$counts, abs(z1), lower.tail=FALSE)
p.new.lower = GetProb_SPA_grouped(obj.null, G1.grouped$values, G1.grouped$counts, -abs(z1), lower.tail=TRUE)

assert_true(
  isTRUE(all.equal(p.new.upper, p.old.upper, tolerance=1e-10, check.attributes=FALSE)),
  "grouped upper-tail SPA should match ungrouped SPA"
)
assert_true(
  isTRUE(all.equal(p.new.lower, p.old.lower, tolerance=1e-10, check.attributes=FALSE)),
  "grouped lower-tail SPA should match ungrouped SPA"
)
assert_true(
  isTRUE(all.equal(p.direct.upper, p.default.grid.upper, tolerance=1e-5, check.attributes=FALSE)),
  "direct upper-tail SPA should agree with the default interpolation grid"
)
assert_true(
  isTRUE(all.equal(p.direct.lower, p.default.grid.lower, tolerance=1e-5, check.attributes=FALSE)),
  "direct lower-tail SPA should agree with the default interpolation grid"
)
