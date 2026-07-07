library(survival)

if(file.exists("R/Library.r")) {
  source("R/Library.r")
} else {
  library(SPACox)
}

assert_true = function(x, msg) {
  if(!isTRUE(x)) stop(msg, call. = FALSE)
}

set.seed(11)

n.subjects = 6
n.intervals = 3

dat = data.frame(
  id = rep(seq_len(n.subjects), each=n.intervals),
  start = rep(0:(n.intervals-1), n.subjects),
  stop = rep(1:n.intervals, n.subjects),
  event = 0,
  x = rnorm(n.subjects*n.intervals)
)
dat$event[seq(n.intervals, n.subjects*n.intervals, by=n.intervals)] = c(1, 0, 1, 0, 1, 0)

g.subject = c(2, 2, 1, 0, 2, 1)
geno = matrix(g.subject, ncol=1)
rownames(geno) = as.character(seq_len(n.subjects))
colnames(geno) = "snp1"

obj.null = SPACox_Null_Model(
  Surv(start, stop, event) ~ x,
  data=dat,
  pIDs=as.character(dat$id),
  gIDs=rownames(geno),
  length.out=200
)

g.expanded = g.subject[obj.null$p2g]
assert_true(length(g.expanded) == length(obj.null$resid), "expanded genotype should match residual length")
assert_true(
  length(g.expanded) - length(which(g.expanded != 0)) >= 0,
  "expanded genotype zero count should not be negative"
)

obj.expanded = obj.null
obj.expanded$p2g = NULL

interval.result = SPACox.one.SNP(
  g.subject,
  obj.null,
  Cutoff=0,
  CovAdj.cutoff=0
)
expanded.result = SPACox.one.SNP(
  g.expanded,
  obj.expanded,
  Cutoff=0,
  CovAdj.cutoff=0
)

assert_true(length(interval.result) == 7, "interval result should have seven output values")
assert_true(all(is.finite(interval.result[5:7])), "interval statistic, variance, and z should be finite")
assert_true(
  isTRUE(all.equal(interval.result, expanded.result, tolerance=1e-10, check.attributes=FALSE)),
  "repeated-ID interval result should match equivalent expanded-genotype result"
)

right.dat = data.frame(
  id = seq_len(n.subjects),
  time = c(1, 2, 3, 4, 5, 6),
  event = c(1, 0, 1, 0, 1, 0),
  x = rnorm(n.subjects)
)
right.obj = SPACox_Null_Model(
  Surv(time, event) ~ x,
  data=right.dat,
  pIDs=as.character(right.dat$id),
  gIDs=rownames(geno),
  length.out=200
)
right.result = SPACox.one.SNP(g.subject, right.obj, Cutoff=0, CovAdj.cutoff=0)

assert_true(length(right.result) == 7, "right-censored result should have seven output values")
assert_true(all(is.finite(right.result[5:7])), "right-censored statistic, variance, and z should be finite")
