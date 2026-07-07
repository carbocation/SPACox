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
SPACox_group_values = get_spacox_fun("SPACox_group_values")
GetProb_SPA = get_spacox_fun("GetProb_SPA")
GetProb_SPA_grouped = get_spacox_fun("GetProb_SPA_grouped")

resid.test = c(0, 0, -0.75, 0, 0.2, 0, 1.1, -0.4, 0)
range.test = c(-3, 3)
length.out.test = 101
new.cgf = SPACox_empirical_CGF(resid.test, range.test, length.out.test)
old.cgf = old_empirical_cgf(resid.test, range.test, length.out.test)

assert_true(
  isTRUE(all.equal(new.cgf$cumul, old.cgf, tolerance=1e-12, check.attributes=FALSE)),
  "zero-compressed empirical CGF should match the old scalar formula"
)
assert_true(new.cgf$n_total == length(resid.test), "CGF total count should match residual length")
assert_true(new.cgf$n_zero == sum(resid.test == 0), "CGF zero count should match exact zero residual count")

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
  length.out=200
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
