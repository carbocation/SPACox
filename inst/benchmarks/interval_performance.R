library(survival)

if(file.exists("R/Library.r")) {
  source("R/Library.r")
} else {
  library(SPACox)
}

args = commandArgs(trailingOnly=TRUE)
n.subjects = if(length(args) >= 1) as.integer(args[1]) else 2000
n.intervals = if(length(args) >= 2) as.integer(args[2]) else 8
n.snps = if(length(args) >= 3) as.integer(args[3]) else 20
length.out = if(length(args) >= 4) as.integer(args[4]) else 500

set.seed(101)

id = rep(seq_len(n.subjects), each=n.intervals)
start = rep(0:(n.intervals-1), n.subjects)
stop = start + 1
event = rep(0, n.subjects*n.intervals)
event[seq(n.intervals, n.subjects*n.intervals, by=n.intervals)] = rbinom(n.subjects, 1, 0.08)

phenotype = data.frame(
  id = id,
  start = start,
  stop = stop,
  event = event,
  x1 = rep(rnorm(n.subjects), each=n.intervals),
  x2 = rep(rbinom(n.subjects, 1, 0.4), each=n.intervals)
)

genotype = matrix(rbinom(n.subjects*n.snps, 2, 0.1), n.subjects, n.snps)
rownames(genotype) = as.character(seq_len(n.subjects))
colnames(genotype) = paste0("snp", seq_len(n.snps))

coxph.time = system.time({
  coxph.fit = coxph(Surv(start, stop, event) ~ x1 + x2, data=phenotype, x=TRUE)
})
null.time = system.time({
  obj.null = SPACox_Null_Model(
    Surv(start, stop, event) ~ x1 + x2,
    data=phenotype,
    pIDs=as.character(phenotype$id),
    gIDs=rownames(genotype),
    length.out=length.out
  )
})
spacox.time = system.time({
  spacox.fit = SPACox(obj.null, genotype)
})

cat("subjects:", n.subjects, "\n")
cat("interval rows:", nrow(phenotype), "\n")
cat("snps:", n.snps, "\n")
cat("length.out:", length.out, "\n")
cat("zero residual rows:", sum(obj.null$resid == 0), "\n")
cat("coxph elapsed:", unname(coxph.time["elapsed"]), "\n")
cat("SPACox_Null_Model elapsed:", unname(null.time["elapsed"]), "\n")
cat("SPACox elapsed:", unname(spacox.time["elapsed"]), "\n")
