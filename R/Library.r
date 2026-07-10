#' Fits a NULL model for SPACox
#'
#' Fits a null Cox proportional hazards model and then calculates the empirical cumulant generation function (CGF) of the martingale residuals
#' @param formula a formula to be passed to function coxph(). For more details, please refer to package survival.
#' @param data a data.frame in which to interpret the variables named in the formula
#' @param pIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order in the formula.
#' @param gIDs a character vector of subject IDs. NOTE: its order should be the same as the subjects order of the Geno.mtx (i.e. the input of the function SPACox()).
#' @param range a two-element numeric vector (default: c(-100,100)) to specify the domain of the empirical CGF.
#' @param length.out a positive integer (default: 10000) for empirical CGF. Larger length.out corresponds to longer calculation time and more accurate estimated empirical CGF.
#' @param cgf.backend implementation used to calculate the empirical CGF. The default "rust" uses AVX2 SIMD when available and otherwise falls back to the original scalar native kernel; "rust-scalar" forces the scalar native kernel; "R" uses the reference implementation.
#' @param cgf.threads number of threads used by the Rust CGF backend. NULL uses the number of threads available to the process.
#' @param y whether to retain the response matrix in the fitted coxph object. The default FALSE reduces the memory retained by the null model.
#' @param ... Other arguments passed to function coxph(). For more details, please refer to package survival.
#' @return an object with a class of "SPACox_NULL_Model".
#' @examples
#' # Please check help(SPACox) for a simulated example.
#' @export
#' @import survival
#' @useDynLib SPACox, .registration=TRUE, .fixes="C_"
SPACox_Null_Model = function(formula,
                             data=NULL,
                             pIDs=NULL,
                             gIDs=NULL,
                             range=c(-100,100),
                             length.out = 10000,
                             cgf.backend = c("rust", "rust-scalar", "R"),
                             cgf.threads = NULL,
                             y = FALSE,
                             ...)
{
  Call = match.call()

  cgf.backend = match.arg(cgf.backend)

  ### Fit a Cox model
  obj.coxph = coxph(formula, data=data, x=TRUE, y=y, ...)

  ### Check input arguments
  obj.check = check_input(pIDs, gIDs, obj.coxph, range)
  p2g = obj.check$p2g
  pIDs = obj.check$pIDs

  ### Get the covariate matrix to adjust for genotype
  mresid = obj.coxph$residuals
  Cova = obj.coxph$x

  X = cbind(1, Cova)
  X.invXX = X %*% solve(t(X)%*%X)
  tX = t(X)

  ### calculate empirical CGF for martingale residuals
  print("Start calculating empirical CGF for martingale residuals...")
  cgf = SPACox_empirical_CGF(mresid, range, length.out,
                             backend=cgf.backend,
                             threads=cgf.threads)

  var.resid = var(mresid)
  row_to_genotype = if(is.null(p2g)) seq_along(pIDs) else as.integer(p2g)
  row_count_by_genotype = tabulate(row_to_genotype, nbins=length(gIDs))
  resid_sum_by_genotype = SPACox_rowsum_vector(mresid, row_to_genotype, length(gIDs))
  tX_by_genotype = SPACox_rowsum_columns(tX, row_to_genotype, length(gIDs))

  re=list(resid = mresid,
          var.resid = var.resid,
          K_org_emp = cgf$K_org_emp,
          K_1_emp = cgf$K_1_emp,
          K_2_emp = cgf$K_2_emp,
          cgf_n_total = cgf$n_total,
          cgf_n_zero = cgf$n_zero,
          cgf_resid_nonzero = cgf$resid_nonzero,
          cgf_backend = cgf$backend,
          cgf_kernel = cgf$kernel,
          cgf_threads = cgf$threads,
          Call = Call,
          obj.coxph = obj.coxph,
          tX = tX,
          X.invXX = X.invXX,
          p2g = p2g,
          row_to_genotype = row_to_genotype,
          row_count_by_genotype = row_count_by_genotype,
          resid_sum_by_genotype = resid_sum_by_genotype,
          tX_by_genotype = tX_by_genotype,
          gIDs = gIDs,
          pIDs = pIDs)

  class(re)<-"SPACox_NULL_Model"
  return(re)
}

SPACox_empirical_CGF = function(mresid,
                                range,
                                length.out,
                                backend=c("rust", "rust-scalar", "R"),
                                threads=NULL)
{
  backend = match.arg(backend)
  threads = SPACox_CGF_threads(threads)

  idx0 = qcauchy(1:length.out/(length.out+1))
  idx1 = idx0 * max(range) / max(idx0)

  cumul = matrix(NA_real_, length.out, 4)
  cumul[,1] = idx1

  n.total = length(mresid)
  resid_nonzero = mresid[mresid != 0]
  n.zero = n.total - length(resid_nonzero)
  kernel = if(backend == "R") "R" else "scalar"

  if(backend == "rust"){
    simd.binding.available = exists("C_spacox_cgf_simd_available",
                                    envir=environment(SPACox_empirical_CGF),
                                    inherits=FALSE)
    if(simd.binding.available && isTRUE(.Call(C_spacox_cgf_simd_available)))
      kernel = "avx2"
  }

  if(length(resid_nonzero) == 0){
    cumul[,2] = 0
    cumul[,3] = 0
    cumul[,4] = 0
  }else if(backend %in% c("rust", "rust-scalar")){
    rust.available = exists("C_spacox_empirical_cgf",
                            envir=environment(SPACox_empirical_CGF),
                            inherits=FALSE)
    if(!rust.available){
      warning("Rust CGF backend is not loaded; using the R reference backend.")
      backend = "R"
      kernel = "R"
      cumul[,2:4] = SPACox_empirical_CGF_R(resid_nonzero,
                                           n.zero,
                                           n.total,
                                           idx1)
    }else{
      cumul[,2:4] = .Call(C_spacox_empirical_cgf,
                           resid_nonzero,
                           idx1,
                           as.double(n.zero),
                           threads,
                           identical(backend, "rust"))
    }
  }else{
    cumul[,2:4] = SPACox_empirical_CGF_R(resid_nonzero,
                                         n.zero,
                                         n.total,
                                         idx1)
  }

  list(cumul = cumul,
       K_org_emp = approxfun(cumul[,1], cumul[,2], rule=2),
       K_1_emp = approxfun(cumul[,1], cumul[,3], rule=2),
       K_2_emp = approxfun(cumul[,1], cumul[,4], rule=2),
       n_total = n.total,
       n_zero = n.zero,
       resid_nonzero = resid_nonzero,
       backend = backend,
       kernel = kernel,
       threads = if(backend %in% c("rust", "rust-scalar")) threads else 0L)
}

SPACox_CGF_threads = function(threads)
{
  if(is.null(threads))
    return(0L)

  if(length(threads) != 1 ||
     !is.numeric(threads) ||
     !is.finite(threads) ||
     threads < 1 ||
     threads != floor(threads) ||
     threads > .Machine$integer.max)
    stop("cgf.threads should be NULL or a positive integer.")

  as.integer(threads)
}

SPACox_empirical_CGF_R = function(resid_nonzero, n.zero, n.total, idx1)
{
    length.out = length(idx1)
    cumul = matrix(NA_real_, length.out, 3)
    max.outer = 5e6
    chunk.size = max(1, min(256, floor(max.outer/length(resid_nonzero))))
    next.print = 1000
    resid.min = min(resid_nonzero)
    resid.max = max(resid_nonzero)
    resid.squared = resid_nonzero^2

    for(chunk.start in seq(1, length.out, by=chunk.size)){
      chunk.end = min(length.out, chunk.start + chunk.size - 1)
      t.chunk = idx1[chunk.start:chunk.end]

      # Shift each column by its largest log weight. This keeps all arguments
      # to exp() non-positive while preserving the empirical CGF exactly.
      log.weights = outer(resid_nonzero, t.chunk, "*")
      shift = ifelse(t.chunk >= 0,
                     t.chunk * resid.max,
                     t.chunk * resid.min)
      if(n.zero != 0)
        shift = pmax(shift, 0)
      weights = exp(sweep(log.weights, 2, shift, "-"))

      zero.weight = n.zero * exp(-shift)
      weight.sum = colSums(weights) + zero.weight
      weighted.mean = colSums(resid_nonzero * weights)/weight.sum
      weighted.second = colSums(resid.squared * weights)/weight.sum
      weighted.var = weighted.second - weighted.mean^2

      # Recalculate tail columns around their weighted mean when the raw
      # second-moment identity is vulnerable to cancellation.
      variance.scale = pmax(weighted.second,
                            weighted.mean^2,
                            .Machine$double.xmin)
      unstable = which(!is.finite(weighted.var) |
                         weighted.var < sqrt(.Machine$double.eps) * variance.scale)
      for(j in unstable){
        centered = resid_nonzero - weighted.mean[j]
        weighted.var[j] =
          (sum(weights[,j] * centered^2) +
             zero.weight[j] * weighted.mean[j]^2)/weight.sum[j]
      }

      cumul[chunk.start:chunk.end, 1] = shift + log(weight.sum) - log(n.total)
      cumul[chunk.start:chunk.end, 2] = weighted.mean
      cumul[chunk.start:chunk.end, 3] = weighted.var

      while(chunk.end >= next.print){
        print(paste0("Complete ",next.print,"/",length.out,"."))
        next.print = next.print + 1000
      }
    }

    cumul
}

SPACox_rowsum_vector = function(x, group, n.groups)
{
  out = numeric(n.groups)
  x.sum = rowsum(matrix(x, ncol=1), group, reorder=FALSE)
  out[as.integer(rownames(x.sum))] = x.sum[,1]
  out
}

SPACox_rowsum_columns = function(x, group, n.groups)
{
  out = matrix(0, nrow=nrow(x), ncol=n.groups)
  x.sum = rowsum(t(x), group, reorder=FALSE)
  out[, as.integer(rownames(x.sum))] = t(x.sum)
  rownames(out) = rownames(x)
  out
}

#' SaddlePoint Approximation implementation of a surival analysis
#'
#' A fast and accurate method for a genome-wide survival analysis on a large-scale dataset.
#' @param obj.null an R object returned from function SPACox_Null_Model()
#' @param Geno.mtx a numeric genotype matrix with each row as an individual and each column as a genetic variant.
#'                 Column names of genetic variations and row names of subject IDs are required.
#'                 Missng genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF. Any SNP with MAF < cutoff will be excluded from the analysis.
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#'               If the test statistic lies within the standard deviation cutoff, its p value is calculated based on a normal distribution approximation,
#'               otherwise, its p value is calculated based on a saddlepoint approximation.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#'                      "fixed" imputes missing genotypes (NA) by assigning the mean genotype value (i.e. 2p where p is MAF).
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#'                       Any variant with missing rate higher than this cutoff will be excluded from the analysis.
#' @param CovAdj.cutoff a numeric value (default: 5e-5). If the p-value is less than this cutoff, then we would use an additional technic to adjust for covariates.
#' @param G.model a character string (default: "Add") to specify the genetic model. Options are "Add", "Dom", and "Rec".
#' @details To run SPACox, the following two steps are required:
#' \itemize{
#'   \item Step 1. Use function SPACox_Null_Model() to fit a null Cox model.
#'   \item Step 2: Use function SPACox() to calculate p value for each genetic variant.
#' }
#'
#' SPACox uses a hybrid strategy with both saddlepoint approximation and normal distribution approximation.
#' Generally speaking, saddlepoint approximation is more accurate than, but a little slower than, the traditional normal distribution approximation.
#' Hence, when the score statistic is close to 0 (i.e. p-values are not small), we use the normal distribution approximation.
#' And when the score statistic is far away from 0 (i.e. p-values are small), we use the saddlepoint approximation.
#' Argument 'Cutoff' is to specify the standard deviation cutoff.
#'
#' To calibrate the score statistics, SPACox uses martingale residuals which are calculated via R package survival.
#' All extentions (such as strata, ties, left-censoring) supported by package survival could also be used in SPACox.
#' Time-varying covariates are also supported by splitting each subject into several observations.
#' Simulation studies and real data analyses indicate that SPACox works well if one subject corresponds to 2~3 observations.
#' While, if there are more than 4 observations for each subject, SPACox has not been fully evaluated and the results should be carefully intepreted.
#'
#' Sometimes, the order of subjects between phenotype data and genotype data are different, which could lead to some errors.
#' To avoid that, we ask users to specify the IDs of both phenotype data (pIDs) and genotype data (gIDs) when fitting the null model.
#' Users are responsible to check the consistency between pIDs and formula, and the consistency between gIDs and Geno.mtx.
#'
#' @return an R matrix with the following columns
#' \item{MAF}{Minor allele frequencies}
#' \item{missing.rate}{Missing rates}
#' \item{p.value.spa}{p value (recommanded) from a saddlepoint approximation.}
#' \item{p.value.norm}{p value from a normal distribution approximation.}
#' \item{Stat}{score statistics}
#' \item{Var}{estimated variances of the score statistics}
#' \item{z}{z values corresponding to the score statistics}
#' @examples
#' \dontrun{
#' # Simulation phenotype and genotype
#' N = 10000
#' nSNP = 1000
#' MAF = 0.1
#' Phen.mtx = data.frame(ID = paste0("IID-",1:N),
#'                       event=rbinom(N,1,0.5),
#'                       time=runif(N),
#'                       Cov1=rnorm(N),
#'                       Cov2=rbinom(N,1,0.5))
#' Geno.mtx = matrix(rbinom(N*nSNP,2,MAF),N,nSNP)
#'
#' # NOTE: The row and column names of genotype matrix are required.
#' rownames(Geno.mtx) = paste0("IID-",1:N)
#' colnames(Geno.mtx) = paste0("SNP-",1:nSNP)
#' Geno.mtx[1:10,1]=NA   # please use NA for missing genotype
#'
#' # Attach the survival package so that we can use its function Surv()
#' library(survival)
#' obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' SPACox.res = SPACox(obj.null, Geno.mtx)
#'
#' # we recommand using column of 'p.value.spa' to associate genotype with time-to-event phenotypes
#' head(SPACox.res)
#'
#' ## missing data in response/indicator variables is also supported. Please do not remove pIDs of subjects with missing data, the program will do it.
#' Phen.mtx$event[2] = NA
#' Phen.mtx$Cov1[5] = NA
#' obj.null = SPACox_Null_Model(Surv(time,event)~Cov1+Cov2, data=Phen.mtx,
#'                              pIDs=Phen.mtx$ID, gIDs=rownames(Geno.mtx))
#' SPACox.res = SPACox(obj.null, Geno.mtx)
#'
#' # The below is an example code to use survival package
#' coxph(Surv(time,event)~Cov1+Cov2+Geno.mtx[,1], data=Phen.mtx)
#' }
#' @export
SPACox = function(obj.null,
                  Geno.mtx,
                  Cutoff = 2,
                  impute.method = "fixed",
                  missing.cutoff = 0.15,
                  min.maf = 0.0001,
                  CovAdj.cutoff = 5e-5,
                  G.model = "Add")
{
  ## check input
  par.list = list(pwd=getwd(),
                  sessionInfo=sessionInfo(),
                  Cutoff=Cutoff,
                  impute.method=impute.method,
                  missing.cutoff=missing.cutoff,
                  min.maf=min.maf,
                  CovAdj.cutoff=CovAdj.cutoff,
                  G.model=G.model)

  check_input1(obj.null, Geno.mtx, par.list)
  print(paste0("Sample size is ",nrow(Geno.mtx),"."))
  print(paste0("Number of variants is ",ncol(Geno.mtx),"."))

  ### Prepare the main output data frame
  n.Geno = ncol(Geno.mtx)
  output = matrix(NA, n.Geno, 7)
  colnames(output) = c("MAF","missing.rate","p.value.spa","p.value.norm","Stat","Var","z")
  rownames(output) = colnames(Geno.mtx)

  ### Start analysis
  print("Start Analyzing...")
  print(Sys.time())

  # Cycle for genotype matrix
  for(i in 1:n.Geno){

    g = Geno.mtx[,i]
    output.one.SNP = SPACox.one.SNP(g,
                                    obj.null,
                                    Cutoff,
                                    impute.method,
                                    missing.cutoff,
                                    min.maf,
                                    CovAdj.cutoff,
                                    G.model)
    output[i,] = output.one.SNP
  }

  print("Analysis Complete.")
  print(Sys.time())
  return(output)
}

#' SaddlePoint Approximation implementation of Cox regression surival analysis (One-SNP-version)
#'
#' One-SNP-version SPACox function. This function is to facilitate users that prefer reading and analyzing genotype line-by-line.
#' @param g a numeric genotype vector. Missing genotype should be coded as NA. Both hard-called and imputed genotype data are supported.
#' @param obj.null an R object returned from function SPACox_Null_Model()
#' @param Cutoff a numeric value (Default: 2) to specify the standard deviation cutoff to be used.
#' @param impute.method a character string (default: "fixed") to specify the method to impute missing genotypes.
#' @param missing.cutoff a numeric value (default: 0.15) to specify the cutoff of the missing rates.
#' @param min.maf a numeric value (default: 0.0001) to specify the cutoff of the minimal MAF.
#' @param CovAdj.cutoff a numeric value (default: 5e-5). If the p-value is less than this cutoff, then we would use an additional technic to adjust for covariates.
#' @param G.model a character string (default: "Add") to specify the genetic model. Options are "Add", "Dom", and "Rec".
#' @return the same as function SPACox.
#' @export
SPACox.one.SNP = function(g,
                          obj.null,
                          Cutoff = 2,
                          impute.method = "fixed",
                          missing.cutoff = 0.15,
                          min.maf = 0.0001,
                          CovAdj.cutoff = 5e-5,
                          G.model = "Add")
{
  g[g==-9]=NA  # since we add plink input
  ## calculate MAF and update genotype vector
  MAF = mean(g, na.rm=T)/2
  n.subjects = length(g)
  pos.na = which(is.na(g))
  missing.rate = length(pos.na)/n.subjects

  if(missing.rate != 0){
    if(impute.method=="fixed")
      g[pos.na] = 2*MAF
  }

  if(MAF > 0.5){
    MAF = 1-MAF
    g = 2-g
  }

  if(G.model=="Add"){}   # do nothing if G.Model is "Add"
  if(G.model=="Dom") g = ifelse(g>=1,1,0)
  if(G.model=="Rec") g = ifelse(g<=1,0,1)

  if(MAF < min.maf || missing.rate > missing.cutoff)
    return(c(MAF, missing.rate, NA, NA, NA, NA, NA))

  use.grouped = SPACox_can_use_grouped_genotype(g, obj.null)

  if(use.grouped){
    n.rows = length(obj.null$resid)
    if(sum(obj.null$row_count_by_genotype) != n.rows)
      stop("sum(obj.null$row_count_by_genotype) should equal length(obj.null$resid).")

    ## Score statistic
    S = sum(g * obj.null$resid_sum_by_genotype)

    ## estimated variance without adjusting for covariates
    G1 = g - 2*MAF   # centered genotype (such that mean=0)
    S.var1 = obj.null$var.resid * sum(obj.null$row_count_by_genotype * G1^2)
    z1 = S/sqrt(S.var1)

    if(abs(z1) < Cutoff){
      pval.norm = pnorm(abs(z1), lower.tail = FALSE)*2
      return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var1, z1))
    }

    G1norm = G1/sqrt(S.var1)  # normalized genotype (such that sd=1)
    G1.grouped = SPACox_group_values(G1norm, obj.null$row_count_by_genotype)

    pval1 = GetProb_SPA_grouped(obj.null, G1.grouped$values, G1.grouped$counts, abs(z1), lower.tail = FALSE)
    pval2 = GetProb_SPA_grouped(obj.null, G1.grouped$values, G1.grouped$counts, -abs(z1), lower.tail = TRUE)
    pval = pval1 + pval2

    if(pval[1] > CovAdj.cutoff)
      return(c(MAF, missing.rate, pval, S, S.var1, z1))

    ## estimated variance after adjusting for covariates
    g.row = g[obj.null$row_to_genotype]
    tXg = obj.null$tX_by_genotype %*% g
    G2 = as.vector(g.row - obj.null$X.invXX %*% tXg)
    S.var2 = obj.null$var.resid * sum(G2^2)
    z2 = S/sqrt(S.var2)

    G2norm = G2/sqrt(S.var2)

    N1set = seq_len(n.rows)
    N0 = 0
    G2N1 = G2norm
    G2N0 = 0   # since N0=0, this value actually does not matter

    pval1 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, abs(z2), lower.tail = FALSE)
    pval2 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, -abs(z2), lower.tail = TRUE)
    pval = pval1 + pval2

    return(c(MAF, missing.rate, pval, S, S.var2, z2))
  }

  if(!is.null(obj.null$p2g))
    g = g[obj.null$p2g]

  n.rows = length(g)
  SPACox_check_row_genotype_length(n.rows, obj.null)

  ## Score statistic
  S = sum(g * obj.null$resid)

  ## estimated variance without adjusting for covariates
  G1 = g - 2*MAF   # centered genotype (such that mean=0)
  S.var1 = obj.null$var.resid * sum(G1^2)
  z1 = S/sqrt(S.var1)

  if(abs(z1) < Cutoff){
    pval.norm = pnorm(abs(z1), lower.tail = FALSE)*2
    return(c(MAF, missing.rate, pval.norm, pval.norm, S, S.var1, z1))
  }

  N1set = which(g!=0)  # position of non-zero genotypes
  N0 = n.rows-length(N1set)

  G1norm = G1/sqrt(S.var1)  # normalized genotype (such that sd=1)

  G1N1 = G1norm[N1set]
  G1N0 = -2*MAF/sqrt(S.var1)   # all subjects with g=0 share the same normlized genotype, this is to reduce computation time

  pval1 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, abs(z1), lower.tail = FALSE)
  pval2 = GetProb_SPA(obj.null, G1N1, G1N0, N1set, N0, -abs(z1), lower.tail = TRUE)
  pval = pval1 + pval2

  if(pval[1] > CovAdj.cutoff)
    return(c(MAF, missing.rate, pval, S, S.var1, z1))

  ## estimated variance after adjusting for covariates

  G2 = g - obj.null$X.invXX %*% (obj.null$tX[,N1set,drop=F] %*% g[N1set])
  S.var2 = obj.null$var.resid * sum(G2^2)
  z2 = S/sqrt(S.var2)

  G2norm = G2/sqrt(S.var2)

  N1set = seq_len(n.rows)
  N0 = 0
  G2N1 = G2norm
  G2N0 = 0   # since N0=0, this value actually does not matter

  pval1 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, abs(z2), lower.tail = FALSE)
  pval2 = GetProb_SPA(obj.null, G2N1, G2N0, N1set, N0, -abs(z2), lower.tail = TRUE)
  pval = pval1 + pval2

  return(c(MAF, missing.rate, pval, S, S.var2, z2))
}

SPACox_can_use_grouped_genotype = function(g, obj.null)
{
  !is.null(obj.null$row_count_by_genotype) &&
    !is.null(obj.null$resid_sum_by_genotype) &&
    !is.null(obj.null$row_to_genotype) &&
    !is.null(obj.null$tX_by_genotype) &&
    length(g) == length(obj.null$row_count_by_genotype)
}

SPACox_check_row_genotype_length = function(n.rows, obj.null)
{
  if(n.rows != length(obj.null$resid))
    stop("length(g) after matching genotype IDs should equal length(obj.null$resid).")
  if(n.rows != nrow(obj.null$X.invXX))
    stop("length(g) after matching genotype IDs should equal nrow(obj.null$X.invXX).")
  if(n.rows != ncol(obj.null$tX))
    stop("length(g) after matching genotype IDs should equal ncol(obj.null$tX).")
}

SPACox_group_values = function(values, counts)
{
  keep = counts != 0
  values = values[keep]
  counts = counts[keep]

  if(length(values) == 0)
    return(list(values=numeric(0), counts=numeric(0)))

  unique.values = unique(values)
  group = match(values, unique.values)
  grouped.counts = as.numeric(tapply(counts, group, sum))

  list(values=unique.values, counts=grouped.counts)
}

GetProb_SPA = function(obj.null, G2NB, G2NA, NBset, N0, q2, lower.tail){

  out = uniroot(K1_adj, c(-20,20), extendInt = "upX",
                G2NB=G2NB, G2NA=G2NA, NBset=NBset,
                N0=N0, q2=q2, obj.null=obj.null)
  zeta = out$root

  k1 = K_org(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)
  k2 = K2(zeta,  G2NB=G2NB, G2NA=G2NA, NBset=NBset, N0=N0, obj.null=obj.null)

  temp1 = zeta * q2 - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  pval.norm = pnorm(q2, lower.tail = lower.tail)

  re = c(pval, pval.norm)
  return(re)
}

GetProb_SPA_grouped = function(obj.null, G2N, counts, q2, lower.tail){

  out = uniroot(K1_adj_grouped, c(-20,20), extendInt = "upX",
                G2N=G2N, counts=counts,
                q2=q2, obj.null=obj.null)
  zeta = out$root

  k1 = K_org_grouped(zeta, G2N=G2N, counts=counts, obj.null=obj.null)
  k2 = K2_grouped(zeta, G2N=G2N, counts=counts, obj.null=obj.null)

  temp1 = zeta * q2 - k1

  w = sign(zeta) * (2 *temp1)^{1/2}
  v = zeta * (k2)^{1/2}

  pval = pnorm(w + 1/w * log(v/w), lower.tail = lower.tail)
  pval.norm = pnorm(q2, lower.tail = lower.tail)

  re = c(pval, pval.norm)
  return(re)
}


K_org = function(t, G2NB, G2NA, NBset, N0, obj.null){

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*obj.null$K_org_emp(t2NA) + sum(obj.null$K_org_emp(t2NB))
  }
  return(out)
}

K_org_grouped = function(t, G2N, counts, obj.null){

  n.t = length(t)
  out = rep(0,n.t)
  for(i in 1:n.t){
    t1 = t[i]
    t2N = t1*G2N
    out[i] = sum(counts*obj.null$K_org_emp(t2N))
  }
  return(out)
}

K1_adj = function(t, G2NB, G2NA, NBset, N0, q2, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA*obj.null$K_1_emp(t2NA) + sum(G2NB*obj.null$K_1_emp(t2NB)) - q2
  }
  return(out)
}

K1_adj_grouped = function(t, G2N, counts, q2, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2N = t1*G2N
    out[i] = sum(counts*G2N*obj.null$K_1_emp(t2N)) - q2
  }
  return(out)
}

K2 = function(t, G2NB, G2NA, NBset, N0, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2NA = t1*G2NA
    t2NB = t1*G2NB
    out[i] = N0*G2NA^2*obj.null$K_2_emp(t2NA) + sum(G2NB^2*obj.null$K_2_emp(t2NB))
  }
  return(out)
}

K2_grouped = function(t, G2N, counts, obj.null)
{
  n.t = length(t)
  out = rep(0,n.t)

  for(i in 1:n.t){
    t1 = t[i]
    t2N = t1*G2N
    out[i] = sum(counts*G2N^2*obj.null$K_2_emp(t2N))
  }
  return(out)
}


check_input = function(pIDs, gIDs, obj.coxph, range)
{
  if(is.null(pIDs) & is.null(gIDs))
    stop("Arguments 'pIDs' and 'gIDs' are required in case of potential errors. For more information, please refer to 'Details'.")

  pIDs = as.character(pIDs)
  gIDs = as.character(gIDs)
  if(!is.null(obj.coxph$na.action)){
    posNA = c(obj.coxph$na.action)
    if(any(posNA > length(pIDs)))
      stop("Number of input data is larger than length(pIDs).")
    pIDsNA = pIDs[posNA]

    print(paste0("Due to missing data in response/indicators, ",length(posNA)," entries are removed from analysis."))
    print("If concerned about the power loss, we suggest users impute data first and then use SPACox package.")
    print(head(cbind(posNA=posNA, pIDsNA=pIDsNA)))

    pIDs = pIDs[-1*posNA]  # remove IDs with missing data
  }

  if(any(!is.element(pIDs, gIDs)))
    stop("All elements in pIDs should be also in gIDs.")

  if(anyDuplicated(gIDs)!=0)
    stop("Argument 'gIDs' should not have a duplicated element.")

  if(range[2]!=-1*range[1])
    stop("range[2] should be -1*range[1]")

  mresid = obj.coxph$residuals

  if(length(mresid)!=length(pIDs))
    stop("length(mresid)!=length(pIDs) where mresid is the martingale residuals from coxph() in survival package.")

  p2g = NULL
  if(length(pIDs)!=length(gIDs)){
    p2g = match(pIDs, gIDs)
  }else{
    if(any(pIDs != gIDs))
      p2g = match(pIDs, gIDs)
  }

  return(list(p2g=p2g,pIDs=pIDs))
}

check_input1 = function(obj.null, Geno.mtx, par.list)
{
  if(class(obj.null)!="SPACox_NULL_Model")
    stop("obj.null should be a returned outcome from SPACox_Null_Model()")

  if(any(obj.null$gIDs != rownames(Geno.mtx))) stop("gIDs should be the same as rownames(Geno.mtx).")
  if(is.null(rownames(Geno.mtx))) stop("Row names of 'Geno.mtx' should be given.")
  if(is.null(colnames(Geno.mtx))) stop("Column names of 'Geno.mtx' should be given.")
  if(!is.numeric(Geno.mtx)|!is.matrix(Geno.mtx)) stop("Input 'Geno.mtx' should be a numeric matrix.")

  if(!is.numeric(par.list$min.maf)|par.list$min.maf<0|par.list$min.maf>0.5) stop("Argument 'min.maf' should be a numeric value >= 0 and <= 0.5.")
  if(!is.numeric(par.list$Cutoff)|par.list$Cutoff<0) stop("Argument 'Cutoff' should be a numeric value >= 0.")
  # if(!is.element(par.list$impute.method,c("none","bestguess","random","fixed"))) stop("Argument 'impute.method' should be 'none', 'bestguess', 'random' or 'fixed'.")
  if(!is.element(par.list$impute.method,c("fixed"))) stop("Argument 'impute.method' should be 'fixed'.")
  if(!is.numeric(par.list$missing.cutoff)|par.list$missing.cutoff<0|par.list$missing.cutoff>1) stop("Argument 'missing.cutoff' should be a numeric value between 0 and 1.")
  if(!is.element(par.list$G.model,c("Add","Dom","Rec"))) stop("Argument 'G.model' should be 'Add', 'Dom' or 'Rec'.")
}

