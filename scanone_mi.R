################################################################################
# Utility functions for DO Genotype HMM.
# Daniel Gatti
# Dan.Gatti@jax.org
# Jan. 12, 2012
################################################################################
convert.pos.to.bp = function(pos) {
  if(max(pos) <= 200) {
    pos = pos * 1e6
  } # if(max(pos] <= 200))
  return(pos)
}
################################################################################
# Get the number of autosomes.  We assume that any chromosome with a number is
# an autosome.
# Arguments: snps: matrix of SNPs with at least 4 columns. SNP ID, chr, Mb and
#                  cM in columns 1:4.
# Returns: integer: the number of autosomes.
get.num.auto = function(snps) {
  # Turn off the warnings for a moment while we get the number of autosomes.
  old.warn = options("warn")$warn
  options(warn = -1)
  num.auto = max(as.numeric(snps[,2]), na.rm = TRUE)
  options(warn = old.warn)
  return(num.auto)
  
} # get.num.auto()
################################################################################
# Get the DO genotype states.
# Returns: character vector with the 36 sorted DO states.
get.do.states = function() {
  states = outer(LETTERS[1:8], LETTERS[1:8], paste, sep = "")
  states = sort(states[upper.tri(states, diag = TRUE)])
  states
} # get.do.states()
qtl.qtlrel.mi = function(pheno, probs, K, addcovar, intcovar, snps) {
ranMtr<- rem(~cage, data=g15.fam)
  pheno = as.matrix(pheno)
  if(!missing(K)) {
    K = as.matrix(K)
  } # if(!missing(K))
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) > 200)
  # Return value.
  retval = NULL
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3],
               snp = snps[,1])
  vTmp = list(AA = NULL, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
              EE = diag(length(pheno)),cage=ranMtr$cage)
  if(!missing(K)) {
    vTmp$AA = 2 * K
  } # if(!missing(K))
  # This tells QTLRel to fit the additive model.
  class(prdat) = c(class(prdat), "addEff") 
  vc = NULL
  if(missing(addcovar)) {
    # No covariates.
    vc = estVC(y = pheno, v = vTmp)
    res = scanOne(y = pheno, prdat = prdat, vc = vc, test = "None", numGeno = TRUE)
  } else {
    vc = estVC(y = pheno, x = addcovar, v = vTmp)
    if(missing(intcovar)) {
      # Additive covariates only.
      res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
            numGeno = TRUE, test = "None")
    } else {
      # Additive and interactive covariates.
      res = scanOne(y = pheno, x = addcovar, prdat = prdat, vc = vc, 
            intcovar = intcovar, numGeno = TRUE, test = "None")
    } # else
  } # else
  # Convert the model coefficients to a matrix.
  coef = matrix(unlist(res$parameters), length(res$parameters),
         length(res$parameters[[1]]), dimnames = list(res$snp,
         names(res$parameters[[1]])), byrow = TRUE)
  # Return the LRS, LOD, p-value and -log10(p-value).
  p = pchisq(q = res$p, df = dim(probs)[[2]] - 1, lower.tail = FALSE)
  return(list(lod = cbind(snps[,1:4], perc.var = res$v, lrs = res$p,
           lod = res$p / (2 * log(10)), p = p,
           neg.log10.p = -log(p, 10)), coef = coef))
} # qtl.qtlrel()



fast.qtlrel.mi = function(pheno, probs, K, addcovar, snps) {
ranMtr<- rem(~cage, data=g15.fam)
  # If the SNPs are in bp, rescale them to Mb.
  if(max(snps[,3], na.rm = TRUE) > 200) {
    snps[,3] = snps[,3] * 1e-6
  } # if(max(snps[,3]) < 200)
  # Return value.
  retval = NULL
  prdat = list(pr = probs, chr = snps[,2], dist = snps[,3], snp = snps[,1])
  class(prdat) = c(class(prdat), "addEff")
  err.cov = diag(length(pheno))
  if(!missing(K)) {
    K = as.matrix(K)
    vTmp = list(AA = 2 * K, DD = NULL, HH = NULL, AD = NULL, MH = NULL,
                EE = err.cov, cage=ranMtr$cage)
    vc = NULL
    if(missing(addcovar)) {
      vc = estVC(y = pheno, v = vTmp)
    } else {
      vc = estVC(y = pheno, x = addcovar, v = vTmp)
    } # else
    err.cov = matrix(0, nrow(K), ncol(K))
    num.cov = 1
    if(!missing(addcovar)) {
      num.cov = ncol(addcovar) + 1
    } # if(missing(addcovar))
    for(j in which(vc$nnl)) {
      err.cov = err.cov + vTmp[[j]] * vc$par[names(vTmp)[j]]
    } # for(j)
    rm(vTmp)
  } # if(!missing(K))
  # Invert the covariance matrix.
  eW = eigen(err.cov, symmetric = TRUE)
  if (min(eW$values) < 0 && abs(min(eW$values)) > sqrt(.Machine$double.eps)) {
    stop("'W' is not positive definite")
  } else {
    eW$values[eW$values <= 0] = Inf
  } # else
  err.cov = diag(eW$values^-0.5) %*% t(eW$vector)
  rm(eW)
  # Add the intercept.
  if(missing(addcovar)) {
    addcovar = matrix(1, nrow = length(pheno), ncol = 1)
  } else {
    addcovar = as.matrix(cbind(rep(1, length(pheno)), addcovar))
  } # else
  colnames(addcovar)[1] = "Intercept"
  # Null model.
  ytmp = err.cov %*% pheno
  xtmp = err.cov %*% addcovar
  qr.null = qr(xtmp)
  ss.null = sum(qr.resid(qr.null, ytmp)^2)
  # Additive model for all SNPs.
  addx = cbind(addcovar, probs[,-1,1])
  ss = rep(0, nrow(snps))
  coef = matrix(0, nrow(snps), ncol(addx), dimnames = list(snps[,1],
         colnames(addx)))
  rng = (ncol(addcovar)+1):ncol(addx)
  for(s in 1:nrow(snps)) {
    addx[,rng] = probs[,-1,s]
    xtmp = err.cov %*% addx
    qr.add = qr(xtmp)
    ss[s] = sum(qr.resid(qr.add, ytmp)^2)
    coef[s,] = qr.coef(qr.add, ytmp)
  } # for(s)
  perc.var = 100 * (1.0 - (ss / ss.null))
  lrs = -length(pheno) * log(ss / ss.null)
  lod = lrs / (2 * log(10))
  # Get the p-value from the LRS.
  p = pchisq(q = -length(pheno) * log(ss / ss.null), 
      df = ncol(addx) - ncol(addcovar), lower.tail = FALSE)
  return(list(lod = cbind(snps[,1:4], perc.var  = perc.var, lrs = lrs, 
         lod = lod, p = p, neg.log10.p = -log(p, 10)), coef = coef))
} # fast.qtlrel
scanone_mi = function (pheno, pheno.col = 1, probs, K, addcovar, intcovar, 
    snps, model = c("additive", "dominance", "full")) 
{
    model = match.arg(model)
    if (is.null(rownames(pheno))) {
        stop("rownames(pheno) is null. The sample IDs must be in rownames(pheno).")
    }
    num.auto = get.num.auto(snps)
    if (is.character(pheno.col)) {
        pheno.col = match(pheno.col, colnames(pheno))
    }
    pheno = pheno[rownames(pheno) %in% dimnames(probs)[[1]], 
        , drop = F]
    probs = probs[dimnames(probs)[[1]] %in% rownames(pheno), 
        , ]
    probs = probs[match(rownames(pheno), dimnames(probs)[[1]]), 
        , ]
    if (any(dim(probs) == 0)) {
        stop(paste("There are no matching samples in pheno and probs. Please", 
            "verify that the sample IDs in rownames(pheno) match the sample", 
            "IDs in dimnames(probs)[[1]]."))
    }
    snps = snps[snps[, 1] %in% dimnames(probs)[[3]], ]
    probs = probs[, , match(snps[, 1], dimnames(probs)[[3]])]
    if (!missing(K)) {
        K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% 
            rownames(pheno)]
        K = K[match(rownames(pheno), rownames(K)), match(rownames(pheno), 
            colnames(K))]
    }
    xchr = which(snps[, 2] %in% "X")
    sex = NULL
    if (length(xchr) > 0) {
        sex.col = grep("sex", colnames(pheno), ignore.case = T)
        if (length(sex.col) == 0) {
            if (!missing(addcovar)) {
                sex.col = grep("sex", colnames(addcovar), ignore.case = T)
            }
            if (length(sex.col) == 0) {
                stop(paste("There is no sex column in the phenotypes or covariates.", 
                  "Please add a sex column for proper mapping on the X chromosome."))
            }
            else {
                sex = addcovar[, sex.col]
            }
        }
        else {
            sex = pheno[, sex.col]
        }
        sex = toupper(sex)
    }
    if (!missing(addcovar)) {
        if (is.null(rownames(addcovar))) {
            stop("rownames(addcovar) is null. The sample IDs must be in rownames(addcovar).")
        }
        addcovar = data.frame(addcovar, stringsAsFactors = T)
        rn = rownames(addcovar)
        addcovar = lapply(addcovar, as.numeric)
        cn = names(addcovar)
        addcovar = matrix(unlist(addcovar), length(addcovar[[1]]), 
            length(addcovar), dimnames = list(rn, cn))
        rownames(addcovar) = rn
        addcovar = addcovar[rownames(addcovar) %in% rownames(pheno), 
            , drop = F]
        addcovar = addcovar[match(rownames(pheno), rownames(addcovar)), 
            , drop = F]
        if (sum(rownames(pheno) %in% rownames(addcovar)) == 0) {
            stop(paste("rownames(pheno) does not contain any sample IDs in", 
                "common with rownames(addcovar). Please make sure that the", 
                "rownames in pheno and addcovar match."))
        }
        pheno = pheno[rownames(pheno) %in% rownames(addcovar), 
            , drop = F]
        addcovar = addcovar[rownames(addcovar) %in% rownames(pheno), 
            , drop = F]
        addcovar = addcovar[match(rownames(pheno), rownames(addcovar)), 
            , drop = F]
        probs = probs[, , dimnames(probs)[[3]] %in% snps[, 1]]
        if (!missing(K)) {
            K = K[rownames(K) %in% rownames(pheno), colnames(K) %in% 
                rownames(pheno)]
        }
        if (is.null(colnames(addcovar))) {
            colnames(addcovar) = paste("addcovar", 1:ncol(addcovar), 
                sep = ".")
        }
        addcovar = as.matrix(addcovar)
    }
    if (!missing(intcovar)) {
        if (is.null(rownames(intcovar))) {
            stop("rownames(intcovar) is null. The sample IDs must be in rownames(intcovar).")
        }
        intcovar = as.matrix(intcovar)
        intcovar = intcovar[rownames(intcovar) %in% rownames(pheno), 
            , drop = F]
        intcovar = intcovar[match(rownames(pheno), rownames(intcovar)), 
            , drop = F]
        if (is.null(colnames(intcovar))) {
            colnames(intcovar) = paste("intcovar", 1:ncol(intcovar), 
                sep = ".")
        }
    }
    retval = as.list(1:length(pheno.col))
    names(retval) = colnames(pheno)[pheno.col]
    index = 1
    for (i in pheno.col) {
        print(colnames(pheno)[i])
        p = pheno[, i]
        names(p) = rownames(pheno)
        keep = which(!is.na(p) & !is.nan(p) & !is.infinite(p))
        auto = which(snps[, 2] %in% 1:num.auto)
        if (missing(addcovar)) {
        	print("using normal fast.qtlrel");
            auto.qtl = fast.qtlrel.mi(pheno = p[keep], probs = probs[keep, 
                , auto], K = K[keep, keep], snps = snps[auto, 
                ])
        }
        else {
            keep = intersect(keep, which(rowSums(is.na(addcovar)) == 
                0 & rowSums(is.nan(addcovar)) == 0 & rowSums(is.infinite(addcovar)) == 
                0))
            if (missing(intcovar)) {
            print("using normal fast.qtlrel");
                auto.qtl = fast.qtlrel.mi(pheno = p[keep], probs = probs[keep, 
                  , auto], K = K[keep, keep], addcovar = addcovar[keep, 
                  , drop = F], snps = snps[auto, ])
            }
            else {
            	print("using normal qtl.qtlrel");
                auto.qtl = qtl.qtlrel.mi(pheno = p[keep], probs = probs[keep, 
                  , auto], K = K[keep, keep], addcovar = addcovar[keep, 
                  , drop = F], intcovar = intcovar[keep, , drop = F], 
                  snps = snps[auto, ])
            }
        }
        auto.qtl = list(lod = list(A = auto.qtl$lod), coef = list(A = auto.qtl$coef))
        xchr = which(snps[, 2] %in% "X")
        if (length(xchr) > 0) {
            females = which(sex == "F" | sex == "0")
            males = which(sex == "M" | sex == "1")
            mfprobs = NULL
            if (length(females) > 0 & length(males) > 0) {
                if (model == "additive") {
                  mfprobs = array(0, c(dim(probs)[1], 2 * dim(probs)[2], 
                    length(xchr)), dimnames = list(dimnames(probs)[[1]], 
                    paste(rep(c("F", "M"), each = dim(probs)[2]), 
                      dimnames(probs)[[2]], sep = "."), dimnames(probs)[[3]][xchr]))
                  for (j in females) {
                    mfprobs[j, 1:dim(probs)[2], ] = probs[j, 
                      , xchr]
                  }
                  for (j in males) {
                    mfprobs[j, (dim(probs)[2] + 1):dim(mfprobs)[2], 
                      ] = probs[j, , xchr]
                  }
                }
                else if (model == "full") {
                  tmp = matrix(unlist(strsplit(dimnames(probs)[[2]], 
                    split = "")), nrow = 2)
                  homo = dimnames(probs)[[2]][which(tmp[1, ] == 
                    tmp[2, ])]
                  mfprobs = array(0, c(dim(probs)[1], dim(probs)[2] + 
                    length(homo), length(xchr)), dimnames = list(dimnames(probs)[[1]], 
                    c(paste(rep("F", each = dim(probs)[2]), dimnames(probs)[[2]], 
                      sep = "."), paste("M", homo, sep = ".")), 
                    dimnames(probs)[[3]][xchr]))
                  for (j in females) {
                    mfprobs[j, 1:dim(probs)[2], ] = probs[j, 
                      , xchr]
                  }
                  for (j in males) {
                    mfprobs[j, (dim(probs)[2] + 1):dim(mfprobs)[2], 
                      ] = probs[j, homo, xchr]
                  }
                }
                mfprobs = mfprobs[, -grep("M.A", dimnames(mfprobs)[[2]]), 
                  ]
            }
            else if (length(females) > 0) {
                mfprobs = probs[, , xchr]
                dimnames(mfprobs)[[2]] = paste("F", dimnames(mfprobs)[[2]], 
                  sep = ".")
            }
            else if (length(males) > 0) {
                mfprobs = probs[, , xchr]
                dimnames(mfprobs)[[2]] = paste("M", dimnames(mfprobs)[[2]], 
                  sep = ".")
            }
            if (missing(addcovar)) {
                if (length(grep("sex", colnames(pheno), ignore.case = T)) > 
                  0) {
                  sex = pheno[, grep("sex", colnames(pheno), 
                    ignore.case = T)]
                  sex = as.matrix(as.numeric(sex == "M"))
                  dimnames(sex) = list(rownames(pheno), "sex")
                }
                x.qtl = qtl.qtlrel.mi(pheno = p[keep, drop = F], 
                  probs = mfprobs[keep, , ], K = K[keep, keep], 
                  addcovar = sex[keep, drop = F], snps = snps[xchr, 
                    , drop = F])
                    print("using normal qtl.qtlrel");
            }
            else {
                keep = intersect(keep, which(rowSums(is.na(addcovar)) == 
                  0))
                if (missing(intcovar)) {
                print("using normal qtl.qtlrel");
                  x.qtl = fast.qtlrel.mi(pheno = p[keep], probs = mfprobs[keep, 
                    , ], K = K[keep, keep], addcovar = addcovar[keep, 
                    , drop = F], snps = snps[xchr, ])
                }
                else {
                  x.qtl = qtl.qtlrel.mi(pheno = p[keep], probs = mfprobs[keep, 
                    , ], K = K[keep, keep], addcovar = addcovar[keep, 
                    , drop = F], intcovar = intcovar[keep, , 
                    drop = F], snps = snps[xchr, ])
                    print("using normal qtl.qtlrel");
                }
            }
            auto.qtl$lod = list(A = auto.qtl$lod$A, X = x.qtl$lod)
            auto.qtl$coef = list(A = auto.qtl$coef$A, X = x.qtl$coef)
        }
        retval[[index]] = auto.qtl
        class(retval[[index]]) = c("doqtl", class(retval[[index]]))
        attr(retval[[index]], "model") = "additive"
        index = index + 1
    }
    if (length(retval) == 1) {
        retval = retval[[1]]
    }
    return(retval)
}

