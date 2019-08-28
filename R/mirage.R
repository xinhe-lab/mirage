#' Mixture model based rare variant analysis
#'
#' This function implements the gene based version of rare variant test
#' with MIRAGE model
#' @param data variant count summary statistics, a 4 column data frame for 1) loci ID, 2) count in group 1, 
#' 3) count in group 2 and 4) functional category for a variant.
#' @param n1 sample size for group 1; either an integer number if all variants have the same sample size, or a vector
#' of length of number of variants to specify sample size for each variant.
#' @param n2 sample size for group 2; either an integer number if all variants have the same sample size, or a vector
#' of length of number of variants to specify sample size for each variant.
#' @param gamma a vector of category specific effect size prior
#' @param sigma prior on proportion of risk genes
#' @param delta ... FIXME
#' @param eta.init ... FIXME (do you have default eta.init? zero?)
#' @return \item{BF}{Bayes factor of genes}
#' \item{delta.est}{...}
#' \item{delta.pvalue}{...}
#' \item{eta.est}{...}
#' \item{pvalue}{...}
#' @examples
#' mirage(...)
#' @importFrom progress progress_bar
#' @export
mirage = function(data, N1, N0, gamma, sigma, eta.init, delta = NULL, max.iter = 10000, 
    tol = 1e-05, verbose = TRUE) {
    # Input check & initialize
    if (!is.null(eta.init)) {
        eta.k[1, ] = eta.init
     } else {
         stop("FIXME: need default input")
     }
    names(data) = c("name", "C1", "C2", "category")
    unique.gene = unique(data$name)
    num.gene = length(unique.gene)
    num.group = length(unique(data$category))
    eta.k = matrix(nrow = max.iter, ncol = num.group)
    BF.gene = matrix(1, nrow = max.iter, ncol = num.gene)
    LoF.BF.gene = matrix(1, nrow = max.iter, ncol = num.gene)
    nonLoF.BF.gene = matrix(1, nrow = max.iter, ncol = num.gene)
    BF.genevar = list()
    delta.est = numeric()
    if (!is.null(delta)) { 
        delta.est[1] = delta
    } else {
        stop("FIXME: need default input")
    }
    # progress bar
    if (verbose) {
        pb = progress_bar$new(format = "[:spin] Initial analysis of unit :unit out of :total :elapsed", 
            clear = FALSE, total = num.gene, show_after = 0.5) 
    } else { 
        pb = null_progress_bar$new()
    }
    # calculate the Bayes factor for variant j in gene i
    for (i in 1:num.gene) {
        var.index.list = which(data$name == unique.gene[i])
        var.BF = numeric()
        if (length(var.index.list) > 0) {
            # calculate Bayes factor for variant (i,j)
            for (j in 1:length(var.index.list)) {
                category = data$category[var.index.list[j]]
                # FIXME: explain what 5 does here
                if (category <= 5) {
                  # FIXME: have to handle N1 N2 being vectors
                  var.BF[j] = BF.var.inte(data$C1[var.index.list[j]], data$C2[var.index.list[j]], 
                    bar.gamma = 6, sig = sigma, N1, N2)
                } else {
                  var.BF[j] = BF.var.inte(data$C1[var.index.list[j]], data$C2[var.index.list[j]], 
                    bar.gamma = gamma.mean, sig = sigma, N1, N2)
                }
                mixed_bf = 1 - eta.k[1, category] + eta.k[1, category] * var.BF[j]
                BF.gene[1, i] = BF.gene[1, i] * mixed_bf
                ################## split BF of LoF and non LoF FIXME: explain what 2 does here
                if (category <= 2) }{ 
                  LoF.BF.gene[1, i] = LoF.BF.gene[1, i] * mixed_bf
                 } else {
                     nonLoF.BF.gene[1, i] = nonLoF.BF.gene[1, i] * mixed_bf
                 }
            }
        }
        BF.genevar[[i]] = cbind(var.index.list, var.BF)
        colnames(BF.genevar[[i]]) = c("idx", "BF")
        pb$tick(tokens = list(total = num.gene, unit = i))
    }
    # EM algorithm
    if (verbose) 
        pb = progress_bar$new(format = "[:spin] Iteration :iteration (diff = :delta) :elapsed", 
            clear = TRUE, total = max.iter, show_after = 0.5)
    for (iter in 2:max.iter) {
        prev_iter = iter - 1
        # E step expectation for variant (i,j), every gene may have varying number of
        # variant
        EUiZij = list()
        # expectation for gene i.
        EUi = numeric()
        total.Zij = matrix(nrow = num.gene, ncol = num.group)
        total.Ui = matrix(nrow = num.gene, ncol = num.group)  # used to estimate eta
        for (i in 1:num.gene) {
            UiZij = numeric()
            if (nrow(BF.genevar[[i]]) > 0) {
                for (j in 1:nrow(BF.genevar[[i]])) {
                  category = data$category[BF.genevar$idx[j]]
                  numer = BF.genevar[[i]]$BF[j] * eta.k[prev_iter, category] * delta.est[prev_iter]
                  denom = (delta.est[prev_iter] + (1 - delta.est[prev_iter])/BF.gene[prev_iter, 
                    i]) * (eta.k[prev_iter, category] * BF.genevar[[i]]$BF[j] + (1 - 
                    eta.k[prev_iter, category]))
                  UiZij[j] = numer/denom
                  mixed_bf = 1 - eta.k[prev_iter, category] + eta.k[prev_iter, 
                    category] * BF.genevar[[i]]$BF[j]
                  BF.gene[iter, i] = BF.gene[iter, i] * mixed_bf
                  ########################### split into LoF and non-LoF two parts
                  if (category <= 2) {
                    LoF.BF.gene[iter, i] = LoF.BF.gene[iter, i] * mixed_bf 
                  } else {
                      nonLoF.BF.gene[iter, i] = nonLoF.BF.gene[iter, i] * mixed_bf
                  }
                }
            }
            EUiZij[[i]] = UiZij
            EUi[i] = delta.est[prev_iter] * bb/(delta.est[prev_iter] * bb + 1 - delta.est[prev_iter])
            ###################### Note here each gene may have multiple annotation groups
            for (g in 1:num.group) {
                total.Zij[i, g] = sum(UiZij[which(data$category[BF.genevar$idx] == 
                  g)], na.rm = TRUE)
                total.Ui[i, g] = sum(sum(UiZij[which(data$category[BF.genevar$idx] == 
                  g)] > 0, na.rm = TRUE) * EUi[i])
            }
        }
        ############## EM algorithm: M step
        delta.est[iter] = sum(EUi)/num.gene
        # delta.est[iter]=delta
        for (g in 1:num.group) {
            if (sum(total.Ui[, g]) != 0) 
                eta.k[iter, g] = sum(total.Zij[, g])/sum(total.Ui[, g])
            if (sum(total.Ui[, g]) == 0) 
                eta.k[iter, g] = 0
        }
        ################ 
        if (num.group > 1) 
            diff = sum(abs(eta.k[iter, ] - eta.k[prev_iter, ]))
        if (diff < tol || iter >= max.iter) {
            pb$tick(max.iter)
            max.iter = iter
            break
        }
        pb$tick(tokens = list(delta = sprintf(diff, fmt = "%#.1e"), iteration = i))
    }  # end of iter
    # FIXME: need to unify this.
    if (num.group > 1) {
        eta.k = eta.k[1:max.iter, ]
    } else {
        eta.k = eta.k[1:max.iter]
    }
    # calculate the LRT statistics and p-value for genes
    if (verbose) 
        cat("Computing gene level LRT statistics and p-values ...\n")
    lkhd = rep(1, num.gene)
    total.lkhd = 0
    teststat = numeric()
    pvalue = numeric()
    for (i in 1:num.gene) {
        if (nrow(BF.genevar[[i]]) > 0) {
            for (j in 1:nrow(BF.genevar[[i]])) {
                category = data$category[BF.genevar$idx[j]]
                lkhd[i] = lkhd[i] * ((1 - eta.k[prev_iter, category]) + eta.k[prev_iter, 
                  category] * BF.genevar$BF[j])
            }
        }
        teststat[i] = 2 * log((1 - delta.est[prev_iter]) + delta.est[prev_iter] * 
            lkhd[i])
        # this is the test statistics of one gene
        total.lkhd = total.lkhd + log((1 - delta.est[prev_iter]) + delta.est[prev_iter] * 
            lkhd[i])
        pvalue[i] = pchisq(teststat[i], 2, lower.tail = F)
    }
    teststat[num.gene + 1] = 2 * total.lkhd
    pvalue[num.gene + 1] = pchisq(teststat[num.gene + 1], 2, lower.tail = F)
    # calculate the LRT statistics and p-value for categories
    if (verbose) 
        cat("Computing LRT statistics and p-values by functional groups ...\n")
    cate.lkhd = rep(1, num.group)
    sum.lkhd = 0
    cate.stat = numeric()
    cate.pvalue = numeric()
    for (g in 1:num.group) {
        total.lkhd = 0
        lkhd.gene = rep(1, num.gene)
        for (i in 1:num.gene) {
            if (nrow(BF.genevar[[i]]) > 0) 
                for (j in 1:nrow(BF.genevar[[i]])) {
                  if (data$category[BF.genevar$idx[j]] == g) {
                    mixed_bf = 1 - eta.k[prev_iter, g] + eta.k[prev_iter, g] * 
                      data$var.BF[j]
                    lkhd.gene[i] = lkhd.gene[i] * mixed_bf
                    cate.lkhd[g] = cate.lkhd[g] * mixed_bf
                  }
                }
            total.lkhd = total.lkhd + log((1 - delta.est[prev_iter]) + delta.est[prev_iter] * 
                lkhd.gene[i])
        }  # end of i
        cate.stat[g] = 2 * total.lkhd
        cate.pvalue[g] = pchisq(2 * total.lkhd, 1, lower.tail = F)
    }
    sum.lkhd = sum(cate.stat)
    ############################################## 
    return(result = list(delta.est = delta.est[prev_iter], delta.pvalue = pvalue[length(pvalue)], 
        eta.est = eta.k[prev_iter, ], eta.stat = cate.stat, eta.pvalue = cate.pvalue, 
        BF.gene = data.frame(Gene = unique.gene, BF = BF.gene[prev_iter, ], LoF.BF = LoF.BF.gene[prev_iter, 
            ], nonLoF.BF = nonLoF.BF.gene[prev_iter, ]), BF.all = BF.genevar, Eui = EUi))
}


intergrand = function(aa, var.case, var.contr, bar.gamma, sig, N1, N0) {
    ff = dbinom(var.case, sum(var.case, var.contr), aa * N1/(aa * N1 + N0)) * dgamma(aa, 
        bar.gamma * sig, sig)
    return(ff)
}
# calculate the bayes factor of a single variant via integration
BF.var.inte = function(var.case, var.contr, bar.gamma, sig, N1, N0) {
    # Under H0: gamma=1
    marglik0.CC <- dbinom(var.case, sum(var.case, var.contr), N1/(N1 + N0))
    # Under H1: gamma~gamma(gamma.mean*sigma, sigma)
    marglik1.CC <- integrate(intergrand, var.case, var.contr, bar.gamma, sig, N1, 
        N0, low = 0, upper = 100, stop.on.error = F)$value
    BF.var <- marglik1.CC/marglik0.CC
    return(BF.var)
}