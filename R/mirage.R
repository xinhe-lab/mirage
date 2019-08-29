#' Mixture model based rare variant analysis
#'
#' This function implements rare variant test with MIRAGE model
#'
#' @param data variant count data, a 4 column data frame for 1) locus ID, 2) count in group 1, 
#' 3) count in group 2 and 4) category for a variant. The 4th column is optional.
#' @param n1 sample size for group 1; either an integer number if all variants have the same sample size, or a vector
#' of length of number of variants to specify sample size for each variant.
#' @param n2 sample size for group 2; either an integer number if all variants have the same sample size, or a vector
#' of length of number of variants to specify sample size for each variant.
#' @param gamma a list of category specific shape parameter for effect size prior, or a numeric value if all category share the same effect size.
#' @param sigma a list of category specific scale parameter for effect size prior, or a numeric value if all category share the same effect size.
#' @param delta.init initial value for prior on proportion of risk genes.
#' @param eta.init initial value for prior on proportion of risk variants in a risk gene.
#' @return \item{BF}{Bayes factor of genes}
#' \item{delta.est}{Estimate for proportion of risk genes}
#' \item{delta.pvalue}{Significant test for delta = 0}
#' \item{eta.est}{Estimate for proportion of risk variants in a risk gene}
#' \item{eta.pvalue}{Significant test for eta = 0}
#' \item{BF.gene}{Bayes factor for each gene}
#' \item{BF.all}{Bayes factor for each variant}
#' @examples
#' # see example at https://xinhe-lab.github.io/mirage/articles/mwe.html
#' @importFrom progress progress_bar
#' @export
mirage = function(data, n1, n2, gamma=3, sigma=1, eta.init=0, delta.init=0, estimate.delta = TRUE, max.iter = 10000, 
    tol = 1e-05, verbose = TRUE) {
    # Input check & initialize
    if (ncol(data) == 3) data = c(data, 0)
    if (ncol(data) != 4) stop("Input data should have 3 or 4 columns!")
    names(data) = c("name", "C1", "C2", "category")
    unique.gene = unique(data$name)
    num.gene = length(unique.gene)
    groups = unique(data$category)
    num.group = length(groups)
    data[,4] = match(data[,4], groups)
    eta.k = matrix(nrow = max.iter, ncol = num.group)
    colnames(eta.k) = groups
    eta.k[1, ] = rep(eta.init, num.group)
    delta.est = delta.init
    BF.gene = matrix(1, nrow = max.iter, ncol = num.gene)
    colnames(BF.gene) = unique.gene
    BF.genevar = list()
    # progress bar
    if (verbose && num.gene > 1) {
        pb = progress_bar$new(format = "[:spin] Initial analysis of unit :unit out of :total :elapsed", 
            clear = FALSE, total = num.gene, show_after = 0.5) 
    } else { 
        pb = null_progress_bar$new()
    }
    # calculate the Bayes factor for variant j in gene i
    for (i in 1:num.gene) {
        var.names = rownames(data)[which(data$name == unique.gene[i])]
        var.BF = numeric()
        if (length(var.names) > 0) {
            # calculate Bayes factor for variant (i,j)
            for (j in 1:length(var.names)) {
                category = data[var.names[j],]$category
                # FIXME: have to handle sample size being vectors
                var.BF[j] = BF.var.inte(data[var.names[j],]$C1, data[var.names[j],]$C2,
                                    ifelse(class(gamma) == "list", gamma[category], gamma), 
                                    ifelse(class(sigma) == "list", sigma[category], sigma), 
                                    n1, n2)
                mixed_bf = 1 - eta.k[1, category] + eta.k[1, category] * var.BF[j]
                BF.gene[1, i] = BF.gene[1, i] * mixed_bf
                # FIXME: need to add a utility function that extract BF from given annotation category
            }
        }
        BF.genevar[[unique.gene[i]]] = data.frame(variant=var.names, BF=var.BF)
        pb$tick(tokens = list(total = num.gene, unit = i))
    }
    # EM algorithm
    if (verbose) {
        pb = progress_bar$new(format = "[:spin] Iteration :iteration (diff = :delta) :elapsed", 
            clear = TRUE, total = max.iter, show_after = 0.5)
    } else { 
        pb = null_progress_bar$new()
    }
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
            BF.var = BF.genevar[[unique.gene[i]]]
            if (nrow(BF.var) > 0) {
                for (j in 1:nrow(BF.var)) {
                  category = data[BF.var$variant[j],]$category
                  numer = BF.var$BF[j] * eta.k[prev_iter, category] * delta.est[prev_iter]
                  denom = (delta.est[prev_iter] + (1 - delta.est[prev_iter])/BF.gene[prev_iter, 
                    i]) * (eta.k[prev_iter, category] * BF.var$BF[j] + (1 - 
                    eta.k[prev_iter, category]))
                  UiZij[j] = numer/denom
                  mixed_bf = 1 - eta.k[prev_iter, category] + eta.k[prev_iter, 
                    category] * BF.var$BF[j]
                  BF.gene[iter, i] = BF.gene[iter, i] * mixed_bf
                }
            }
            EUiZij[[i]] = UiZij
            EUi[i] = delta.est[prev_iter] * BF.gene[iter, i]/(delta.est[prev_iter] * BF.gene[iter, i] + 1 - delta.est[prev_iter])
            ###################### Each gene may contain multiple annotation groups
            for (g in 1:num.group) {
                total.Zij[i, g] = sum(UiZij[which(data[BF.var$variant,]$category == 
                  g)], na.rm = TRUE)
                total.Ui[i, g] = sum(sum(UiZij[which(data[BF.var$variant,]$category == 
                  g)] > 0, na.rm = TRUE) * EUi[i])
            }
        }
        ############## EM algorithm: M step
        delta.est[iter] = ifelse(estimate.delta, sum(EUi)/num.gene, delta)
        for (g in 1:num.group) {
            if (sum(total.Ui[, g]) != 0) 
                eta.k[iter, g] = sum(total.Zij[, g])/sum(total.Ui[, g])
            if (sum(total.Ui[, g]) == 0) 
                eta.k[iter, g] = 0
        }
        ################ 
        diff = sum(abs(eta.k[iter, ] - eta.k[prev_iter, ]))
        if (diff < tol || iter >= max.iter) {
            pb$tick(max.iter)
            max.iter = iter
            break
        }
        pb$tick(tokens = list(delta = sprintf(diff, fmt = "%#.1e"), iteration = i))
    }  # end of iter
    eta.k = eta.k[1:max.iter, , drop=FALSE]
    # calculate the LRT statistics and p-value for genes
    if (verbose) 
        cat("Computing gene level LRT statistics and p-values ...\n")
    lkhd = rep(1, num.gene)
    total.lkhd = 0
    teststat = numeric()
    pvalue = numeric()
    for (i in 1:num.gene) {
        BF.var = BF.genevar[[unique.gene[i]]]
        if (nrow(BF.var) > 0) {
            for (j in 1:nrow(BF.var)) {
                category = data[BF.var$variant[j],]$category
                lkhd[i] = lkhd[i] * ((1 - eta.k[max.iter, category]) + eta.k[max.iter, 
                  category] * BF.var$BF[j])
            }
        }
        teststat[i] = 2 * log((1 - delta.est[max.iter]) + delta.est[max.iter] * 
            lkhd[i])
        # this is the test statistics of one gene
        total.lkhd = total.lkhd + log((1 - delta.est[max.iter]) + delta.est[max.iter] * 
            lkhd[i])
        pvalue[i] = pchisq(teststat[i], 2, lower.tail = F)
    }
    teststat[num.gene + 1] = 2 * total.lkhd
    pvalue[num.gene + 1] = pchisq(teststat[num.gene + 1], 2, lower.tail = F)
    # calculate the LRT statistics and p-value for categories
    if (verbose) 
        cat("Computing LRT statistics and p-values by categories ...\n")
    cate.lkhd = rep(1, num.group)
    sum.lkhd = 0
    cate.stat = numeric()
    cate.pvalue = numeric()
    for (g in 1:num.group) {
        total.lkhd = 0
        lkhd.gene = rep(1, num.gene)
        for (i in 1:num.gene) {
            BF.var = BF.genevar[[unique.gene[i]]]
            if (nrow(BF.var) > 0) 
                for (j in 1:nrow(BF.var)) {
                  if (data[BF.var$variant[j],]$category == g) {
                    mixed_bf = 1 - eta.k[max.iter, g] + eta.k[max.iter, g] * BF.var$BF[j]
                    lkhd.gene[i] = lkhd.gene[i] * mixed_bf
                    cate.lkhd[g] = cate.lkhd[g] * mixed_bf
                  }
                }
            total.lkhd = total.lkhd + log((1 - delta.est[max.iter]) + delta.est[max.iter] * 
                lkhd.gene[i])
        }  # end of i
        cate.stat[g] = 2 * total.lkhd
        cate.pvalue[g] = pchisq(2 * total.lkhd, 1, lower.tail = F)
    }
    sum.lkhd = sum(cate.stat)
    ############################################## 
    return(result = list(delta.est = delta.est[max.iter], delta.pvalue = pvalue[length(pvalue)], 
        eta.est = eta.k[max.iter, ], eta.pvalue = cate.pvalue, 
        BF.gene = data.frame(Gene = unique.gene, BF = BF.gene[max.iter, ]), BF.all = BF.genevar, Eui = EUi))
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