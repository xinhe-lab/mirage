#' mirage: MIxture model based Rare variant Analysis on GEnes
#'
#' This function implements rare variant test with full MIRAGE model with both variant and gene level 
#'
#' @param data variant count data, a 5 column data frame for 1) locus ID, 2) Gene  3) count in cases, 
#' 4) count in controls and 5) variant category index for a variant. The 1st column is optional.
#' @param n1 sample size in cases.
#' @param n2 sample size in controls.
#' @param gamma a list of category specific hyper prior shape parameter in Beta distribution  for effect size, or a numeric value if all category share the same effect size.
#' @param sigma a list of category specific hyper prior scale parameter in Beta distribution for effect size, or a numeric value if all category share the same effect size.
#' @param delta.init initial value for prior on proportion of risk genes. Must be a positive number between 0 and 1. 
#' @param eta.init initial value for prior on proportion of risk variants in a variant set.
#' @param estimate.delta When TRUE delta is to be estimated and FALSE delta is fixed at delta.init
#' @param max.iter maximum number of iterations enforcing EM algorithm to stop 
#' @param tol threshold of parameter estimate difference to determine the convergence of EM algorithm  
#' 
#' @return \item{BF.PP.gene}{Bayes factor and posterior probability  of genes}
#' \item{delta.est}{Estimate for proportion of risk genes}
#' \item{delta.pvalue}{Significant test for delta = 0}
#' \item{eta.est}{Estimate for proportion of risk variants in a variant group}
#' \item{eta.pvalue}{Significant test for eta = 0}
#' \item{BF.all}{a list of Bayes factor for all variants in a gene}
#' @examples
#' # see example at https://xinhe-lab.github.io/mirage/articles/mwe.html
#' @importFrom progress progress_bar
#' @export
# format of input data column 1: variant ID 2: Gene ID 3 No.variant in cases 4 No.variant in control 5 variant group index 
# n1: sample size in cases n2: sample size in control
mirage=function(data, n1, n2, gamma=3, sigma=2, eta.init=0.1, delta.init=0.1, estimate.delta = TRUE, max.iter = 10000, tol = 1e-05, verbose = TRUE)
{
  # Input check & initialize
  if (ncol(data) == 4) data = cbind(seq(1, nrow(data)), data)
  if (ncol(data) != 5) stop("Input data should have 4 or 5 columns!")
  names(data) = c("ID", "Gene", "No.case", "No.contr", "category")
  # data=data[order(data$category, decreasing = F),]
  # 
  # ################# re-index orignal group index to new natural consecutive index, 1, 2, 3,... 
  # original.group.index=unique(data$category)
  # for (i in 1:length(original.group.index))
  #   data[data$category==original.group.index[i],]$category=i
  # #################
  
  data$category_factor <- factor(data$category, levels = unique(data$category))
  data$category <- as.numeric(data$category_factor)
  original.group.index <- as.numeric(as.character(data[!duplicated(data$category),]$category_factor)) ###
  
  gene.list=data$Gene
  unique.gene = unique(gene.list)
  num.gene = length(unique.gene)
  groups = unique(data$category)
  num.group = length(groups)
  data[,5] = match(data[,5], groups)
  eta.k = matrix(nrow = max.iter, ncol = num.group) 
  colnames(eta.k) = groups
  eta.k[1, ] = rep(eta.init, num.group)
  delta.est = delta.init
  BF.gene = matrix(1, nrow = max.iter, ncol = num.gene)
  BF.genevar = list()
  full.info.genevar=list()
  ########################
  
  if (verbose && num.gene > 1) {
    pb = progress_bar$new(format = "[:spin] Initial analysis of unit :unit out of :total :elapsed", 
                          clear = FALSE, total = num.gene, show_after = 0.5) 
  } else { 
    pb = null_progress_bar$new()
  }
  
  # calculate the Bayes factor for variant (i,j) and gene i as initials.
  for (i in 1:num.gene)
  {
    var.index.list=which(gene.list==unique.gene[i]) # extract all variant index for a gene 
    indi.gene=data[var.index.list,]
    indi.gene.BF=1; var.BF=numeric(); original.index=numeric()
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {
        var.index=var.index.list[j]
        category=data$category[var.index]
        original.index[j]=original.group.index[category]
        var.BF[j]=BF.var.inte(data$No.case[var.index], data$No.contr[var.index], ifelse(length(gamma)>1, gamma[original.index[j]], gamma), ifelse(length(sigma)>1, sigma[original.index[j]], sigma), n1, n2) # use uniform sigma/bar.gamma or category specific
        indi.gene.BF=indi.gene.BF*((1-eta.k[1, category])+eta.k[1, category]*var.BF[j])
        
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF, original.index)
    indi.gene.BF=ifelse(indi.gene.BF==Inf, 3*10^300, indi.gene.BF) # set to the limit value when overflow
    BF.gene[1, i]=indi.gene.BF
    
    pb$tick(tokens = list(total = num.gene, unit = i))
  }  # end of i 
  ########################## EM algorithm
  ######################
  if (verbose) {
    pb = progress_bar$new(format = "[:spin] Iteration :iteration (diff = :delta) :elapsed", 
                          clear = TRUE, total = max.iter, show_after = 0.5)
  } else { 
    pb = null_progress_bar$new()
  }
  for (iter in 2:max.iter)
  {
    prev_iter=iter-1
    ############## EM algorithm: E step
    EUiZij=list() # expectation for variant (i,j), every gene may have varying number of variant
    EUi=numeric() # expectation for gene i.
    total.Zij=matrix(nrow=num.gene, ncol=num.group); total.Ui=matrix(nrow=num.gene, ncol=num.group)  # used to estimate beta
    
    for (i in 1:num.gene)
    {
      info.single.gene=full.info.genevar[[i]] # this is a small matrix for that single gene. each row is one variant
      bb=1
      UiZij=numeric()
      if (nrow(info.single.gene)>0)
        for (j in 1:nrow(info.single.gene))
        {
          category=info.single.gene$category[j]  # category of variant j in gene i 
          numer=info.single.gene$var.BF[j]*eta.k[prev_iter, category]*delta.est[prev_iter]
          denom=(delta.est[prev_iter]+(1-delta.est[prev_iter])/BF.gene[prev_iter,i])*(eta.k[prev_iter, category]*info.single.gene$var.BF[j]
                                                                                      +(1-eta.k[prev_iter, category]))
          UiZij[j]=numer/denom
          bb=bb*((1-eta.k[prev_iter, category])+eta.k[prev_iter, category]*info.single.gene$var.BF[j])
          
        }
      EUiZij[[i]]=UiZij
      bb=ifelse(bb==Inf, 3*10^300, bb) # set to the limit value when overflow
      BF.gene[iter,i]=bb
      EUi[i]=delta.est[prev_iter]*bb/(delta.est[prev_iter]*bb+1-delta.est[prev_iter])
      ######################
      # Note here each gene may have multiple annotation groups
      tt=EUiZij[[i]]
      tt[is.na(tt)]=0
      for (g in 1:num.group)
      {
        total.Zij[i, g]=sum(tt[which(info.single.gene$category==g)])
        total.Ui[i, g]=sum(sum(tt[which(info.single.gene$category==g)]>0)*EUi[i])
      }
      
    }  # end of i
    
    ############## EM algorithm: M step
    delta.est[iter]=ifelse(estimate.delta==T, sum(EUi)/num.gene, delta.init) # either estimate delta when estimate.delta=T or use fixed delta.init when estimate.delta=F
    
    for (g in 1:num.group)
    {
      if (sum(total.Ui[,g])!=0)
        eta.k[iter, g]=sum(total.Zij[,g])/sum(total.Ui[,g])
      if (sum(total.Ui[,g])==0)
        eta.k[iter, g]=0
    }
    ################
    diff = sum(abs(eta.k[iter, ] - eta.k[prev_iter, ]))
    if (diff < tol || iter >= max.iter) {
      pb$tick(max.iter)
      max.iter = iter
      break
    }
    pb$tick(tokens = list(delta = sprintf(diff, fmt = "%#.1e"), iteration = iter))
  } # end of iter
  ######################################################################################################
  ######################################################################################################
  eta.k = eta.k[1:max.iter, , drop=FALSE]
  ################## calculate the likelihood ratio test statistics and p value
  if (verbose) 
    cat("Computing gene level LRT statistics and p-values ...\n")
  
  lkhd=rep(1,num.gene); total.lkhd=0
  teststat=numeric(); pvalue=numeric()
  for (i in 1:num.gene)
  {
    single.gene=full.info.genevar[[i]]
    if (nrow(single.gene)>0)
      for (j in 1:nrow(single.gene))
      {
        category=single.gene$category[j]
        lkhd[i]=lkhd[i]*((1-eta.k[max.iter, category])+eta.k[max.iter, category]*single.gene$var.BF[j])
      }  
    
    teststat[i]=2*log((1-delta.est[max.iter])+delta.est[max.iter]*lkhd[i]); # this is the test statistics of one gene
    total.lkhd=total.lkhd+log((1-delta.est[max.iter])+delta.est[max.iter]*lkhd[i])
    
    pvalue[i]=pchisq(teststat[i], 2, lower.tail=F)
  } # end of i
  teststat[num.gene+1]=2*total.lkhd
  pvalue[num.gene+1]=pchisq(teststat[num.gene+1], (1+num.group), lower.tail=F)
  
  ##################
  # calculate the LRT statistics and p-value for categories
  if (verbose) 
    cat("Computing LRT statistics and p-values by categories ...\n")
  cate.lkhd=rep(1,num.group); cate.stat=numeric()
  cate.pvalue=numeric(num.group); sum.lkhd=0
  for (g in 1:num.group)
  { # g=2
    total.lkhd=0; lkhd.gene=rep(1, num.gene)
    for (i in 1:num.gene)
    {
      single.gene=full.info.genevar[[i]]
      if (nrow(single.gene)>0)
        for (j in 1:nrow(single.gene))
          if (single.gene$category[j]==g)
          {
            lkhd.gene[i]=lkhd.gene[i]*((1-eta.k[max.iter, g])+eta.k[max.iter, g]*single.gene$var.BF[j])
            cate.lkhd[g]=cate.lkhd[g]*((1-eta.k[max.iter, g])+eta.k[max.iter, g]*single.gene$var.BF[j])
          }
      
      total.lkhd=total.lkhd+log((1-delta.est[max.iter])+delta.est[max.iter]*lkhd.gene[i])
    } # end of i
    cate.stat[g]=2*total.lkhd
    cate.pvalue[g]=pchisq(2*total.lkhd, 1, lower.tail=F)
  } # end of g
  sum.lkhd=sum(cate.stat)
  ##############################################
  ## map back the original index 
  colnames(eta.k)=original.group.index
  ##############
  
  return(result = list(delta.est = delta.est[max.iter], delta.pvalue = pvalue[length(pvalue)], 
                       eta.est = eta.k[max.iter, ], eta.pvalue = cate.pvalue, 
                       BF.PP.gene = data.frame(Gene = unique.gene, BF = BF.gene[max.iter, ], 
                                               post.prob=(delta.est[max.iter]* BF.gene[max.iter, ])/(delta.est[max.iter]* BF.gene[max.iter, ]+1-delta.est[max.iter])),
                       BF.all = full.info.genevar, Eui = EUi))
}



BF.var.inte = function(var.case, var.contr, bar.gamma, sig, N1, N0) {
  # Under H0: gamma=1
  log.marglik0.CC <- dbinom(var.case, sum(var.case, var.contr), N1/(N1 + N0),log = T)
  # Under H1: gamma~gamma(gamma.mean*sigma, sigma)
  log.marglik1.CC <- log(integrate(intergrand_log, var.case, var.contr, bar.gamma, sig, N1, 
                                   N0, low = 0, upper = 100, stop.on.error = F)$value)
  log.BF.var <- log.marglik1.CC - log.marglik0.CC
  BF.var <- exp(log.BF.var)
  return(BF.var)
}


intergrand_log = function(aa, var.case, var.contr, bar.gamma, sig, N1, N0) {
  log_ff = dbinom(var.case, sum(var.case, var.contr), aa * N1 / (aa * N1 + N0), log = TRUE) +
    dgamma(aa, bar.gamma * sig, sig, log = TRUE)
  return(exp(log_ff))  # Exponentiating to avoid log of negative values in integration
}
