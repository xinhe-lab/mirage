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
#' @return \item{BF.gene}{Bayes factor of genes}
#' \item{delta.est}{Estimate for proportion of risk genes}
#' \item{delta.pvalue}{Significant test for delta = 0}
#' \item{eta.est}{Estimate for proportion of risk variants in a variant group}
#' \item{eta.pvalue}{Significant test for eta = 0}
#' \item{BF.all}{a list of Bayes factor for all variants in a gene}
#' @examples
#' # see example at https://xinhe-lab.github.io/mirage/articles/mwe.html
#' @importFrom progress progress_bar
#' @export
# format of input data column 1: variant ID 2: Gene ID 3 No.variant in cses 4 No.variant in control 5 variant group index 
# n1: sample size in cases n2: sample size in control
mirage=function(data, n1, n2, gamma=3, sigma=2, eta.init=0.1, delta.init=0.1, estimate.delta = TRUE, max.iter = 10000, tol = 1e-05, verbose = TRUE)
{
  # Input check & initialize
  if (ncol(data) == 4) data = cbind(seq(1, nrow(data)), data)
  if (ncol(data) != 5) stop("Input data should have 4 or 5 columns!")
  
  names(data) = c("ID", "Gene", "No.case", "No.contr", "category")
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
    indi.gene.BF=1; var.BF=numeric()
    if (length(var.index.list)>0) # calculate Bayes factor for variant (i,j)
      for (j in 1:length(var.index.list))
      {
        var.index=var.index.list[j]
        category=data$category[var.index]
        var.BF[j]=BF.var.inte(data$No.case[var.index], data$No.contr[var.index], ifelse(length(gamma)>1, gamma[category], gamma), ifelse(length(sigma)>1, sigma[category], sigma), n1, n2) # use uniform sigma/bar.gamma or category specific
        indi.gene.BF=indi.gene.BF*((1-eta.k[1, category])+eta.k[1, category]*var.BF[j])
        
      }
    full.info.genevar[[i]]=cbind(indi.gene, var.BF)
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
  pvalue[num.gene+1]=pchisq(teststat[num.gene+1], 2, lower.tail=F)
  
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
  return(result = list(delta.est = delta.est[max.iter], delta.pvalue = pvalue[length(pvalue)], 
                       eta.est = eta.k[max.iter, ], eta.pvalue = cate.pvalue, 
                       BF.PP.gene = data.frame(Gene = unique.gene, BF = BF.gene[max.iter, ], 
                                               post.prob=(delta.est[max.iter]* BF.gene[max.iter, ])/(delta.est[max.iter]* BF.gene[max.iter, ]+1-delta.est[max.iter])),
                       BF.all = full.info.genevar, Eui = EUi))
}

###############################################################
#' mirage: MIxture model based Rare variant Analysis on GEnes
#'
#' This function implements rare variant test with MIRAGE model for variant set only without gene level information  
#'
#' @param data variant count data, a 4 column data frame for 1) locus ID  2) variant count in cases, 
#' 3) variant count in control and 4) variant category index for a variant. The 1st column is optional.
#' @param n1 sample size in cases.
#' @param n2 sample size in controls.
#' @param gamma a list of category specific hyper prior shape parameter in Beta distribution  for effect size, or a numeric value if all category share the same effect size.
#' @param sigma a list of category specific hyper prior scale parameter in Beta distribution for effect size, or a numeric value if all category share the same effect size.
#' @param eta.init initial value for prior on proportion of risk variants in a variant set.
#' @param max.iter maximum number of iterations enforcing EM algorithm to stop 
#' @param tol threshold of parameter estimate difference to determine the convergence of EM algorithm
#' 
#' 
#' @return \item{full.info}{Bayes factor of individual variant}
#' \item{eta.est}{Estimate for proportion of risk variants in a variant group}
#' \item{eta.pvalue}{Significant test for eta = 0}
#' @examples
#' # see example at https://xinhe-lab.github.io/mirage/articles/mwe.html
#' @importFrom progress progress_bar
#' @export
# format of input data column 1: variant ID 2: No.variant in cses 3 No.variant in control 4 variant group index 
# n1: sample size in cases n2: sample size in control
# this is for variant sets (VS) analysis which may be from multiple variant groups 
mirage_vs=function(data, n1, n2, gamma=3, sigma=2, eta.init=0.1, max.iter = 10000, tol = 1e-05, verbose = TRUE)
{
  # Input check & initialize
  if (ncol(data) == 3) data = cbind(seq(1, nrow(data)), data)
  if (ncol(data) != 4) stop("Input data should have 3 or 4 columns!")
  
  names(data) = c("ID", "No.case", "No.contr", "category")
  groups = unique(data$category)
  num.group = length(groups)
  eta.k = matrix(nrow = max.iter, ncol = num.group)
  data[,4] = match(data[,4], groups)
  colnames(eta.k) = groups
  eta.k[1, ] = rep(eta.init, num.group)
  full.info.var=list()
  num.var=nrow(data)
  var.BF=numeric()
  ########################
  if (verbose && num.var > 1) {
    pb = progress_bar$new(format = "[:spin] Initial analysis of unit :unit out of :total :elapsed", 
                          clear = FALSE, total = num.var, show_after = 0.5) 
  } else { 
    pb = null_progress_bar$new()
  }
  
  
  #########################################
  ########################
  # calculate the Bayes factor for variant j as initials.
  var.index.list=data$category
  if (length(var.index.list)>0) # calculate Bayes factor for variant j
    for (j in 1:length(var.index.list))
    {
      category=var.index.list[j]  # category for variant j
      var.BF[j]=BF.var.inte(data$No.case[j], data$No.contr[j], ifelse(length(gamma)>1, gamma[category], gamma), ifelse(length(sigma)>1, sigma[category], sigma), n1, n2) # use uniform sigma/bar.gamma or category specific
    
    }
  full.info.var=cbind(data, var.BF)
  #########################################
  ########################## EM algorithm
  if (verbose) {
    pb = progress_bar$new(format = "[:spin] Iteration :iteration (diff = :delta) :elapsed", 
                          clear = TRUE, total = max.iter, show_after = 0.5)
  } else { 
    pb = null_progress_bar$new()
  }
  for (iter in 2:max.iter)
  {
    prev_iter=iter-1
    ############################
    ############## EM algorithm: E step
    EZj=numeric() # expectation for variant j
    if (nrow(full.info.var)>0)
      for (j in 1:nrow(full.info.var))
      {
        category=full.info.var$group.index[j]  # category index for variant j 
        
        if (num.group>1)
        {
          numer=full.info.var$var.BF[j]*eta.k[prev_iter, category]
          denom=full.info.var$var.BF[j]*eta.k[prev_iter, category]+(1-eta.k[prev_iter, category])
        }
        
        if (num.group==1)
        {
          numer=full.info.var$var.BF[j]*eta.k[prev_iter]
          denom=full.info.var$var.BF[j]*eta.k[prev_iter]+(1-eta.k[prev_iter])
        }
        
        EZj[j]=numer/denom
      }
    
    ############ EM algorithm: M step
    for (g in 1:num.group)
    {
      var.in.group.index=which(data$category==g)
      if (length(var.in.group.index)>0)
        eta.k[iter, g]=sum(EZj[var.in.group.index])/length(var.in.group.index)
      if (length(var.in.group.index)==0)
        eta.k[iter, g]=0
    }
    ################
    if (num.group>1)
      diff=sum(abs(eta.k[iter,]-eta.k[prev_iter,]))
    if (num.group==1)
      diff=sum(abs(eta.k[iter]-eta.k[prev_iter]))
    
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
    cat("Computing LRT statistics and p-value for every variant and all as a whole ...\n")
  lkhd=rep(1,num.var); total.lkhd=0
  teststat=numeric(); pvalue=numeric()
  pp=numeric() # pp:posterior probability
  
  if (nrow(full.info.var)>0)
    for (j in 1:nrow(full.info.var))
    {
      category=full.info.var$group.index[j]
      if (num.group>1)
      {
        lkhd[j]=lkhd[j]*((1-eta.k[max.iter, category])+eta.k[max.iter, category]*full.info.var$var.BF[j])
        pp[j]=(eta.k[max.iter, category]*full.info.var$var.BF[j])/(eta.k[max.iter, category]*full.info.var$var.BF[j]+1-eta.k[max.iter, category])
      }  
      if (num.group==1)
      {
        lkhd[j]=lkhd[j]*((1-eta.k[max.iter])+eta.k[max.iter]*full.info.var$var.BF[j])
        pp[j]=(eta.k[max.iter]*full.info.var$var.BF[j])/(eta.k[max.iter]*full.info.var$var.BF[j]+1-eta.k[max.iter])
      }  
      
      teststat[j]=2*log(lkhd[j]); # this is the test statistics of one gene
      total.lkhd=total.lkhd+log(lkhd[j])
      
      pvalue[j]=pchisq(teststat[j], num.group, lower.tail=F)
    }
  teststat[num.var+1]=2*total.lkhd
  pvalue[num.var+1]=pchisq(teststat[num.var+1], num.group, lower.tail=F) # 
  
  ##################
  ############################################## calculate category specific test statistics and p value
  if (verbose) 
    cat("Computing LRT statistics and p-values by categories ...\n")
  ##################
  cate.lkhd=rep(1,num.group); cate.stat=numeric()
  cate.pvalue=numeric(num.group); sum.lkhd=0
  if (num.group>1)
    for (g in 1:num.group)
    { # g=2
      if (nrow(full.info.var)>0)
        for (j in 1:nrow(full.info.var))
          if (full.info.var$group.index[j]==g)
            cate.lkhd[g]=cate.lkhd[g]*((1-eta.k[max.iter, g])+eta.k[max.iter, g]*full.info.var$var.BF[j])
          
          cate.stat[g]=2*log(cate.lkhd[g])
          cate.pvalue[g]=pchisq(cate.stat[g], 1, lower.tail=F)
    } # end of g
  if (num.group==1)
  {
    if (nrow(full.info.var)>0)
      for (j in 1:nrow(full.info.var))
        cate.lkhd[1]=cate.lkhd[1]*((1-eta.k[max.iter])+eta.k[max.iter]*full.info.var$var.BF[j])
      
      cate.stat[1]=2*log(cate.lkhd)
      cate.pvalue[1]=pchisq(cate.stat, 1, lower.tail=F)
      
  }
  return(result=list(eta.est=eta.k[max.iter,], full.info=full.info.var, eta.pvalue=cate.pvalue, post.prob=pp))
  
}

#################################### two self defined functions 
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
