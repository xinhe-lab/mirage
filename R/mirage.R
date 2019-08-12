#'@title Mixture model based rare variant analysis
#'@export
mirage = function(...) {
    
    while (stop.cond == 0) {
        iter = iter + 1
        ############## EM algorithm: E step
        EUiZij = matrix(0, nrow = num.gene, ncol = m)
        EUi = numeric()
        for (i in 1:num.gene) {
            data = all.data[[i]]
            bb = 1
            if (ncol(data$geno) > 0) 
                for (j in 1:ncol(data$geno)) {
                  orig.var.index = data$var.index[j]
                  for (k in 1:num.group) {
                    from = split.ratio[k] * m + 1
                    to = split.ratio[k + 1] * m
                    if (orig.var.index >= from & orig.var.index <= to) 
                      gp.index = k
                  }
                  numer = BF.var[i, orig.var.index] * beta.k[(iter - 1), gp.index] * delta.est[iter - 1]
                  denom = (delta.est[iter - 1] + (1 - delta.est[iter - 1])/BF.gene[(iter - 1), i]) * (beta.k[(iter - 1), gp.index] * BF.var[i, orig.var.index] + (1 - beta.k[(iter - 1), gp.index]))
                  EUiZij[i, orig.var.index] = numer/denom
                  bb = bb * ((1 - beta.k[(iter - 1), gp.index]) + beta.k[(iter - 1), gp.index] * BF.var[i, orig.var.index])
                }
            BF.gene[iter, i] = bb
            EUi[i] = delta.est[iter - 1] * bb/(delta.est[iter - 1] * bb + 1 - delta.est[iter - 1])
        }
        ############## EM algorithm: M step
        delta.est[iter] = sum(EUi)/num.gene
        # delta.est[iter]=delta
        for (k in 1:num.group) {
            from = split.ratio[k] * m + 1
            to = split.ratio[k + 1] * m
            beta.k[iter, k] = sum(EUiZij[, (from:to)])/sum(actu.no.var[, k] * EUi)
        }
        
        diff = sum(abs(beta.k[iter, ] - beta.k[(iter - 1), ])) + abs(delta.est[iter] - delta.est[iter - 1])
        if (diff < thrshd || iter > (max.iter - 1)) 
            stop.cond = 1
        cat(iter, "th iteration is running", "\t", run, "th run", "\n")
    }
    # 
    if (iter < max.iter) 
        beta.k = beta.k[complete.cases(beta.k), ]
    ################## calculate the likelihood ratio test statistics
    lkhd = matrix(1, nrow = num.gene, ncol = num.group)
    total.lkhd = rep(0, num.group)
    for (i in 1:num.gene) {
        # i=33
        data = all.data[[i]]
        
        for (j in 1:ncol(data$geno)) {
            orig.var.index = data$var.index[j]
            # find the variant group for variant (i,j)
            for (k in 1:num.group) {
                from = split.ratio[k] * m + 1
                to = split.ratio[k + 1] * m
                if (orig.var.index >= from & orig.var.index <= to) 
                  gp.index = k
            }
            lkhd[i, gp.index] = lkhd[i, gp.index] * ((1 - beta.k[(iter - 1), gp.index]) + beta.k[(iter - 1), gp.index] * BF.var[i, orig.var.index])
        }  # end of j 
        total.lkhd = total.lkhd + log((1 - delta.est[iter - 1]) + delta.est[iter - 1] * lkhd[i, ])
    }  # end of i
    all.teststat[run, ] = 2 * total.lkhd
}
