intergrand = function(aa, indi.var, pheno, bar.gamma, sig) {
    ff = dbinom(sum(indi.var[pheno == 1]), sum(indi.var), aa * N1/(aa * N1 + N0)) * dgamma(aa, bar.gamma * sig, sig)
    return(ff)
}
# calculate the bayes factor of a single variant via integration
BF.gene.inte = function(geno.var, pheno, bar.gamma, sig) {
    marglik0.CC <- dbinom(sum(geno.var[pheno == 1]), sum(geno.var), N1/(N1 + N0))  # Under H0: gamma=1
    
    marglik1.CC <- integrate(intergrand, geno.var, pheno, bar.gamma, sig, low = 0, upper = 100, stop.on.error = F)$value  # Under H1: gamma~gamma(gamma.mean*sigma, sigma) 
    BF.var <- marglik1.CC/marglik0.CC
    
    return(BF.var)
}
