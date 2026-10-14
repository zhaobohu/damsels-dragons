seed = 67
seqfile = ajun_concat.phy
treefile = ajun_calib.tre
mcmcfile = mcmc.txt
outfile = out.txt

ndata = 3 * three codon positions
seqtype = 0 * nucleotides
usedata = 2 * 0 no likelihood; 1 seq; 2 approximation; 3 out.BV (in.BV)
clock = 2 * 1 global; 2 independent (Newton et al. used this); 3 correlated rates
RootAge = 'B(2.284, 3.199, 0.025, 0.025)' * fallback if no root calibration on tree

model = 4 * 0:JC69, 1:K80, 2:F81, 3:F84, 4:HKY85 (closest to what Newton et al. used for Hessian)
alpha = 0.5 * gamma shape for rate variation across sites
ncatG = 5 * No. categories in discrete gamma

cleandata = 0 * remove sites with ambiguity data (1:yes, 0:no)?

BDparas = 1 1 0.1 C * birth, death, sampling
kappa_gamma = 6 2 * gamma prior on transition/transversion ratio
alpha_gamma = 1 1 * gamma prior on alpha

rgene_gamma = 2 20 1    * gamma-Dirichlet prior on mean substitution rate - mean rate ~ 2/20 = 0.1 substitutions per site per 100 Ma
sigma2_gamma = 1 10 1   * gamma-Dirichlet prior on rate variance across branches

finetune = 1: .1 .1 .1 .1 .1 .1 * auto (0 or 1) : times, rates, mixing...

print = 1 * 0: no mcmc sample; 1: everything except branch 2: ev...
burnin = 200000
sampfreq = 100
nsample = 20000