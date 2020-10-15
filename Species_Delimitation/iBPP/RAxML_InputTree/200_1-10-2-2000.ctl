* iBPP control file written by BPPmulti.pl

seed     = -1
	
seqfile  = ../../L200.txt
Imapfile = ../../imap.txt
traitfile = ../../morpho_ibpp.txt
outfile  = std.theta1-tau2.Tree2.iBPP.intgr.out.txt
mcmcfile = std.theta1-tau2.Tree2.iBPP.intgr.mcmc.txt

speciesdelimitation = 1 1 2 1 0 1
uniformrootedtrees  = 1

species&tree = 9 A B C D E F G H I
 20 1 5 1 1 1 31 1 26
 (A,(((B,C),(D,E)),((I,(H,F)),G)));

usedata      = 1
useseqdata   = 1
usetraitdata = 1
nloci        = 200
ntraits      = 15
nindT        = 320
cleandata    = 0

thetaprior   = 1 10
tauprior     = 2 2000 1
nu0          = 0
kappa0       = 0

heredity   =  ../../
locusrate  =  ../../

finetune     = 1: .03 .04 .05 .06 .07 .01 .01 .01

print        = 1
burnin       = 150000
sampfreq     = 10
nsample      = 1000000
