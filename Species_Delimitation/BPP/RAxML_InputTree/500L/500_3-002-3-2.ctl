			seed =  -67812

       seqfile = L500.txt
      Imapfile = imap_red.txt
       outfile = svd_500_3-002-3-2.txt
      mcmcfile = svd_500_3-002-3-2.out

* speciesdelimitation = 0 * fixed species tree
* speciesdelimitation = 1 0 2   * species delimitation rjMCMC algorithm0 and finetune(e)
 speciesdelimitation = 1 1 2 1  * species delimitation rjMCMC algorithm1 finetune (a m)
         speciestree = 0  0.4 0.2 0.1   * speciestree pSlider ExpandRatio ShrinkRatio

   speciesmodelprior = 1  * 0: uniform LH; 1:uniform rooted trees; 2: uniformSLH; 3: uniformSRooted

  species&tree = 9      F   G   C   D   B  I   H   E   A     
                        1   2   2   1   1  2   1   1   2
                       (A,(((B,C),(D,E)),((I,(H,F)),G)));
        diploid =   1  1  1  1  1  1  1  1  1
                  
       usedata = 1  * 0: no data (prior); 1:seq like
         nloci = 500  * number of data sets in seqfile

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

    thetaprior = 3 0.002   # invgamma(a, b) for theta
      tauprior = 3 0.2    # invgamma(a, b) for root tau & Dirichlet(a) for other tau's

*     heredity = 1 4 4
*    locusrate = 1 5

      finetune = 1: .03 .04 .05 .06 .07 .01 .01  # finetune for GBtj, GBspr, theta, tau, mix, locusrate, seqerr

         print = 1 0 0 0   * MCMC samples, locusrate, heredityscalars, Genetrees
        burnin = 50000
      sampfreq = 10
       nsample = 75000