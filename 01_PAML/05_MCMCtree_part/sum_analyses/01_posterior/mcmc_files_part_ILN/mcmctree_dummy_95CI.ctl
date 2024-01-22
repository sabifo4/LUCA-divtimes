          seed = -1
       seqfile = /mnt/c/Users/sandr/Sync/00_Collabs/Ed_LUCA_v2/01_PAML/05_MCMCtree_part/dummy_aln/dummy_aln.aln
      treefile = LUCAdup_dummy_cal.tree
      mcmcfile = mcmc.txt
       outfile = out.txt

         ndata = 5
       seqtype = 2    * 0: nucleotides; 1:codons; 2:AAs
       usedata = 0
                      * 2:approximate likelihood; 3:out.BV (in.BV)
         clock = 1

         model = 3    *     0:poisson, 1:proportional,2:Empirical,3:Empirical+F
                      *     6:FromCodon, 8:REVaa_0, 9:REVaa(nr=189)
         alpha = 0.5  * alpha for gamma rates at sites
         ncatG = 4    * No. categories in discrete gamma

     cleandata = 0    * remove sites with ambiguity data (1:yes, 0:no)?

       BDparas = 1 1 0.1    * birth, death, sampling

   rgene_gamma = 2 2.5 * gammaDir prior for rate for genes
  sigma2_gamma = 2 10  * gammaDir prior for sigma^2     (for clock = 1

         print = -1       * 0: no mcmc sample; 1: everything except branch rates 2: everything
        burnin = 100000
      sampfreq = 1000 
       nsample = 20000
