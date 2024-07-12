# NOTES ON OUR MODIFIED VERSION OF [`bs_inBV`](https://github.com/evolbeginner/bs_inBV)

>> **DISCLAIMER**: The property of this code belongs to Sishuo Wang ([`@evolbeginner`](https://github.com/evolbeginner/bs_inBV)). We cloned the repository [`bs_inBV`](https://github.com/evolbeginner/bs_inBV) in February 2024 and modified the code so that we could run it in a different environment than the one pre-defined by the developer and maintainer of the code back then. What you can find in this directory is the cloned repository with our changes as of 24/02/25. We raised [an issue on March 2024](https://github.com/evolbeginner/bs_inBV/issues/1) to suggest the changes you shall read below, which were then added onto a [pull request on April 2024](https://github.com/evolbeginner/bs_inBV/pull/2). Our changes were accepted and merged with `main`, and so both the issue and the pull request are closed.

## Code changes

Given that we could not run this tool with the scripts available in the [`bs_inBV` GitHub repository](https://github.com/evolbeginner/bs_inBV) as of 24/02/25, we have kept track of the various changes we had to apply to the original code and a list of the dependencies we had to install to make it work. The procedure we followed is detailed below:

* Install not only `parallel` and `colorize`, but also all the ruby gems that are required in the ruby scripts available as part of this tool:
  * parallel
  * colorize
  * find
  * getoptlong
  * time
  * fileutils
  * tmpdir
  * bio
  * bio-nwk
* Command `require_relative` would not work. Consequently, we could not successfully load the functions written in the scripts available within the [`lib`](lib/) directory. To that end, we had to include `$LOAD_PATH << './lib'` as current line 12 in script [`create_hessian_by_bootstrapping.rb`](create_hessian_by_bootstrapping.rb) and replaced commands `require_relative` with `require` to load the ruby scripts in the [`lib`](lib/) directory: [`Dir.rb`](lib/Dir.rb), [`processbar.rb`](lib/processbar.rb), and [`do_mcmctree.rb`] (i.e., current lines 14-16). We also had to include extension `rb` when calling these scripts.
* We had issues when installing `newick-utils` from [the source code](https://github.com/tjunier/newick_utils). Therefore, we searched for pre-compiled binaries and found them [on this website](https://web.archive.org/web/20210409163921/http://cegg.unige.ch/newick_utils). Specifically, we downloaded and installed [Newick Utilities v1.6](https://web.archive.org/web/20210409163921/http://cegg.unige.ch/pub/newick-utils-1.6-Linux-x86_64-disabled-extra.tar.gz).
* R is also required so, if users have not yet installed it via the terminal (i.e., WSL users will need to install R via the terminal as R installed on Windows does not work), please follow the instruction on the `CRAN` to get it installed. Command `Rscript` is also required. Package `MASS`, which is used by script [`reorder_node.rb](reorder_node.rb), is also required.
* We modified the path to PAML programs in line 40 in script [`do_mcmctree.rb`](lib/do_mcmctree.rb). If users have an alias for `MCMCtree`, they may have to change line 43 accordingly. If users have not exported `IQ-TREE` to their PATH or have an alias different from the one specified in line 26 in script [`create_hessian_by_bootstrapping.rb`](create_hessian_by_bootstrapping.rb), they will need to change that too. The same applies to lines 27-28 for `newick-utils` programs and line 29 for `MCMCtree` in that same script. Perhaps the developer could modify the scripts so that these paths were automatically found or add an argument for users to pass their own paths.

Below, you can find the original content of the `README.md` file in the [`bs_inBV` GitHub repository](https://github.com/evolbeginner/bs_inBV).

---

# bs_inBV
A tool to help generate the file in.BV when using "complex" substitution models for MCMCTree. Particularly, the hessian is approximated by a bootstrap approach.

# Idea #
People have increasingly recognized the importance of using more complex substitution models in modeling sequence evolution, particularly for deep-time evolution and in case the sequence divergence is large. MCMCTree supports only a number of substitution models. These models, however, may be available in other software specifically used for tree construction, e.g., IQ-Tree and RAxML.

To work with the approximate likelihood method of MCMCTree (usedata=2), one needs to get the MLE (maximum likelihood estimate) of the branch length, as well as the gradient and hessian evaluated at the MLE. This tool takes the advantages of the abundant substitution models of IQ-Tree. Briefly, the tool uses IQ-Tree's estimate of the branch length, set the gradient to all zeros, and importantly, **approximate the hessian by calculating the negative inverse of the bootstrap covariance matrix of branch length estimates**.
![image](https://github.com/evolbeginner/bs_inBV/assets/8715751/6b7ae95a-f018-4331-8812-720601f637ed)


# Installation
Make sure [RUBY](https://www.ruby-lang.org/en/) is installed.

Many scripts included in this product require the following RUBY packages. If any is not installed, please install it by `gem install package_name`.
* parallel
* colorize

This product also employs several other computational tools. Please ensure that you have them installed.
* [PAML](https://github.com/abacus-gene/paml)
* [IQ-Tree](http://www.iqtree.org/)
* [newick_utilities](https://github.com/tjunier/newick_utils)

# Usage #
1. Perform a regular MCMCTree analysis. To save time, you can run with only the prior (e.g., set usedata = 0 in mcmctree.ctl) and very few iterations (set burnin=1, sampfreq=1, nsample=1 in mcmctree.ctl).
  The reason for this step is to make sure the order of the tips in the file in.BV generated by MCMCTree is the same as that of the species.tree specified by the user.

2. Extract a "reference" tree from the file in.BV generated by MCMCTree in the previous step.
  
    `sed '4!d' in.BV > ref.tre`

3. Run the following

    `ruby create_hessian_by_bootstrapping.rb --ali alignment --calibrated_tree species_tree --outdir C20 --ref ref.tre --force -b 1000 --cpu 8 -m LG+G+C20 --mcmctree_ctl mcmctree.ctl --run_mcmctree --pmsf`
    
    Arguments:
      * `--ali`: alignment file
      * `--calibrated_tree`: the species tree used in regular MCMCTree
      * `--outdir`: output directory
      * `--force`: tells the program to overwrite the output if it exists
      * `-b`: the no. of bootstraps (note that only traditional bootstrapping is allowed as IQ-TRee's UFB cannot be used for a species tree with a fixed topology)
      * `--cpu`: no. of cores used in IQ-Tree
      * `-m`: the model passed to IQ-Tree
      * `--mcmctree_ctl`: the control file for MCMCTree
      * `--pmsf`: use the PMSF approximation

4. Output
  Go to the folder "Sample". Run `ruby create_hessian_by_bootstrapping.rb --ali combined.phy --calibrated_tree species.trees --outdir C20 --ref ref.tre --force -b 100 --cpu 8 -m LG+G+C20 --mcmctree_ctl mcmctree.ctl --run_mcmctree --pmsf`. Note that as detailed above, the files "ref.tre" and "mcmctree.ctl" are obtained by running a regular MCMCTree analysis prior to using create_hessian_by_bootstrapping. The example data are simulated using LG+C20+G assuming a root age of 4.0 Ga. For the details of simulation, see [link](https://sishuowang2022.wordpress.com/2023/06/20/substitution-models-lgcxx-and-cxxlg-differ-in-iq-tree/).
    * C20/split:  IQ-Tree output
    * C20/mcmctree:  MCMCTree results. The branch length estimated by the specified model (in the example "-m LG+G+C20"), gradients, and the approximated hessian by bootstrap in the in.BV

# Notes #
1. With `--run_mcmctree`, MCMCTree will be run directly after generating in.BV. In case you want only the file in.BV based on your specified model say LG+G+C60, please do not use `--run_mcmctree`.
2. Without `--pmsf` IQ-Tree is much slower because it is the traditional Cxx model that will be applied in IQ-Tree.
3. In case of problems due to invertible variance-covariance matrix of the branch lengths, try either increasing bootstrap no. or avoiding too closely related species.
4. For more details, please see **Note S4 and Figs. S9-S10** in the original refenrence (see below).

# How to cite
Dating the bacterial tree of life based on ancient symbiosis Sishuo Wang, Haiwei Luo bioRxiv 2023.06.18.545440; doi: https://doi.org/10.1101/2023.06.18.545440.

You may also need to cite corresponding papers for the use of PAML, IQ-Tree, and Newick_Utilities.
