# Bayesian inference of species divergences | partitioned dataset + cross-bracing A (all nodes with same speciation events are mirrored)

## 1. Setting the file structure to run `MCMCtree`

Given that we already generated the calibrated tree, the partitioned alignment file, and we have just generated the `in_5parts.BV` file... We have everything we need to run `MCMCtree`!

First, we will create the file structure required for the timetree inference analyses using the following code snippet:

```sh
# Run from `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
i=1 # num_aln
num_chains=16 # num chains we will run
mkdir -p pipelines_MCMCtree_part/{GBM,ILN}
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree_part/$j/{GBM,ILN}/$i
done
```

The `LUCAdup_arcsin` directory  will now have these two extra directories with the corresponding subdirectories:

```text
LUCAdup_arcsin
  |- MCMCtree_part
  |    |- [1-16] #16 chains
  |         |- [GBM|ILN]
  |              |- 1 # Only 1 partitioned dataset, only 1 dir
  |                
  |- pipelines_MCMCtree_part
       |- [GBM|ILN]/
```

>**IMPORTANT NOTE**: When sampling from the posterior, the likelihood is being calculated or approximated, depending on the `userdata` option you set in the control file to run `MCMCtree`. In other words, the larger the dataset, the more time it will take for `MCMCtree` to finish.

Now, we will transfer the partitioned alignment (i.e., five alignment blocks, one per gene alignment) to the server, which is the file that `MCMCtree` will use for timetree inference:

```sh
# Run from `00_data_formatting/01_inp_data`
rsync -avz --copy-links LUCAdup_246sp_5parts_aln.phy <uname>@<logdetails>:<path>/LUCAdup_arcsin/alignments
```

To prepare the script to run `MCMCtree`, we will use our in-house script, which should already be in the HPC -- they are the same scripts we used for the analyses with the concatenated dataset! You can run the bash script mentioned above using the commands below:

```sh
# Run from the `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
i=1 # num alignment
printf "Generating job array for dir "$i" and both clocks ... ...\n\n"
# Arg 1: Number of gene alignments being analysed.
# Arg 2: Clock model (e.g., "GBM", "ILN", or "CLK).
# Arg 3: Number of partitions in the alignment file.
# Arg 4: Path to the pipeline directory.
# Arg 5: Command to execute MCMCtree (e.g., `mcmctree`, `mcmctree_2000`, etc.)
# Arg 6: Number of chains run.
# Arg 7: Name of working directory (e.g., `LUCAdup_arcsin`)
# Arg 8: Boolean, enable duplication option? 0: no, 1: yes
# 
# Please modify if you are running different analyses
num_chains=16
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree_part mcmctree4.10.7 $num_chains LUCAdup_arcsin 1
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree_part mcmctree4.10.7 $num_chains LUCAdup_arcsin 1
# Now, we need to modify the `ndata` variable as we will be using a partitioned alignment (i.e., 5 alignment blocks, `ndata = 5`)
cd ../pipelines_MCMCtree_part
sed -i 's/ndata\\ \\\=\\ 1/ndata\\ \\\=\\ 5/' */*sh
# Now, change the alignment path too!
sed -i 's/alignments\/\$dir/alignments/' */*sh
# Modify `MCMCtree` paths to save results 
sed -i 's/MCMCtree/MCMCtree\_part/' */*sh
# Modify path to `Hessian_part` and in.BV name
sed -i 's/Hessian\/\$dir/Hessian\_part/' */*sh
sed -i 's/echo in\.BV/echo in\_5parts\.BV/' */*sh
# Check that it worked
grep 'MCMCtree' */*sh
grep 'alignments' */*sh
grep 'Hessian' */*sh
grep 'in_5parts' */*sh
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior. In this way, we can verify whether the user-specified priors constraining some node ages in our fixed tree topology (i.e., probability distributions that we want the dating program to use) and the effective prior (i.e., probability distributions the dating program will use after building the joint prior distribution, which considers not only the user-specified priors but also additional priors such as the birth-death process for uncalibrated nodes, the rate prior, etc.). Oftentimes, truncation issues may arise when the effective priors are in disagreement with the user-specified priors (see an extensive study about this effect in [dos Reis et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26603774/)). To assess the effect of truncation prior to analysing our dataset, we will run `MCMCtree` without the data (i.e., `MCMCtree` will be sampling from the prior distribution, and thus sequence data will not be used).

First, we will generate a directory where `MCMCtree` will run when sampling from the prior:

```sh
# Run from `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
i=1 # num_aln
num_chains=6
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree_part_prior/$j/CLK/$i
done
```

>**IMPORTANT NOTE**: When sampling from the prior, the likelihood is not being calculated or estimated. In other words, the most time-consuming part of the MCMC does not take place. In that way, you should be able to gather enough samples with fewer runs than those needed when sampling from the posterior.

Then, we will copy the directory `pipelines_MCMCtree` and will generate a copy called `pipelines_MCMCtree_prior` :

```sh
# Run from `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
cp -R pipelines_MCMCtree_part pipelines_MCMCtree_part_prior
cd pipelines_MCMCtree_part_prior
```

We will modify the bash script that will be submitted as a job array so that the `userdata` option in the control file is equal to `0`, which enables `MCMCtree` to sample from the prior instead of sampling from the posterior (i.e., the alignment file is ignored). In that way, the rest of the setting concerning the evolutionary model will not be enabled by `MCMCtree` as the sequence data are not being used. Last, we will change the path to where the results will be stored:

```sh
# Run from `pipelines_MCMCtree_part_prior` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
# 
# Prepare directories to sample from the prior,
# only one needed as nucleotide subsitution models 
# are not used.
rm -r ILN 
mv GBM CLK
cd CLK
mv pipeline_GBM.sh pipeline_CLK.sh

# Modify bash script: options `usedata` and `model`
sed -i "s/..*usedata..*/sed \-i \'s\/usedata\.\.\*\/usedata \= 0\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"\nsed \-i \'s\/model\.\.\*\/model \= 3\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"/" *sh
# Modify clock model 
sed -i "s/^fi/fi\nsed \-i \'s\/clock\.\.\*\/clock \= 1\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"/" *sh

# Modify path to save results 
sed -i 's/MCMCtree\_part/MCMCtree\_part\_prior/' *sh
sed -i 's/\/GBM\//\/CLK\//' *sh
# Comment soft link
sed -i 's/ln \-s/\#ln \-s/' *sh
# Change number of chains to 6!
sed -i 's/-t 1-16/-t 1-6/' *sh
```

Now, you can check that the lines have been correctly modified:

```sh
# Run from `pipelines_MCMCtree_part_prior` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
cd ../
grep 'trees' */*sh
grep 'alignments' */*sh
grep 'usedata' */*sh
grep 'model'  */*sh
grep 'ndata' */*sh
grep 'MCMCtree_part_prior' */*sh
grep '#$ -t' */*sh
cd CLK
```

## 2. Analyses with `MCMCtree` when sampling from the prior

### Submit jobs in an HPC (prior)

Now, we will be able to run `MCMCtree` first when sampling from the prior (i.e., no data used!) using the code snippet below:

```sh
# Run from `pipelines_MCMCtree_part_prior/CLK` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
chmod 775 *sh
qsub pipeline_CLK.sh
```

### Setting the file structure to analyse `MCMCtree` output - prior

We will now create a `sum_analyses` directory to analyse the `MCMCtree` output. Nevertheless, we first need to transfer the data from the cluster to the corresponding directory on our local PC for further analyses:

```sh
# Run everything from `LUCAdup_arcsin` in your HPC
num_chains=6
mkdir -p tmp_to_transfer_part/{00_prior,01_posterior} # We create a `01_posterior` for later!
cd tmp_to_transfer_part
for i in `seq 1 $num_chains`
do
mkdir -p 00_prior/CLK/$i 
# Now, copy the files that are required for sum stats.
# We have run 6 chains for analyses sampling from the prior
printf "\n[[ Copying run "$i" for analyses sampling from the prior ]]\n\n"
cp ../MCMCtree_part_prior/$i/CLK/1/mcmc.txt 00_prior/CLK/$i
cp ../MCMCtree_part_prior/$i/CLK/1/*ctl 00_prior/CLK/$i
cp ../MCMCtree_part_prior/$i/CLK/1/SeedUsed 00_prior/CLK/$i
done
grep 'Species tree for FigTree' -A1 ../MCMCtree_part_prior/1/CLK/1/out.txt | sed -n '2,2p' > 00_prior/node_tree.tree
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/05_MCMCtree_part` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands. If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
mkdir sum_analyses
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_arcsin/tmp_to_transfer_part/00_prior .
rsync -avz --copy-link <uname>@<logdetails>:<path>/LUCAdup_arcsin/pipelines_MCMCtree_part_prior .
# Remove blank output files
rm pipelines_MCMCtree_part_prior/*/*sh.o*
```

### MCMC diagnostics - prior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run MCMC diagnostics!

We are going to run the R script [`MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R) and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same node sin the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filters.
5. Generate a new convergence plot with those chains that passed the filters.
6. Calculate Rhat, tail-ESS, and bulk-ESS to check whether chain convergence has been reached with the chains that have passed filters.

The MCMC diagnostics did not find any of the chains problematic after running [our in-house R script `MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R). Therefore, we used [our in-house bash script `Combine_MCMC.sh`](scripts/Combine_MCMC.sh) to concatenate all the `mcmc.txt` files for the 6 chains in a unique file.

```sh
# Run from `05_MCMCtree_part/scripts`
cp Combine_MCMC.sh ../sum_analyses/00_prior
# One argument taken: number of chains
cd ../sum_analyses/00_prior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
dataset=$( echo CLK )
./Combine_MCMC.sh $dataset mcmc_files_part_CLK "`seq 1 6`" CLK 20000
```

The script above will generate a directory called `mcmc_files_part_CLK` inside the `00_prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. A template script to generate the `FigTree.tre` file with this `mcmc.txt` has been saved inside the [`dummy_ctl_files`](dummy_ctl_files) directory.

We will now create a dummy alignment with only 2 AAs to generate the `FigTree` files using the concatenated `mcmc.txt` files. In order to do that, we can run the [`Generate_dummy_aln.R`](scripts/Generate_dummy_aln.R). Once you run it, a new directory called `dummy_aln` will be created, which will contain the dummy alignment.

We have also generated dummy control file with option `print = -1`, which will not run an MCMC but, instead, will use the input files (file with the dummy alignment, calibrated tree file, and concatenated `mcmc.txt` file) to generate a `FigTree.tre` file with the mean estimated divergence times and the corresponding mean CIs using all the samples collected during all the MCMCs.

```sh
# Run from `sum_analyses/00_prior`
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls LUCAdup_246sp_allcb_calib_MCMCtree.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_part_CLK
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_part_CLK_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree LUCAdup_dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ LUCAdup\_dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" LUCAdup_dummy_cal.tree
mcmctree49h_sum95CI *ctl mcmctree_dummy_95CI.ctl
printf "\n"
mv FigTree.tre FigTree_part_CLK_95CI.tree
cd $base_dir
```

The next step is to plot the user-specified prior VS the effective prior. We used our in-house R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R) to generate these plots. If you are to run this script with other datasets, however, make sure that your "hard bounds" are not `0.000` in the `Calibnodes_*csv` files and, instead, they are `1e-300` (i.e., while 1e-300 is rounded to `0.000` in the `MCMCtre` output, which can be used to generate the csv files aforementioned, we need `1e-300` to plot distributions in R). To make sure this was not affecting our csv files, we ran the following code snippet:

```sh
# Run from `05_MCMCtree_part/calib_files`
sed -i 's/0\.000/1e\-300/g' *csv
```

Once this script has finished, you will see that a new directory `plots/effVSuser/LUCAdup_part` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots).

Now, once the MCMC diagnostics have finished, you can extract the final data that you can use to write a manuscript as it follows:

```sh
# Run from `05_MCMCtree_part`
mkdir sum_files_prior
cp -R sum_analyses/00_prior/mcmc_files_part_CLK/*CLK*tree sum_files_prior/
cp -R sum_analyses/00_prior/CLK/*CLK*/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/effVSuser sum_files_prior/
mkdir sum_files_prior/dupnodes
cp -R plots/dupnodes*pdf sum_files_prior/dupnodes
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

### Submit jobs in an HPC (posterior)

Now that we have verified that there are no issues between the user-specified prior and the effective prior, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# Now, go to directory `LUCAdup_arcsin/pipelines_MCMCtree/GBM` dir on your HPC
# and run the following command. Please change directories until
# you are there.
chmod 775 *sh
qsub pipeline_GBM.sh

# Now, go to directory `pipelines_MCMCtree/ILN` 
# and run the following command. Please change directories until
# you are there.
chmod 775 *sh
qsub pipeline_ILN.sh
```

### Setting the file structure to analyse `MCMCtree` output - posterior

We will now go to the previously created `sum_analyses` directory to analyse the `MCMCtree` output. First, we need to come back to the `01_PAML/01_MCMCtree` directory and run the following code snippet:

```sh
# Go to your HPC and copy the files that are required for sum stats.
# We have run 16 chains for analyses sampling from the posterior.
# Therefore, `i` will go form 1 to 16
# Run from `LUCAdup_arcsin`
cd tmp_to_transfer_part
num_chains=16
# The `01_posterior` should already exist from the previous
# analyses. If not, it will be created during the `for` loop
for i in `seq 1 $num_chains` 
do
mkdir -p 01_posterior/{GBM,ILN}/$i
printf "\n[[ Copying run "$i" for analyses sampling from the posterior ]]\n\n"
cp ../MCMCtree_part/$i/GBM/1/mcmc.txt 01_posterior/GBM/$i
cp ../MCMCtree_part/$i/GBM/1/*ctl* 01_posterior/GBM/$i
cp ../MCMCtree_part/$i/GBM/1/SeedUsed 01_posterior/GBM/$i
cp ../MCMCtree_part/$i/ILN/1/mcmc.txt 01_posterior/ILN/$i
cp ../MCMCtree_part/$i/ILN/1/*ctl 01_posterior/ILN/$i
cp ../MCMCtree_part/$i/ILN/1/SeedUsed 01_posterior/ILN/$i
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/05_MCMCtree_part` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands. If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_arcsin/tmp_to_transfer_part/01_posterior .
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_arcsin/pipelines_MCMCtree_part .
# Remove blank output files
rm pipelines_MCMCtree_part/*/*sh.o*
```

### MCMC diagnostics - posterior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the R script [`MCMC_diagnostics_posterior.R`](scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones used when analysing the chains when sampling from the prior. Given that no problems have been found with any of the chains we ran, we are ready to concatenate the parameter values sampled across the 16 independent chains we ran:

```sh
# Run from `05_MCMCtree_part/scripts`
cp Combine_MCMC.sh ../sum_analyses/01_posterior
# One argument taken: number of chains
cd ../sum_analyses/01_posterior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
dirname_1=GBM
dirname_2=ILN
./Combine_MCMC.sh $dirname_1 mcmc_files_part_GBM "`seq 1 16`" GBM 20000
./Combine_MCMC.sh $dirname_2 mcmc_files_part_ILN "`seq 1 16`" ILN 20000
```

Once the scripts above have finished, a new directory called `mcmc_files_part_[GBM|ILN]` will be created inside `01_posterior/`, respectively. To map the mean time estimates with the filtered chains, we need to copy a control file, the calibrated Newick tree, and the dummy alignment we previously generated when analysing the results when sampling from the prior:

```sh
# Run from `sum_analyses_prot/01_posterior` directory.
# Please change directories until
# you are there. Then run the following
# commands.
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls LUCAdup_246sp_allcb_calib_MCMCtree.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_part_GBM
printf "[[ Generating tree file for concatenated \"mcmc.txt\" for GBM ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
sed -i 's/ndata\ \=\ 1/ndata\ \=\ 5/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_GBM_part_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree LUCAdup_dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ LUCAdup\_dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" LUCAdup_dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_GBM_part_95CI.tree
printf "\n"
cd $base_dir/mcmc_files_part_ILN
printf "[[ Generating tree file for concatenated \"mcmc.txt\" in "$data"/"$i" for ILN ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
sed -i 's/ndata\ \=\ 1/ndata\ \=\ 5/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_ILN_part_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree LUCAdup_dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ LUCAdup\_dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" LUCAdup_dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_ILN_part_95CI.tree
printf "\n"
cd $base_dir
```

Now, once the MCMC diagnostics have finished, we can run our [in-house R script](scripts/Check_priors_VS_posteriors.R) to plot the posterior distributions against the prior distributions, which can help to better assess how informative the data are and whether there are any serious contradictions between the prior and the posterior distributions.

Last, you can extract the final data that we used to write our manuscript as it follows:

```sh
# Run from `05_MCMCtree_part`
mkdir sum_files_post
cp -R sum_analyses/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/
cp -R sum_analyses/01_posterior/*/*/*all_mean*tsv sum_files_post/
mkdir sum_files_post/plots_cbnodes
cp -R plots/priorVSpost*pdf sum_files_post/plots_cbnodes
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```
