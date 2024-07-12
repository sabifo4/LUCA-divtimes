# Bayesian inference of species divergences | concatenated dataset + cross-bracing B (only fossil-calibrated nodes are cross-braced)

## 1. Setting the file structure to run `MCMCtree`

We have already established our working directory, so we just need to create our file structure for the new analysis using the next code snippet:

```sh
# Run from `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
i=1 # num_aln
num_chains=16 # num chains we will run
mkdir -p pipelines_MCMCtree_fossbraced/{GBM,ILN}
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree_fossbraced/$j/{GBM,ILN}/$i
done
```

The `LUCAdup_arcsin` directory will now have these two additional directories with the corresponding subdirectories:

```text
LUCAdup_arcsin
  |- MCMCtree_fossbraced
  |    |- [1-16] # 16 chains
  |         |- [GBM|ILN]
  |              |- 1 # Only 1 dataset, only 1 dir
  |                
  |- pipelines_MCMCtree_fossbraced
       |- [GBM|ILN]/
```

To prepare the script to run `MCMCtree`, we will use again the bash scripts that we had previously uploaded on the HPC:

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
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree_fossbraced mcmctree4.10.7 $num_chains LUCAdup_arcsin 1
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree_fossbraced mcmctree4.10.7 $num_chains LUCAdup_arcsin 1
# Update parts of the script so that paths to directories are upated
cd ../pipelines_MCMCtree_fossbraced
sed -i 's/trees\/calibrated/trees\/calibrated\_fossbraced/' */*sh
sed -i 's/\/MCMCtree\//\/MCMCtree\_fossbraced\//g' */*sh
# Check that it worked
grep 'calibrated' */*sh
grep 'MCMCtree' */*sh
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior (see [the previous `README.md` file for more details](../02_MCMCtree/README.md#1-setting-the-file-structure-to-run-mcmctree)). Firstly, we will generate a directory where `MCMCtree` will run when sampling from the prior:

```sh
# Run from `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
i=1 # num_aln
num_chains=6
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree_fossbraced_prior/$j/CLK/$i
done
```

Then, we will copy the directory `pipelines_MCMCtree_fossbraced` and generate a copy called `pipelines_MCMCtree_fossbraced_prior` :

```sh
# Run from `LUCAdup_arcsin` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
cp -R pipelines_MCMCtree_fossbraced pipelines_MCMCtree_fossbraced_prior
cd pipelines_MCMCtree_fossbraced_prior
```

We will modify the bash script that will be submitted as a job array following the same approach specified in [the previous `README.md` file](../02_MCMCtree/README.md#1-setting-the-file-structure-to-run-mcmctree):

```sh
# Run from `pipelines_MCMCtree_fossbraced_prior` dir on your HPC.
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
sed -i 's/MCMCtree\_fossbraced/MCMCtree\_fossbraced\_prior/' *sh
sed -i 's/\/GBM\//\/CLK\//' *sh
# Comment soft link
sed -i 's/ln \-s/\#ln \-s/' *sh
# Change number of chains to 6!
sed -i 's/-t 1-16/-t 1-6/' *sh
```

Now, you can check that the lines have been correctly modified:

```sh
# Run from `pipelines_MCMCtree_fossbraced_prior` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
cd ../
grep 'usedata' */*sh
grep 'model'  */*sh
grep 'MCMCtree_fossbraced' */*sh
grep '#$ -t' */*sh
cd CLK
```

## 2. Analyses with `MCMCtree` when sampling from the prior

### Submit jobs in an HPC (prior)

Now, we will be able to run `MCMCtree` first when sampling from the prior (i.e., no data used!) using the code snippet below:

```sh
# Run from `pipelines_MCMCtree_fossbraced_prior/CLK` dir on your HPC.
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
mkdir -p tmp_to_transfer_fossbraced/00_prior
cd tmp_to_transfer_fossbraced
for i in `seq 1 $num_chains`
do
mkdir -p 00_prior/CLK/$i 
# Now, copy the files that are required for sum stats.
# We have run 6 chains for analyses sampling from the prior
printf "\n[[ Copying run "$i" for analyses sampling from the prior ]]\n\n"
cp ../MCMCtree_fossbraced_prior/$i/CLK/1/mcmc.txt 00_prior/CLK/$i
cp ../MCMCtree_fossbraced_prior/$i/CLK/1/*ctl 00_prior/CLK/$i
cp ../MCMCtree_fossbraced_prior/$i/CLK/1/SeedUsed 00_prior/CLK/$i
done
grep 'Species tree for FigTree' -A1 ../MCMCtree_fossbraced_prior/1/CLK/1/out.txt | sed -n '2,2p' > 00_prior/node_tree.tree
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/04_MCMCtree_fossbraced` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
mkdir sum_analyses
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_arcsin/tmp_to_transfer_fossbraced/00_prior .
rsync -avz --copy-link <uname>@<logdetails>:<path>/LUCAdup_arcsin/pipelines_MCMCtree_fossbraced_prior .
# Remove blank output files
rm pipelines_MCMCtree_fossbraced_prior/*/*sh.o*
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

> [!NOTE]
>
> Please note that you will only have the `mcmc.txt` files if (i) you have run the analyses following our guidelines or (ii) you have downloaded [the archive we stored at the University of Bristol data repository (~71Gb)](https://data.bris.ac.uk/datasets/405xnm7ei36d2cj65nrirg3ip/LUCAdivtimes.zip). If the latter, please make sure you follow the instructions in the [main `README.md` file of this repository with regards to how to store them following our file structure](../../README.md#paml-analyses). Once you have the `mcmc.txt` files, you will be able to run the R script aforementioned!

The MCMC diagnostics did not find any of the chains problematic after running [our in-house R script `MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R). Therefore, we used an in-house bash script, [`Combine_MCMC.sh`](scripts/Combine_MCMC.sh), to concatenate all the `mcmc.txt` files for the 6 chains in a unique file.

```sh
# Run from `04_MCMCtree_fossbraced/scripts`
cp Combine_MCMC.sh ../sum_analyses/00_prior
# One argument taken: number of chains
cd ../sum_analyses/00_prior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too.
dataset=$( echo CLK )
./Combine_MCMC.sh $dataset mcmc_files_fosscb_CLK "`seq 1 6`" CLK 20000 Y fosscb_CLK
```

The script above will generate a directory called `mcmc_files_fosscb_CLK` inside the `00_prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. In addition, a directory called `mcmcf4traces_fosscb_CLK` will also be generated so that formatted MCMC files compatible with programs such as `Tracer` can be used to check for chain convergence. A template script to generate the `FigTree.tre` file with this `mcmc.txt` has been saved inside the [`dummy_ctl_files`](dummy_ctl_files) directory.

> [!NOTE]
>
> Please note that we have not been able to upload to this GitHub repository the `mcmcf4traces_fosscb_CLK` directory nor the large concatenated `mcmc.txt` inside `mcmc_files_fosscb_CLK`. If you have not generated them by following our guidelines and/or want to check out results, please download [the archive we stored at the University of Bristol data repository (~71Gb)](https://data.bris.ac.uk/datasets/405xnm7ei36d2cj65nrirg3ip/LUCAdivtimes.zip), which follows the same file structure as in this repository.

We will now create a dummy alignment with only 2 AAs, which is required to generate the `FigTree` files with the mean time estimates obtained when using the concatenated `mcmc.txt` files. In order to do that, we can run the [`Generate_dummy_aln.R`](scripts/Generate_dummy_aln.R). Once you run this script, a new directory called `dummy_aln` will be created, which will contain the input dummy alignment.

We have also generated dummy control file to read the dummy alignment. Additionally, we have enabled option `print = -1`. This print setting lets `MCMCtree` know that an MCMC is not to be run. Instead, `MCMCtree` is told to read the input files (file with the dummy alignment, the calibrated tree file, and the concatenated `mcmc.txt` file) and summarise the samples in the `mcmc.txt` (those that were collected from those chains that passed the filters!). The final mean estimated divergence times and the corresponding CIs will be written in the output `FigTree.tre` file.

```sh
# Run from `sum_analyses/00_prior`
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls LUCAdup_246sp_fosscb_calib_MCMCtree.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_fosscb_CLK
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
printf "\n"
mv FigTree.tre FigTree_fosscb_CLK_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree LUCAdup_dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ LUCAdup\_dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" LUCAdup_dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
printf "\n"
mv FigTree.tre FigTree_fosscb_CLK_95CI.tree
cd $base_dir
```

The next step is to plot the calibration densities VS the marginal densities. We used our in-house R script [`Check_priors_margVScalib.R`](scripts/Check_priors_margVScalib.R) to generate these plots. If you are to run this script with other datasets, however, make sure that your "hard bounds" are not `0.000` in the `Calibnodes_*csv` files (i.e., those input calibration files that are used to match the node labels that are calibrated with the corresponding calibrations). Instead, they should all be changed to `1e-300` (i.e., 1e-300 is used as an equivalent of a hard bound to plot distributions in R; `0` values cannot be used). To make sure this was not affecting our csv files, we ran the following code snippet:

```sh
# Run from `04_MCMCtree_fossbraced/calib_files`
sed -i 's/0\.000/1e\-300/g' *csv
```

Once this script has finished, you will see that a new directory `plots/margVScalib/LUCAdup_fosscb` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots).

Now, once the MCMC diagnostics have finished, you can extract the final data that you can use to write a manuscript as it follows:

```sh
# Run from `04_MCMCtree_fossbraced`
mkdir sum_files_prior
cp -R sum_analyses/00_prior/mcmc_files_*/*tree sum_files_prior/
cp -R sum_analyses/00_prior/CLK/*CLK*/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/margVScalib sum_files_prior/
mkdir sum_files_prior/dupnodes
cp -R plots/dupnodes*pdf sum_files_prior/dupnodes
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

### Submit jobs in an HPC (posterior)

Now that we have verified that there are no issues between the calibration densities and the marginal densities, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# Now, go to directory `LUCAdup_arcsin/pipeline_MCMCtree_fossbraced/GBM` dir on your HPC
# and run the following command. Please change directories until
# you are there.
chmod 775 *sh
qsub pipeline_GBM.sh

# Now, go to directory `pipeline_MCMCtree_fossbraced/ILN` 
# and run the following command. Please change directories until
# you are there.
chmod 775 *sh
qsub pipeline_ILN.sh
```

### Setting the file structure to analyse `MCMCtree` output - posterior

We will now create a directory inside the `sum_analyses` directory to analyse the `MCMCtree` output. Nevertheless, we first need to transfer the data from the cluster to the corresponding directory on our local PC for further analyses

```sh
# Go to your HPC and copy the files that are required for sum stats.
# We have run 16 chains for analyses sampling from the posterior.
# Therefore, `i` will go form 1 to 16
# Run from `LUCAdup_arcsin`
cd tmp_to_transfer_fossbraced
num_chains=16
for i in `seq 1 $num_chains` 
do
mkdir -p 01_posterior/{GBM,ILN}/$i
printf "\n[[ Copying run "$i" for analyses sampling from the posterior ]]\n\n"
cp ../MCMCtree_fossbraced/$i/GBM/1/mcmc.txt 01_posterior/GBM/$i
cp ../MCMCtree_fossbraced/$i/GBM/1/*ctl* 01_posterior/GBM/$i
cp ../MCMCtree_fossbraced/$i/GBM/1/SeedUsed 01_posterior/GBM/$i
cp ../MCMCtree_fossbraced/$i/ILN/1/mcmc.txt 01_posterior/ILN/$i
cp ../MCMCtree_fossbraced/$i/ILN/1/*ctl 01_posterior/ILN/$i
cp ../MCMCtree_fossbraced/$i/ILN/1/SeedUsed 01_posterior/ILN/$i
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/04_MCMCtree_fossbraced` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_arcsin/tmp_to_transfer_fossbraced/01_posterior .
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_arcsin/pipelines_MCMCtree_fossbraced .
# Remove blank output files
rm pipelines_MCMCtree_fossbraced/*/*sh.o*
```

### MCMC diagnostics - posterior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the R script [`MCMC_diagnostics_posterior.R`](scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones used when analysing the chains when sampling from the prior.

> [!NOTE]
>
> Please note that you will only have the `mcmc.txt` files if (i) you have run the analyses following our guidelines or (ii) you have downloaded [the archive we stored at the University of Bristol data repository (~71Gb)](https://data.bris.ac.uk/datasets/405xnm7ei36d2cj65nrirg3ip/LUCAdivtimes.zip). If the latter, please make sure you follow the instructions in the [main `README.md` file of this repository with regards to how to store them following our file structure](../../README.md#paml-analyses). Once you have the `mcmc.txt` files, you will be able to run the R script aforementioned!

Given that no problems have been found with any of the chains we ran, we are ready to concatenate the parameter values sampled across the 16 independent chains we ran:

```sh
# Run from `04_MCMCtree_fossbraced/scripts`
cp Combine_MCMC.sh ../sum_analyses/01_posterior
# One argument taken: number of chains
cd ../sum_analyses/01_posterior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too.
dirname_1=GBM
dirname_2=ILN
./Combine_MCMC.sh $dirname_1 mcmc_files_fosscb_GBM "`seq 1 16`" GBM 20000 Y fosscb_GBM
./Combine_MCMC.sh $dirname_2 mcmc_files_fosscb_ILN "`seq 1 16`" ILN 20000 Y fosscb_ILN
```

Once the scripts above have finished, new directories called `mcmc_files_fosscb_[GBM|ILN]` and `mcmcf4traces_fosscb_[GBM|ILN]` will be created inside `01_posterior/`.

> [!NOTE]
>
> Please note that we have not been able to upload to this GitHub repository the `mcmcf4traces_fosscb_[GBM|ILN]` directory nor the large concatenated `mcmc.txt` inside `mcmc_files_fosscb_[GBM|ILN]`. If you have not generated them by following our guidelines and/or want to check out results, please download [the archive we stored at the University of Bristol data repository (~71Gb)](https://data.bris.ac.uk/datasets/405xnm7ei36d2cj65nrirg3ip/LUCAdivtimes.zip), which follows the same file structure as in this repository.

To map the mean time estimates with the filtered chains, we need to use a control file with `print = -1` that reads the calibrated Newick tree, the dummy alignment we previously generated when analysing the results when sampling from the prior, and the concatenated `mcmc.txt` file with the samples from all chains that passed the filters:

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
name_tt=`ls LUCAdup_246sp_fosscb_calib_MCMCtree.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_fosscb_GBM
printf "[[ Generating tree file for concatenated \"mcmc.txt\" for GBM ... ... ]]\n"
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
mv FigTree.tre FigTree_fosscb_GBM_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree LUCAdup_dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ LUCAdup\_dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" LUCAdup_dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_fosscb_GBM_95CI.tree
printf "\n"
cd $base_dir/mcmc_files_fosscb_ILN
printf "[[ Generating tree file for concatenated \"mcmc.txt\" in "$data"/"$i" for ILN ... ... ]]\n"
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
mv FigTree.tre FigTree_fosscb_ILN_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree LUCAdup_dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ LUCAdup\_dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" LUCAdup_dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_fosscb_ILN_95CI.tree
printf "\n"
cd $base_dir
```

Now, once the MCMC diagnostics have finished, we can run our [in-house R script](scripts/Check_priors_VS_posteriors.R) to plot the posterior distributions against the prior distributions, which can help to better assess how informative the data are and whether there are any serious contradictions between the prior and the posterior distributions.

Lastly, you can extract the final data that we used to write our manuscript as it follows:

```sh
# Run from `04_MCMCtree_fossbraced`
mkdir sum_files_post
cp -R sum_analyses/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/
cp -R sum_analyses/01_posterior/*/*/*all_mean*tsv sum_files_post/
mkdir sum_files_post/fosscb_nodes
cp -R plots/priorVSpost*pdf sum_files_post/fosscb_nodes
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```
