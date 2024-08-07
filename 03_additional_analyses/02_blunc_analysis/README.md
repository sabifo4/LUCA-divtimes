# Assessing branch length impact on timetree inference

## Using `bs_inBV` to generate `in.BV`

We have decided to use the `bs_inBV` approach [developed by Sishuo Wang](https://github.com/evolbeginner/bs_inBV) to assess the impact that branch lengths may have under a complex model. After installing various dependencies and modifying the tool following the comments in [section "Notes SAC" in their `README.md` file](bs_inBV/README.md#notes-sac), we managed to make it work. You can run the commands below to reproduce our analyses:

```sh
# Run from `02_blunch_analysis`
mkdir 00_inp_data
cd 00_inp_data
cp ../../../00_data_formatting/01_inp_data/LUCAdup_246sp_aln.phy .
cp ../../../00_data_formatting/01_inp_data/LUCAdup_246sp_allcb_calib_MCMCtree.tree .
cp ../../../01_PAML/00_CODEML/control_files/* .
# Update the name of the control file so that the `inBV` tool works
mv prepcodeml.ctl mcmctree.ctl
# Extract a reference tree from an already computed `in.BV` file as required
# by the `bsinBV` tool
sed '4!d' ../../../01_PAML/00_CODEML/out_CODEML/Hessian/1/in.BV > ref.tre
# Call the program after having modified some paths
# We do not want the PMSF approximation nor running `MCMCtree` after obtaining the `in.BV`.
# Therefore, we are not using options `--pmsf` nor `--run mcmctree`.
cd ../bs_inBV
nohup nice -n 19 ruby create_hessian_by_bootstrapping.rb --ali ../00_inp_data/LUCAdup_246sp_aln.phy --calibrated_tree ../00_inp_data/LUCAdup_246sp_allcb_calib_MCMCtree.tree --outdir ../01_out_bsinBV --ref ../00_inp_data/ref.tre --force -b 1000 --cpu 4 -m LG+F+G+C60 --mcmctree_ctl ../00_inp_data/mcmctree.ctl --pmsf
```

Once the analyses have finished, the `in.BV` file can be found inside [01_out_bsinBV/mcmctree/in.BV](01_out_bsinBV/mcmctree/in.BV). Now, we can prepare the directory to run `MCMCtree` on the HPC:

```sh
# Run from `02_PAML`
mkdir -p HPC/LUCAdup_bsinBV
cd HPC/LUCAdup_bsinBV 
i=1
mkdir -p alignments/$i
cp ../../../00_inp_data/LUCAdup_246sp_aln.phy alignments/1
mkdir -p Hessian/$i/
cp ../../../01_out_bsinBV/mcmctree/in.BV Hessian/1/
mkdir control_files
cp ../../../00_inp_data/lg.dat control_files/
cp ../../../00_inp_data/mcmctree.ctl control_files/prepcodeml.ctl # needed for the pipeline we subsequently use!
mkdir -p trees/calibrated
cp ../../../00_inp_data/LUCAdup_246sp_allcb_calib_MCMCtree.tree trees/calibrated
mkdir scripts
cp ../../../../../01_PAML/02_MCMCtree/scripts/*MCMCtree*sh scripts/
# Now, we can transfer this directory to the cluster
cd ../
rsync -avz --copy-links LUCAdup_bsinBV <uname>@<logdetails>:<path_to_dir_where_to_save_wd>
```

## Running `MCMCtree` for timetree inference

### 1. Setting the file structure to run `MCMCtree`

Given that we already have the calibrated tree, the individual gene alignment file, and the corresponding `in.BV` file in the HPC... We have everything we need to run `MCMCtree`!

Firstly, we will create the file structure required for the timetree inference analyses using the following code snippet:

```sh
# Run from `LUCAdup_bsinBV` dir on your HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
i=1 # num_aln
num_chains=16 # num chains we will run
mkdir -p pipelines_MCMCtree/{GBM,ILN}
for j in `seq 1 $num_chains`
do
mkdir -p MCMCtree/$j/{GBM,ILN}/$i
done
# NOTE: We will also have the executable files in this working
# directory!
```

The `LUCAdup_bsinBV` directory  will now have these two additional directories with the corresponding subdirectories:

```text
LUCAdup_bsinBV
  |- MCMCtree
  |    |- [1-16] # 16 chains
  |         |- [GBM|ILN]
  |              |- 1 # Only 1 dataset, only 1 dir
  |                
  |- pipelines_MCMCtree
       |- [GBM|ILN]/
```

You can now run the [`generate_job_MCMCtree.sh` script](scripts/generate_job_MCMCtree.sh), one of our in-house bash scripts mentioned above, using the commands below:

```sh
# Run from the `LUCAdup_bsinBV` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
i=1 # num alignment 
num_chains=16 # num chains
printf "Generating job array for dir "$i" and both clocks ... ...\n\n"
# Arg 1: Number of gene alignments being analysed.
# Arg 2: Clock model (e.g., "GBM", "ILN", or "CLK).
# Arg 3: Number of partitions in the alignment file.
# Arg 4: Path to the pipeline directory.
# Arg 5: Command to execute MCMCtree (e.g., `mcmctree`, `mcmctree_2000`, etc.)
# Arg 6: Number of chains run.
# Arg 7: Name of working directory (e.g., `LUCAdup_bsinBV`)
# Arg 8: Boolean, enable duplication option? 0: no, 1: yes
# 
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree4.10.7 $num_chains LUCAdup_bsinBV 1
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree4.10.7 $num_chains LUCAdup_bsinBV 1
# Check paths were correctly incorporated!
cd ../pipelines_MCMCtree
grep 'MCMCtree' */*sh
grep 'alignments' */*sh
grep 'Hessian' */*sh
grep 'ctl_dir' */*sh
```

Given that there are no changes to the tree topology that we used in our previous analyses, there is no need to run the analyses when sampling from the prior (i.e., data are not used, and hence we would be sampling from the same target distribution from which we have already collected samples from).

### 2. Analyses with `MCMCtree` when sampling from the posterior

#### Submit jobs in an HPC (posterior)

We already verified that there are no issues between the calibration densities and marginal densities, so we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# Now, go to directory `LUCAdup_bsinBV/pipelines_MCMCtree/GBM` dir on your HPC
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

#### Setting the file structure to analyse `MCMCtree` output - posterior

We will now create a directory inside the `sum_analyses` directory to analyse the `MCMCtree` output. Nevertheless, we first need to transfer the data from the cluster to the corresponding directory on our local PC for further analyses:

```sh
# Go to your HPC and copy the files that are required for sum stats.
# We have run 16 chains for analyses sampling from the posterior.
# Therefore, `i` will go form 1 to 16
# Run from `LUCAdup_bsinBV`
mkdir tmp_to_transfer
cd tmp_to_transfer
num_chains=16
for i in `seq 1 $num_chains` 
do
mkdir -p 01_posterior/{GBM,ILN}/$i
printf "\n[[ Copying run "$i" for analyses sampling from the posterior ]]\n\n"
cp ../MCMCtree/$i/GBM/1/mcmc.txt 01_posterior/GBM/$i
cp ../MCMCtree/$i/GBM/1/*ctl* 01_posterior/GBM/$i
cp ../MCMCtree/$i/GBM/1/SeedUsed 01_posterior/GBM/$i
cp ../MCMCtree/$i/ILN/1/mcmc.txt 01_posterior/ILN/$i
cp ../MCMCtree/$i/ILN/1/*ctl 01_posterior/ILN/$i
cp ../MCMCtree/$i/ILN/1/SeedUsed 01_posterior/ILN/$i
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `02_PAML/MCMCtree` dir on your local
# PC. Please change directories until
# you are there. Then run the following
# commands.
mkdir sum_analyses
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_bsinBV/tmp_to_transfer/01_posterior .
rsync -avz --copy-links <uname>@<logdetails>:<path>/LUCAdup_bsinBV/pipelines_MCMCtree .
# Remove blank output files
rm pipelines_MCMCtree/*/*sh.o*
```

#### MCMC diagnostics - posterior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the R script [`MCMC_diagnostics_posterior.R`](scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones used when analysing the chains when sampling from the prior.

> [!NOTE]
>
> Please note that you will only have the `mcmc.txt` files if (i) you have run the analyses following our guidelines or (ii) you have downloaded [the archive we stored at the University of Bristol data repository (~71Gb)](https://data.bris.ac.uk/datasets/405xnm7ei36d2cj65nrirg3ip/LUCAdivtimes.zip). If the latter, please make sure you follow the instructions in the [main `README.md` file of this repository with regards to how to store them following our file structure](../../README.md#paml-analyses). Once you have the `mcmc.txt` files, you will be able to run the R script aforementioned!

Given that no problems have been found with any of the chains we ran, we are ready to concatenate the parameter values sampled across the 16 independent chains we ran:

```sh
# Run from `01_MCMCtree/scripts`
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
./Combine_MCMC.sh $dirname_1 mcmc_files_GBM "`seq 1 16`" GBM 20000 Y bsinBV_GBM
./Combine_MCMC.sh $dirname_2 mcmc_files_ILN "`seq 1 16`" ILN 20000 Y bsinBV_ILN
```

Once the scripts above have finished, new directories called `mcmc_files_[GBM|ILN]` and `mcmcf4traces_bsinBV_[GBM|ILN]` will be created inside `01_posterior/`.

> [!NOTE]
>
> Please note that we have not been able to upload to this GitHub repository the `mcmcf4traces_bsinBV_[GBM|ILN]` directory nor the large concatenated `mcmc.txt` inside `mcmc_files_[GBM|ILN]`. If you have not generated them by following our guidelines and/or want to check out results, please download [the archive we stored at the University of Bristol data repository (~71Gb)](https://data.bris.ac.uk/datasets/405xnm7ei36d2cj65nrirg3ip/LUCAdivtimes.zip), which follows the same file structure as in this repository.

To map the mean time estimates with the filtered chains, we need to use a control file with `print = -1` that reads the calibrated Newick tree, the dummy alignment we previously generated when analysing the results when sampling from the prior, and the concatenated `mcmc.txt` file with the samples from all chains that passed the filters:

```sh
# Run from `sum_analyses_prot/01_posterior` directory.
# Please change directories until
# you are there. Then run the following
# commands.
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls LUCAdup_246sp_allcb_calib_MCMCtree.tree`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_GBM
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
mv FigTree.tre FigTree_GBM_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_GBM_95CI.tree
printf "\n"
cd $base_dir/mcmc_files_ILN
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
mv FigTree.tre FigTree_ILN_95HPD.tree
cp $tt_dir/LUCAdup_246sp_uncalib.tree dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'U\(4\.520\,1e\-300\)\'\;/" dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_ILN_95CI.tree
printf "\n"
cd $base_dir
```

Now, once the MCMC diagnostics have finished, we can run our [in-house R script](scripts/Check_priors_VS_posteriors.R) to plot the posterior distributions against the prior distributions, which can help to better assess how informative the data are and whether there are any serious contradictions between the prior and the posterior distributions.

Lastly, you can extract the final data that we used to write our manuscript as it follows:

```sh
# Run from `MCMCtree`
mkdir sum_files_post
cp -R sum_analyses/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/
cp -R sum_analyses/01_posterior/*/*/*all_mean*tsv sum_files_post/
mkdir sum_files_post/plots_cbnodes
cp -R plots/priorVSpost*pdf sum_files_post/plots_cbnodes
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```
