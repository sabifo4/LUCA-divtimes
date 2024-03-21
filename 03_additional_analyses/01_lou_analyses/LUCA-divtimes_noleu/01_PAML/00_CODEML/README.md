# `CODEML` analysis - concatenated alignment

## 1. Pick rate prior

We will use a vague gamma distribution centered on a mean evolutionary rate estimated by considering the tree height (molecular distance in substitutions per site) and the mean age for the root of our phylogeny (time unit = 100 Mya). As the [rooted phylogeny that we inferred with `IQ-TREE 2`](../../00_data_formatting/00_raw_data/trees/00_IQTREE/LUCAdup_topo_bl_rooted.tree) has information about the branch lengths, we can USE [our R in-house script](scripts/calculate_rateprior.R) to calculate the tree height. We also have a calibration to constrain the root age, which we will use as a rough estimation of the age of the root of our phylogeny based on a geological event: the moon forming impact. Specfically, this calibration is an upper bound that constrains the maximum age of the root to be 4,520 Ma (i.e., 4.520 (time unit = 1 Ga = 100 Mya)), which we will use as an approximate age for the root of our phylogeny to estimate the mean evolutionary rate.

By setting a vague shape ($\alpha=2$) for the gamma distribution that we will use as a rate prior, we can account for the uncertainty surrounding our mean rate estimate. If we had more knowledge on the mean rate, however, we should use a narrower prior with a larger $\alpha$ that better represents our prior information.

Now, we have all the information we need to calculate the $\beta$ parameter for the gamma distribution. We have written the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R) to carry out all the tasks mentioned above. You can open this file in RStudio to find out the value of $\beta$ and plot the final prior on the rates. A summary of what you will find in the script is described below:

```text
First, we know that the molecular distance (tree height, distance from the root to present time) is equal to the mean evolutionary rate (in substitutions per site per time unit) times the age of the divergence time at the root (in time unit, which we can define later). If we have estimated our phylogeny, and therefore have estimated the value of the branch lengths, we will be able to calculate the tree height. The units of the tree height will be the following:

tree_height = rate * root_age --> units_tree_height = subst/site/time * time = subst/site

There are various functions we can use to calculate the tree heigt. We have chosen the R function `phytools::nodeHeights`. The maximum height calculated by this function corresponds to the length from the root to the heighest tip.

After calculating the tree height of our phylogeny (in subst/site) and considering the age of the root based on the fossil record or geological events (time unit = 1 Ga = 100 Mya = 1e9 years), we can get a rough estimate for the mean rate:

Time unit = 1 Ga (mean root age in Ga) --> mean_rate = tree_height / root_age = (subst/site) / (Ga) = subst/site/Ga (time unit = 1 Ga = 1e9 years) 

We also know that the mean of the gamma distribution that we want to use as rate prior is our parameter of interest: the mean evolutionary rate. Therefore:

mean_Gamma = mean_rate = alpha / beta 
Time unit = 1 Ga: mean_rate = alpha / beta --> beta = alpha / mean_rate = 2 / mean_rate

The calibrated tree needs to incorporate the age constraints in the same time unit used to infer the mean evolutionary rate and establish the rate prior (i.e., do not forget to scale the calibrations accordingly if needed!). 
```

If you run the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R), you will see how all the steps described above take place and a new PDF file with the prior distribution to be used will be generated in a new directory called `out_RData`.

A [template control file](control_files/prepcodeml.ctl) with the $\alpha$ and $\beta$ parameters (as defined using the R script above) for the gamma distribution as a prior on the rates is also generated. Note that several options will be subsequently modified to fit the analysis with this dataset (i.e., you will see some options that have flags in capital letters, which will be replaced with the correct value for said option). Given how shallow this tree is, the clock may be expected to be seriously violated, and thus we have fixed a mean for the `sigma2` paramter (i.e., variation in the clock) as 0.1 using a gamma prior with $\alpha=1$ and $\beta=10$: `sigma2_gamma 1 10​`.

## 2. Set up the file structure

Before running `MCMCtree` using the approximate likelihood calculation to speed up timetree inference, we first need to calculate the vector of branch lengths, the gradient (vector), and the Hessian (matrix). We will use `CODEML` for this purpose as our dataset is an amino acid alignment.

The file structure we will use is the following:

```text
LUCAdup_noleu/
  |
  |- alignments/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments
  |       
  |- control_files/ # Pre-defined control file with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- X # Directory for alignment X -- have as many directories as alignments
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `CODEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `CODEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree for `MCMCtree`
      |- uncalibrated # Directory with the uncalibrated tree for `CODEML`
```

To create the `LUCAdup_noleu` file structure, we run the following commands from the PC before transferring to the HPC:

```sh
# Run the following commands from 
# directory `00_CODEML`
mkdir -p HPC/LUCAdup_noleu
cd HPC/LUCAdup_noleu 
i=1
mkdir -p alignments/$i
mkdir -p Hessian/$i/prepare_codeml
mkdir -p pipelines_Hessian
mkdir control_files
mkdir -p trees/{uncalibrated,calibrated}
mkdir scripts
```

Once the file structure is created, we can now populate it with the input files we have generated some minutes ago: alignment file, tree files, and control file. We will also add the `lg.dat` file, which has the matrix to enable `CODEML` to use the LG substitution model. You can transfer the files to the HPC as you prefer (below, we show an example of how to do this with `rsync`, a system that we will keep as an example throughout the tutorial):

```sh
# Run from `HPC/LUCAdup_noleu`
# Copy alignment
cp ../../../../00_data_formatting/01_inp_data/LUCAdup_243sp_aln.phy alignments/1/
# Now, transfer the trees
cp ../../../../00_data_formatting/01_inp_data/LUCAdup_243sp_allcb_calib_MCMCtree.tree trees/calibrated/
cp ../../../../00_data_formatting/01_inp_data/LUCAdup_243sp_uncalib.tree trees/uncalibrated/
# Next, copy control files
cp ../../control_files/* control_files/
# Last, copy the in-house bash scripts with our pipeline
cp ../../scripts/*sh scripts/
# Once everything is ready, you can transfer this directory to your cluster!
# One way of doing this is by using `rsync`, but you may use other approaches.
# Below, you will find an example of the `rsync` commands you should run once
# you replace the tags with your own credentials.
# First, move one dir back so you are inside `HPC`
cd ../
rsync -avz --copy-links LUCAdup_noleu <uname>@<server>:<path_to_your_wd_in_HPC>
```

----

**IMPORTANT INFORMATION REGARDING THE PAML VERSION INSTALLED ON OUR HPC SERVER**
We have compiled the PAML programs available for the latest version of this software in our HPC server. We have saved the executable files to launch `MCMCtree` and `CODEML` inside our `LUCAdup_noleu` working directory so that they are launched using relative paths to the following executable files: `mcmctree4.10.7` and `codeml4.10.7`.

These are the programs that we will use during all the steps of timetree inference given that the latest PAML version (v4.10.7) has implemented cross-bracing. All inference analyses are therefore run by calling these programs via relative paths. We have decided to work with this system given that we wanted a specific version of PAML not yet available through `conda` and because we use other PAML versions in the HPC server we used.

Alternatively, you can do one of the following:

1. Set up a `conda` environment, install the latest PAML version when available through `conda`, use this `conda` environment to run the subsequent inference analyses. The commands to run `MCMCtree` and `CODEML` are `mcmctree` and `codeml`, respectively.
2. If you do not have other versions of PAML and/or do not want to use `conda`, you can also install the latest PAML version in the HPC server under your user account, export it to your PATH. The commands to run `MCMCtree` and `CODEML` are `mcmctree` and `codeml`, respectively.

----

Now, we need to generate other input files to estimate the Hessian and the gradient: the input control files for `CODEML`. To do this in a reproducible manner, you can use the [script `generate_prepcodeml.sh`](scripts/generate_prepcodeml.sh), which you can find in the [`01_PAML/00_CODEML/scripts`](01_PAML/00_CODEML/scripts) and which you should have just transferred to your HPC. Now, connect to your server and run the next code snippet, where you will execute this script. Specifically, the [`generate_prepcodeml.sh` script](scripts/generate_prepcodeml.sh) needs one argument: the amount of alignment files to be analysed. As we only have one alignment file, we will use `1` as the argument:

```sh
# Run from `LUCAdup_noleu/scripts` in the HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
# In this case, there is only one 
# alignment and hence no need to run the 
# `generate_prepcodeml.sh` in a for` loop
num_aln=1
./generate_prepcodeml.sh $num_aln
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `LUCAdup_noleu/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_codeml/*ctl
grep 'treefile' */prepare_codeml/*ctl
grep 'aaRatefile' */prepare_codeml/*ctl
```

## 3. Run `CODEML`

### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `CODEML` (control file) in our HPC server, we will be manually running `MCMCtree` inside each `prepare_codeml` directory (see file structure above) in a special mode that launches `CODEML` for the sole purpose want: to infer the vectors and matrix required to approximate the likelihood calculation.

```sh
# Run `MCMCtree` from
# `LUCAdup_noleu/Hessian/1/prepare_codeml`
# dir on the HPC. 
# Please change directories until
# you are in there.
# The first command to change directories 
# will work if you are still in 
# `main/Hessian`, otherwise ignore and 
# move to such directory with the command
# that best suits your current directory.
cd 1/prepare_codeml
../../../mcmctree4.10.7 prepcodeml.ctl
```

First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your cluster!):

```text
*** Locus 1 ***
running codeml tmp0001.ctl

AAML in paml version 4.9h, March 2018
ns = 243        ls = 1123
Reading sequences, sequential format..
Reading seq #243: Volvox_carteri_mito_2                          
Sequences read..

1123 site patterns read, 1135 sites
Counting frequencies..

   235224 bytes for distance
   359360 bytes for conP
    35936 bytes for fhK
  5000000 bytes for space
```

As soon as you see the last line, you will see that various `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal tjat you have used to run such command. Once you have done this, you can check that the control file you will later need has been created:

```sh
# Run from the `LUCAdup_noleu/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 1
```

Note that, when we ran the commands above, we were not interested in running `CODEML` or `MCMCtree`. We just wanted to execute `MCMCtree` with option `usedata = 3` so that it generates the `tmp000*` files that `CODEML` will later need to estimate the branch lengths, the gradient, and the Hessian. We do this analysis in two steps given that there are restrictions in the HPC we are using that do not allow us to run `CODEML` + `MCMCtree` in one unique job within a reasonable amount of time. In addition, we want to modify some settings in the control file that is automatically generated when enabling `usedata = 3` so that they match what we want to do for our inference. In a nutshell, this is what you will be doing:

1. Run `MCMCtree` to generate the `tmp000*` files.
2. Modify the `tmp0001.ctl` file according to the settings we want to enable to analyse our dataset with `CODEML`.
3. Run `CODEML` using the `tmp000*` files so that it estimates the branch lengths, the gradient, and the Hessian and saves them in a file called `rst2`.
4. Generate the final `in.BV` file for our dataset, which will be later used by `MCMCtree`.

Once all `tmp000*` files are generated for all alignments, we need to make sure that the following options have been enabled:

1. The (absolute or relative) path to the `lg.dat` file with the protein substitution model should be the argument of option `aaRatefile` in the `tmp000*.ctl`.
2. Following the `Tutorial 4: Approximate likleihood with protein data`, one of the sections in the [`MCMCtree Tutorials` document](http://abacus.gene.ucl.ac.uk/software/MCMCtree.Tutorials.pdf), we set the template control file to use gamma rates among sites instead of the default model, which uses no gamma rates. As we had already defined these options in the template control file, these are already enabled in the `tmp000*.ctl` file generated after running `MCMCtree` (i.e., the control file should already have (i) `fix_alpha = 0` and `alpha = 0.5` to enable the estimation of alpha for the gamma distribution for variable substitution rates across sites by starting the search of the value of $\alpha$ at 0.5 and (ii) `ncatG = 4` with the number of categories for this distribution equal to 4). We can double check to make sure that everything is fine.
3. According to the settings above, we will be using the "LG+F+G4" empirical model for amino acid data, which is enabled by setting `model = 3`, `ncatG = 4`, `aaRatefile = <path_to_lg_matrix>/lg.dat`.
4. In addition, you need to make sure that option `method = 1` is enabled, which will speed up the computation of the Hessian and the gradient.

We can run the next code snippet to very that the four requirements aforementioned are met:

```sh
# Run from the `LUCAdup_noleu/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as datasets you have
grep 'alpha' */*/tmp0001.ctl   # You should see `fix_alpha = 0` and `alpha = 0.5`
grep 'lg.dat' */*/tmp0001.ctl  # You should see the absolute path to the `lg.dat` file in your system
grep 'ncatG' */*/tmp0001.ctl   # You should see `ncatG = 4`
grep 'model' */*/tmp0001.ctl   # You should see `model = 3` (i.e., empirical+F model)
```

### Executing `CODEML`

We can now run `CODEML` given that we have the control file ready as well as all the required input files!

We have created a template bash script with flags (i.e., see script `pipeline_Hessian_CODEML_template.sh` in the [`scripts` directory](01_PAML/00_Hessian/scripts/pipeline_Hessian_CODEML_template.sh)), which will be replaced with the appropriate values by another bash script (`generate_job_CODEML.sh`, also saved in the [`scripts` directory](01_PAML/00_Hessian/scripts/generate_job_CODEML.sh)). Please note that the second bash script will edit the template bash script according to the data alignment/s that will be analysed. We had already transferred these scripts to the HPC server when setting up our file structure. Therefore, we just need to execute the following code snippet there:

```sh
# Run from `LUCAdup_noleu` dir on your HPC. Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=1
# Arg1: Number of alignments
# Arg2: Path to the pipeline directory
# Arg3: Name of the working directory (i.e., `LUCAdup_noleu` in this analysis)
./generate_job_CODEML.sh $num_aln $home_dir/pipelines_Hessian LUCAdup_noleu
```

Next, we will go to the `pipelines_Hessian` directory and run the script that will have been generated using the commands above:

```sh
# Run from `LUCAdup_noleu/pipelines_Hessian` dir on your HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipeline_Hessian.sh`
ll *
# Now, execute this bash script
chmod 775 *sh
qsub pipeline_Hessian.sh
```

Once `CODEML` finishes, we are ready to generate the `in.BV` file that we will later use when running `MCMCtree` to approximate the likelihood calculation:

```sh
# Run from dir `LUCAdup_noleu/Hessian/1` dir on your HPC
# Please change directories until
# you are there. Then, run the following
# commands.
printf "\nGenerating in.BV files for dir 1  ... ...\n\n"
cp rst2 in.BV
```

Next, we can transfer the output generated by `CODEML` to our HPC so that we can keep a backup:

```sh
# Run from `00_CODEML` in your PC
mkdir out_CODEML
cd out_CODEML
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/LUCAdup_noleu/Hessian .
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/LUCAdup_noleu/pipelines_Hessian .
# Remove unnecessary empty output files
rm pipelines_Hessian/*sh.o*
```

We can now proceed to timetree inference with `MCMCtree` using the concatenated dataset while enabling cross-bracing across all possible mirrored nodes, regardless a fossil calibrations is constraining their age. [You can click this link to move to the next `README.md` file](../02_MCMCtree/README.md)!
