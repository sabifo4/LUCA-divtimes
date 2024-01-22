# `CODEML` analysis - partitioned dataset

## 1. Pick rate prior

We will use the same vague gamma distribution that we used when analysing the concatenated alignment: `G(2,2.5)`. For more details on how the `beta` value was calculate, you can take a look at [the R in-house script we previously ran](../00_CODEML/scripts/calculate_rateprior.R) and at the section `Pick rate prior` in the corresponding [`00_CODEML/README.md` file](../06_CODEML/README.md#1-pick-rate-prior).

## 2. Set up the file structure

We will now add the following directories to the already pre-defined file structure under our working directory `LUCAdup`:

```text
LUCAdup/
  |
  |- alignments_part/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments
  |       
  |- Hessian_part/
  |    |- X # Directory for alignment X -- have as many directories as alignments
  |          
  |- pipelines_Hessian_part # Directory where the pipeline to run `CODEML` will be executed

```

To create the new directories to be saved under the `LUCAdup` file structure, we will run the following commands from our PC before transferring them to the HPC:

```sh
# Run the following commands from 
# `01_CODEML`
mkdir -p HPC/LUCAdup
cd HPC/LUCAdup 
i=5
for i in `seq 1 $i`
do
mkdir -p alignments_part/$i
mkdir -p Hessian_part/$i/prepare_codeml
done
mkdir pipelines_Hessian_part
mkdir scripts_part
```

Once the file structure is created, we can now populate it with the partitioned alignment file we generated a while ago and new in-house scripts for this analysis with a partitioned dataset:

```sh
# Run from `HPC/LUCAdup`
# Copy alignment
i=5
for i in `seq 1 $i`
do
cp ../../../../00_data_formatting/01_inp_data/ind_aln/*p$i".phy" alignments_part/$i
done
# Last, copy the in-house bash scripts with our pipeline
cp ../../scripts/*sh scripts_part/
# Once everything is ready, you can transfer these two directories to your cluster!
# One way of doing this is by using `rsync`, but you may use other approaches.
# Below, you will find an example of the `rsync` commands you should run once
# you replace the tags with your own credentials.
# First, move one dir back so you are inside `HPC`
rsync -avz --copy-links *part <uname>@<server>:<path_to_your_wd_in_HPC>/LUCAdup
```

Now, we need to generate the input control files for `CODEML`. To do this in a reproducible manner, you can use our in-house bash scripts that you will in the [`scripts` directory](01_PAML/00_CODEML/scripts), which you should have just transferred to your HPC. Now, connect to your server and run the next code snippet:

```sh
# Run from `LUCAdup/scripts_part` in the HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
# In this case, there is only one 
# alignment and hence no need to run the 
# `generate_prepcodeml.sh` in a for` loop
num_aln=5
for i in `seq 1 $num_aln`
do
./generate_prepcodeml.sh $i
done
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `LUCAdup/Hessian_part` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_codeml/*ctl
grep 'treefile' */prepare_codeml/*ctl
grep 'aaRatefile' */prepare_codeml/*ctl
```

## 3. Run `CODEML`

### Preparing input files

We will now be manually running `MCMCtree` inside each `prepare_codeml` directory (see file structure above) as we did with the concatenated dataset:

```sh
# Run `MCMCtree` from
# `LUCAdup/Hessian_part/[1-5]/prepare_codeml`
# dirs on the HPC. 
# Please change directories until
# you are in there.
#
# Alternatively, you can run the following
# `for` loop from the `Hessian_part` directory
# and keep pressing "ctrl+C" every
# time that you see the `tmp*` files being
# generated to move on to the next iteration
base_dir=$( pwd )
for i in `seq 1 5`
do
cd $base_dir/$i/prepare_codeml
## IMPORTANT TO USE 4.10.7 WHEN PARTIITONED DATA!
../../../mcmctree4.10.7 prepcodeml.ctl # you will need to kill these processes as soon as you have the files you need!
done
```

First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your PC!):

```text
*** Locus 1 ***
running codeml tmp0001.ctl

AAML in paml version 4.9h, March 2018
ns = 246        ls = 161
Reading sequences, sequential format..
Reading seq #246: Volvox_carteri_mito_2         5_2
Sequences read..

161 site patterns read, 161 sites
Counting frequencies..

   241080 bytes for distance
    51520 bytes for conP
     5152 bytes for fhK
  5000000 bytes for space
```

As we did with the concatenated dataset, as you see the last line, you will see that the `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal where you have run such command. Once everything is done, you can check that the control files you will later need for each individual alignment have been created:

```sh
# Run from the `LUCAdup/Hessian_part` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 5
```

As we did with out concatenated dataset, we need to make sure that the settings in the control files correspond to those we want to enable to analyse our dataset:

```sh
# Run from the `LUCAdup/Hessian_part` dir on your local
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

We are now ready to run `CODEML`!

We will run our in-house bash scripts to prepare the bash script that we will submit to our HPC running the following code snippet:

```sh
# Run from `LUCAdup` dir on your HPC. Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
cd scripts_part
chmod 775 *sh
num_aln=5
# Arg1: Number of alignments
# Arg2: Path to the pipeline directory
# Arg3: Name of the working directory (i.e., `LUCAdup` in this analysis)
for i in `seq 1 $num_aln`
do
./generate_job_CODEML.sh $i $home_dir/pipelines_Hessian_part LUCAdup
done
```

Now, we just need to go to the `pipelines_Hessian_part` directory and run the script that will have been generated using the commands above:

```sh
# Run from `LUCAdup/pipelines_Hessian_part` dir on your HPC.
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

Once `CODEML` finishes, we are ready to generate the `in.BV` files that we will later use when running `MCMCtree` using the approximate likelihood calculation:

```sh
# Run from dir `LUCAdup/Hessian_part` dir on your HPC
# Please change directories until
# you are there. Then, run the following
# commands.
num_aln=5
for i in `seq 1 $num_aln`
do
printf "\nGenerating in.BV files for dir $i  ... ...\n\n"
if [[ $i -eq 1 ]]
then
cat $i/rst2 > in_5parts.BV
printf "\n" >> in_5parts.BV
else
cat $i/rst2 >> in_5parts.BV
printf "\n" >> in_5parts.BV
fi
done
```

Next, we can transfer the output generated by `CODEML`:

```sh
# Run from `01_CODEML` in your PC
mkdir out_CODEML
cd out_CODEML
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/LUCAdup/Hessian_part .
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/LUCAdup/pipelines_Hessian_part .
# Delete blank output files
rm pipelines_Hessian_part/*sh.o*
```

We can now proceed to timetree inference with `MCMCtree` using the partitioned dataset while enabling cross-bracing across all possible mirrored nodes, regardless a fossil calibrations is constraining their age. [You can click this link to move to the next `README.md` file](../05_MCMCtree_part/README.md)!
