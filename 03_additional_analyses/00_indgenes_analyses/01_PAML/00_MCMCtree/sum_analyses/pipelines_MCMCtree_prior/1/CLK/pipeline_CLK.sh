#!/bin/bash
#$ -cwd                    # Run the code from the current directory
#$ -V                      # Export environment to job
#$ -j y                    # Merge the standard output and standard error
#$ -l h_rt=240:00:00       # Limit each task to 10 days
#$ -l h_vmem=1G            # Requested RAM
#$ -t 1-6                

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# --------------------------------------- #
# Creating file structure to run MCMCtree_prior #
# --------------------------------------- # 

# 1. First arg is dataset
dir=$( echo 1 )       # 1, 2, 3... According to the dataset being analysed

# 2. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/\/LUCAdup_indgenes\/..*/\/LUCAdup_indgenes\//' )
cd $main_dir/Hessian/$dir
codeml_dir=$( pwd )
name_inBV=`ls *BV`
path_inBV=$( echo $codeml_dir/$name_inBV )
cd $main_dir/MCMCtree_prior/$SGE_TASK_ID/CLK/$dir
home_dir=$( pwd )
cd $main_dir/control_files
ctl_dir=$( pwd )
path_lgdat=$( echo $ctl_dir/lg.dat | sed 's/\//\\\//g' | sed 's/\_/\\\_/g' | sed 's/\./\\\./g' )
cd $main_dir/alignments/$dir 
aln_dir=$( pwd )
name_aln=`ls *phy`
path_aln=$( echo $aln_dir/$name_aln | sed 's/\_/\\\_/g' | sed 's/\//\\\//g' | sed 's/\./\\\./g' )
cd $main_dir/trees/calibrated
tree_dir=$( pwd )
name_tree=`ls *tree`
path_tree=$( echo $tree_dir/$name_tree | sed 's/\_/\\\_/g' | sed 's/\//\\\//g' | sed 's/\./\\\./g' )

# 3. Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.MCMCtree_prior.dir$dir"_r"$SGE_TASK_ID".txt"; done;) 2>&1
start=`date`
echo Job starts":" $start

# 4. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Copy control file
cp $ctl_dir/prepcodeml.ctl $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
sed -i 's/ALN/'${path_aln}'/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
sed -i 's/TREE/'${path_tree}'/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
sed -i 's/usedata..*/usedata = 0/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
sed -i 's/model..*/model = 3/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
sed -i 's/\ ndata..*/\ ndata\ \=\ 1/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
sed -i 's/lg\.dat/'${path_lgdat}'/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
if [[ GBM =~ "GBM" ]]
then 
sed -i 's/clock..*/clock\ \=\ 3/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
elif [[ GBM =~ "ILN" ]]
then 
sed -i 's/clock..*/clock\ \=\ 2/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
fi
sed -i 's/clock..*/clock = 1/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl" 

# Soft link the in.BV file here
#ln -s $path_inBV $home_dir/in.BV

# 5. Enable duplication setting for inequality constrained if required
# and run MCMCtree_prior
duplication=1
if [[ $duplication -eq 1 ]]
then 
printf "\nduplication = 1\n" >> $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"
printf "\nRunning MCMCtree_prior for divergence times estimation ... ...\n"
cd $home_dir
$main_dir/mcmctree4.10.7 *.ctl
else
printf "\nRunning MCMCtree_prior for divergence times estimation ... ...\n"
cd $home_dir
$main_dir/mcmctree4.10.7 *.ctl
fi
sed -i 's/clock..*/clock = 1/' $home_dir/mcmctree_$dir"_r"$SGE_TASK_ID".ctl"

# 6. Close
printf "\n"
echo MCMCtree_prior FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end

# 7. Return to pipedir
cd $pipeline_dir

