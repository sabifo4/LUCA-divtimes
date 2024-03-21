# Data formatting

Before proceeding with timetree inference, we need to make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

## Input files

For this additional analysis, we will be analysing each gene individually. As we have already formatted these gene alignments, we just need to move them to a file structure similar to the one we have been using in our previous analyses:

```sh
# Run from `00_HPC_dat`
# Create file structure
mkdir LUCAdup_indgenes
cd LUCAdup_indgenes
num_aln=5
for i in `seq 1 $num_aln`
do
mkdir -p alignments/$i
mkdir -p Hessian/$i/
mkdir -p trees/calibrated
mkdir control_files
done
mkdir scripts
# Now, populate `Hessian` with the `in.BV` files that were
# previously generated for each individula gene alignment,
# which will be saved under `alignments` 
for i in `seq 1 5`
do
cp ../../../../00_data_formatting/01_inp_data/ind_aln/*p$i".phy" alignments/$i
cp ../../../../01_PAML/01_CODEML/out_CODEML/Hessian_part/$i/rst2 Hessian/$i/in_p$i".BV"
done
# We will reuse the calibrated tree file under strategy A 
# (i.e., (all nodes with same speciation events are mirrored) 
cp ../../../../00_data_formatting/01_inp_data/LUCAdup_246sp_allcb_calib_MCMCtree.tree trees/calibrated
# Copy template control file to run MCMCtree and `lg.dat` from
# previous analyses
cp ../../../../01_PAML/00_CODEML/control_files/*MCMCtree*sh control_files/
```

The next step will consist of sending this file structure with the input data to the cluster, where we will proceed to run our timetree inference analyses:

```sh
# Run from `00_HPC_data/`
rsync -avz --copy-links LUCAdup_indgenes <uname>@<login_details>:<path_to_main_dir>
```

Now, we have everything we need to run `MCMCtree` for timetree inference using the approximate likelihood calculation! You can move to the [`01_PAML` directory](../01_PAML) now and access directory [00_MCMCtree](../01_PAML/00_MCMCtree/README.md) to follow the step-by-step guidelines for the next analyses!
