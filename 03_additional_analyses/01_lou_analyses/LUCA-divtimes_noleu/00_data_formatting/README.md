# Data formatting

Before proceeding with timetree inference, we need to make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

## Alignment files

### Parse individual gene alignments before concatenation

If you look at [the individual alignment files inside `alignment/partitioned`](../../../../00_data_formatting/00_raw_data/alignment/partitioned/), you will see that the alignments are already in a one-line FASTA format, which makes it easier to parse. First, we will remove all those sequences that are not included in the concatenated alignment previously parsed (i.e., 246 taxa after removing _E. coli_). We will use the previous FASTA alignment that we used in our main analyses, which includes the ATP alignment block, so that we can obtain a list of all the taxa:

```sh
# Run from `00_raw_data/alignment`
grep '>' ../../../../../../00_data_formatting/00_raw_data/alignment/concat5_one_line.fa | sed 's/>//g' > all_taxa.txt
# Remove E.coli!
sed -i '/^Escherichia_coli_2$/d' all_taxa.txt
```

> NOTE: We realised that 5 taxa had a sequence full of gaps: `Chlamydomonas_reinhardtii_chloro_2` (total: 6 copies, including `chloro_2`), `Cyanophora_paradoxa_1` (total: 2 copies, and so we will also remove `Cyanophora_paradoxa_2`).Therefore, we decided to remove these taxa from all the alignments:

```sh
# Remove 3 taxa
sed -i '/^Chlamydomonas_reinhardtii_chloro_2$/d' all_taxa.txt
sed -i '/^Cyanophora_paradoxa_1$/d' all_taxa.txt
sed -i '/^Cyanophora_paradoxa_2$/d' all_taxa.txt
```

Now, we can run [our in-house R script `Parse_partitioned_alignments.R`](scripts/Parse_partitioned_alignments.R) to generate our filtered alignments containing all the taxa present in the concatenated alignment. Note that, whenever a taxon is missing, the filtered individual FASTA alignment will include such taxon and a sequence of gaps as long as the rest of aligned sequences. The output alignments will be generated inside directory [`00_raw_data/alignment/partitioned`](00_raw_data/alignment/partitioned/). Nevertheless, they are still in FASTA format, and so we need to convert them into PHYLIP format:

```sh
# Run from `00_raw_data/alignment/partitioned`
# Uncomment the next line with '##' if you are following this tutorial
# to automatically change to the desired directory
## cd partitioned
count=0
mkdir ind_aln
for i in *fasta
do
num=$( grep '>' $i | wc -l )
len=$( sed -n '2,2p' $i | sed 's/\r//' | sed 's/\n//' | wc -L )
name=$( echo $i | sed 's/\_..*//' )
# Remove weird characters
sed -i 's/\r//g' $i
perl ../../../../../../../src/FASTAtoPHYL.pl $i $num $len
mv log_lenseq.txt log_lenseq_$name".txt"
# Move PHYLIP format to input data
count=$(( count + 1 ))
mv *phy ind_aln/LUCAdup_aln_$name"_p"$count".phy"
done
# Create `01_inp_data` if not yet created
if [ ! -d "../../../01_inp_data" ]
then
mkdir ../../../01_inp_data
fi
mv ind_aln/ ../../../01_inp_data/
```

The individual gene alignments in PHYLIP format will be saved inside `01_inp_data/ind_aln`. You will also find four log files called `log_lenseq*.txt` inside [the `00_raw_data/alignment/partitioned` directory](00_raw_data/alignment/partitioned), one for each alignment that has been converted into PHYLIP format.

The partitioned alignment is now in the correct format, so we can start to generate the concatenated alignment file!

### Concatenation

Now that we have the individual gene alignment in the correct PHYLIP format, we just need to concatenate them in the same order as they were in the previous analyses without the ATP gene alignment: `ATP_filt.fasta`, `EF_filt.fasta`, `SRP_filt.fasta`, and `Tyr_filt.fasta`. We will use the the [`fasta-phylip-partitions` pipeline](https://github.com/sabifo4/fasta-phylip-partitions) to concatenate these alignments. To make things easier, this tool has been compressed and saved in the [`src`](../../src/) directory. You just need to run the next commands to concatenate the dataset:

```sh
# Run from `01_lou_analyses/src` to uncompress the file
# Run these commands only if you have not already done that with 
# other analyses
tar -xvf fasta-phylip-partitions.tar.gz
chmod 775 fasta-phylip-partitions/src/*sh
chmod 775 fasta-phylip-partitions/src/Tools/*
# Now, change to `00_raw_data/00_raw_data/alignment` and run the following
# code to prepare the input data for the pipeline
cd ../LUCA-divtimes_noatp/00_data_formatting/00_raw_data/alignment/
mkdir concatenated
cd concatenated
cp ../all_taxa.txt species_names.txt
cp ../partitioned/*fasta .
# Now, run the pipeline. In essence, the first argument is the current 
# directory ("."), the second argument the tag ID for the job
# ("cytb"), and then a tag specifying if the alignment 
# needs to be partitioned into CPs, which we do now want, and hence use 
# "partN".
# If you do not have this pipeline on your system, you can
# run the code below as it is using the source code that has already
# been added in the `src` directory. For more details about 
# how to use this pipeline, you can find this in the
# following link: https://github.com/sabifo4/fasta-phylip-partitions/blob/main/README.md
# NOTE: If you are running this code in a directory that you have synched to Dropbox or another 
# cloud storage and you have issues, just move the folder out of the synched directory and run the 
# command below again
../../../../../src/fasta-phylip-partitions/src/Run_tasks.sh . LUCAdup partN
# Move now the concatenated alignment to the `01_inp_data` directory
mv phylip_format/02_concatenated_alignments/LUCAdup_concat.aln ../../../01_inp_data/LUCAdup_243sp_aln.phy
# Remove fasta files and text files no longer needed
rm *fasta *txt
```

We have now both partitioned and unpartitioned alignments in the correct format, so we can start to parse the tree file!

## Tree files

We re-ran `IQ-TREE 2` but with the [PHYLIP alignment file we generated following the steps indicated above](01_inp_data/LUCAdup_243sp_aln.phy). In addition, we fixed the re-arranged tree topology we used in our previous analyses but now without the 3 taxa that we had to additional removed (see [`LUCAdup_topo.tree`](00_raw_data/trees/00_IQTREE/LUCAdup_topo.tree)). The command we used to run `IQ-TREE 2` was the following:

```sh
# If ran locally, then uncomment the following line commented with '##' and run the corresponding command.
# Choose `iqtree`, `iqtree2`, or other aliases you may have to run your version of `IQ-TREE` on your machine:
##iqtree2 -s ../../../01_inp_data/LUCAdup_243sp_aln.phy -g LUCAdup_topo.tree -nt AUTO -m LG+C20+F+G -pre BLE_LUCA_DUP_NOLEU
# If not, then adapt the following command to your file structure and any other aliases you may have!
iqtree2 -s LUCAdup_243sp_aln.phy -g LUCAdup_topo.tree -nt AUTO -m LG+C20+F+G -pre BLE_LUCA_DUP_NOLEU
```

You can find the output files and the file with the constrained topology in the [`00_IQTREE` directory](00_raw_data/trees/00_IQTREE). Then, we loaded the output phylogeny by `IQ-TREE 2`, [`BLE_LUCA_DUP_NOLEU.treefile`](00_raw_data/trees/00_IQTREE/BLE_LUCA_DUP_NOLEU.treefile) in `FigTree`. Then, we rooted the tree and saved it in directory [`00_IQTREE`](00_raw_data/trees/00_IQTREE/LUCAdup_topo_bl_rooted.tree).

We will be using the same tree files that we used when the ATP gene alignment was present, but now we have removed 3 additional taxa -- and, consequently, some flags will disappear (i.e., some cross-braced calibrations are no longer part of the tree topology when these taxa are removed). Specifically, the flags that are no longer used to constrain some node ages are the following:

* `#36`: when removing `Chlamydomonas_reinhardtii_chloro_2`. There are 5 flags left.
* `#31`: when removing `Cyanophora_paradoxa_1` and `Cyanophora_paradoxa_2`. There are none flags left.

In that way, we just copied and edited the input files from our previous analyses ("strategy A" only) so that these taxa and the aforementioned calibrations are removed. You can find the final edited files under directory [`01_inp_data`](01_inp_data):

* [`LUCAdup_243sp_allcb_calib_MCMCtree.tree`](01_inp_data/LUCAdup_243sp_allcb_calib_MCMCtree.tree): calibrated tree file under strategy "cross-bracing A".
* [`LUCAdup_243sp_uncalib.tree`](01_inp_data/LUCAdup_243sp_uncalib.tree): uncalibrated tree used to run `CODEML`.

---

Now, we have all the input files correctly formatted to be analysed with PAML programs! We will subsequently run `CODEML` to calculate the branch lengths, the gradient, and the Hessian for the concatenated alignment so that `MCMCtree` can use them to approximate the likelihood calculation during timetree inference! You can move to the [`01_PAML` directory](../01_PAML) now and access directories [00_CODEML](../01_PAML/00_CODEML/README.md) to follow the step-by-step guidelines for the next analyses!
