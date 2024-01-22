# Data formatting

Before proceeding with timetree inference, we need to make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

## Alignment files

### Unpartitioned

If you open [the unpartitioned alignment file](00_raw_data/alignment/concat5.fa), you will see that each aligned sequence is not in a unique line, which makes it harder to parse the file to convert it into PHYLIP format.

We will first generate a one-line FASTA file with [our in-house bash script](scripts/one_line_fasta.pl). Then, we will use the output file as input to [our second in-house bash script](scripts/FASTAtoPHYL.pl) to obtain a protein alignment in PHYLIP format:

```sh
# Run from `00_raw_data/alignment`
# First, convert into a one-line FASTA file
perl ../../scripts/one_line_fasta.pl concat5.fas 
# Now, format to PHYLIP
num=$( grep '>' *one_line.fa | wc -l )
len=$( sed -n '2,2p' *one_line.fa | sed 's/\r//' | sed 's/\n//' | wc -L )
perl ../../scripts/FASTAtoPHYL.pl *one_line.fa $num $len
# Move PHYLIP format to input data
mkdir ../../01_inp_data
mv concat5_one_line.phy ../../01_inp_data/LUCAdup_aln.phy
```

You will now see a new directory called `01_inp_data` inside directory [`00_data_formatting`](README.md). If you navigate to this newly created `01_inp_data` directory, you will find the alignment in PHYLIP format (i.e., the input file we need!). You will also find a log file called `log_lenseq.txt` inside [the `00_raw_data/alignment` directory](00_raw_data/alignment) where you can read how many taxa were parsed and the length of the sequence. Nevertheless, given that we modified the resulting inferred tree topology to implement the available fossil calibrations (see above) and following consensus views of species-level relationships, we decided to remove taxon _E. coli_. In that way, we also need to remove the corresponding sequence from the alignments:

```sh
# Run from `01_inp_dat`
# Uncomment the next line with '##' if you are following this tutorial
# to automatically change to the desired directory
## cd ../../01_inp_dat
cp LUCAdup_aln.phy LUCAdup_246sp_aln.phy
#cp LUCAdup_5parts_aln.phy LUCAdup_246sp_5parts_aln.phy
sed -i '/^Escherichia_coli_2..*$/d' LUCAdup_246sp_aln.phy
sed -i 's/247   /246   /g' LUCAdup_246sp_aln.phy
```

The unpartitioned alignment is now in the correct format, so we can start to generate the partitioned alignment file!

### Partitioned

If you look at [the individual alignment files inside `alignment/partitioned`](00_raw_data/alignment/partitioned), you will see that the alignments are already in a one-line FASTA format, which makes it easier to parse. First, we will remove all those sequences that are not included in the concatenated alignment previously parsed (i.e., 246 taxa after removing _E. coli_). We will first generate a list of taxa that need to be included:

```sh
# Run from `00_raw_data/alignment
# Uncomment the next line with '##' if you are following this tutorial
# to automatically change to the desired directory
## cd ../00_raw_data/alignment
grep '>' concat5_one_line.fa | sed 's/>//g' > all_taxa.txt
# Remove E.coli!
sed -i '/^Escherichia_coli_2$/d' all_taxa.txt
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
perl ../../../scripts/FASTAtoPHYL.pl $i $num $len
mv log_lenseq.txt log_lenseq_$name".txt"
# Move PHYLIP format to input data
count=$(( count + 1 ))
mv *phy ind_aln/LUCAdup_aln_$name"_p"$count".phy"
done
mv ind_aln/ ../../../01_inp_data/
# Generate partitioned alignment
cd ../../../01_inp_data/
count=0
for i in ind_aln/*
do
count=$(( count + 1 ))
if [[ $count -eq 1 ]]
then
cat $i > LUCAdup_246sp_5parts_aln.phy
printf "\n" >> LUCAdup_246sp_5parts_aln.phy
else
cat $i >> LUCAdup_246sp_5parts_aln.phy
printf "\n" >> LUCAdup_246sp_5parts_aln.phy
fi
done
```

The final partitioned alignment, `LUCAdup_246sp_5parts_aln.phy`, will be saved inside directory [`01_inp_data`](01_inp_data), while the individual gene alignments in PHYLIP format will be saved inside `01_inp_data/ind`. You will also find five log files called `log_lenseq*.txt` inside [the `00_raw_data/alignment/partitioned` directory](00_raw_data/alignment/partitioned), one for each alignment that has been converted into PHYLIP format.

We have now both partitioned and unpartitioned alignments in the correct format, so we can start to parse the tree file!

## Tree files

### Phylogeny inference

#### IQ-TREE 2 (round 1)

We ran `IQ-TREE 2` with our [main alignment file in FASTA format](00_raw_data/alignment/concat5.fas), the topological constraints specified in a file called [`Minimal_constraint_pLUCA.tre`](00_raw_data/trees/01_round2/IQTREE/Minimal_constraint_pLUCA.tre), and under the LG+C20+F+G protein model. The command used was the following:

```sh
# Note that this command was executed in an HPC server, and thus all the files were saved
# in a unique directory (no relative paths related to the file structure in this repository)
# Note that the command `iqtree` executes the following IQ-TREE version:
# >> IQ-TREE multicore version 2.1.4-beta COVID-edition for Linux 64-bit built Jun 24 2021
iqtree -s concat5.fas -g Minimal_constraint_pLUCA.tre -nt 11 -m LG+C20+F+G -pre BLE_LUCA_DUP_CONCAT5
```

You can find the output files and the file with the constrained topology in the [`00_IQTREE/round1` directory](00_raw_data/trees/00_IQTREE/round1/). Once we analysed the inferred best-scoring maximum-likelihood tree (i.e., see [`BLE_LUCA_DUP_CONCAT5.treefile`](00_raw_data/trees/00_IQTREE/round1/BLE_LUCA_DUP_CONCAT5.treefile)), we modified the resulting tree topology so that we could implement the available fossil calibrations (see later in the tutorial) while following consensus views of species-level relationships. In addition, as aforementioned, we decided to remove taxon _E. coli_.

#### IQ-TREE 2 (round 2)

We re-ran `IQ-TREE 2` but, this time, with the [PHYLIP alignment file we generated following the steps indicated above](01_inp_data/LUCAdup_246sp_aln.phy). In addition, we fixed the re-arranged tree topology specified in file [`LUCAdup_allcb_topo.tree`](00_raw_data/trees/02_round3/IQTREE/LUCAdup_allcb_topo.tree). The command we used to run `IQ-TREE 2` was the following:

```sh
# If ran locally, then uncomment the following line commented with '##' and run the corresponding command.
# Choose `iqtree`, `iqtree2`, or other aliases you may have to run your version of `IQ-TREE` on your machine:
##iqtree2 -s ../../../../01_inp_data/LUCAdup_246sp_aln.phy -g LUCAdup_topo.tree -nt AUTO -m LG+C20+F+G -pre BLE_LUCA_DUP_NEWTOP
# If not, then adapt the following command to your file structure and any other aliases you may have!
iqtree -s LUCAdup_246sp_aln.phy -g LUCAdup_topo.tree -nt AUTO -m LG+C20+F+G -pre BLE_LUCA_DUP_NEWTOPO
```

You can find the output files and the file with the constrained topology in the [`00_IQTREE/round2` directory](00_raw_data/trees/00_IQTREE/round2).

### Tree calibration

#### Cross-bracing all duplicated nodes that match the same speciation event

First, we loaded the output phylogeny by `IQ-TREE 2` (round 2, [`BLE_LUCA_DUP_NEWTOPO.treefile`](00_raw_data/trees/00_IQTREE/round2/BLE_LUCA_DUP_NEWTOPO.treefile)) in `FigTree`. Then, we rooted tree and saved it in directory [`00_IQTREE/round2`](00_raw_data/trees/00_IQTREE/round2/LUCAdup_topo_bl_rooted.tree).

Then, we manually incorporated the labels to identify those nodes which ages were going to be constrained (i.e., both using traditional node calibrations as well as cross-braced calibrations) using the `Annotate` feature in `FigTree`. In other words, the various mirrored nodes correspond to the same speciation event and fall into one of these two categories:

* **Equality calibrations with fossil calibrations**: mirrored nodes which age can be constrained based on geological evidence or the fossil record. There is a minimum/maximum age constraint that can be applied to these nodes that correspond to the same speciation event and appear in the phylogeny as duplication events. In total, we modelled 13 uniform densities with either hard and soft bounds based on the geological events and the fossil record to constrain the age of mirrored nodes corresponding to the following clades: LUCA (2 mirrored nodes), TG-OXYPHOTOBACTERIA (2 duplication events), LECA (2 duplication events), CG-OXYPHOTOBACTERIA (2 duplication events), TG-EUKARYA-MITO (2 duplication events), TG-EUKARYA-ARCH (2 duplication events), ARCHAEPLASTIDA (6 duplication events), EMBRYOPHYTA (6 duplication events), EUDICOT-MONOCOT (5 duplication events), CG-FORAMINIFERA (3 duplication events), FUNGI (4 duplication events), EUMETAZOA (4 duplication events), METAZOA (4 duplication events). We used these names as flags/labels to identify the mirrored nodes using `FigTree`.
* **Equality calibrations without fossil calibrations**: mirrored nodes which age can be constrained because we know that such nodes correspond to the same speciation event. Therefore, they appear as duplication events in this phylogeny but we want to constrain their posterior time densities to be the same. We identified a total of 64 speciation events that had mirrored nodes across the phylogeny. We used either the name of the specific clade or random names for the flags/labels to identify the mirrored nodes using `FigTree`.

Under both scenarios, the identified mirrored nodes for the same speciation event will have the same posterior time densities. In other words, the sampling during the MCMC will not be independent for those mirrored nodes that are cross-braced, and therefore the same posterior time densities will be assigned to such nodes given that they represent the same speciation event.

Once we finished to annotate our tree topology, we exported the annotated file in NEXUS format and saved it as [`LUCAdup_allcb_outFigTree`](00_raw_data/trees/01_calibrations/LUCAdup_allcb_outFigTree). The next step consists of extract the tree in Newick format with labels within this NEXUS file, which can be easily parsed in subsequent steps:

```sh
# Run from `00_raw_data/trees/01_calibrations`
grep 'tree tree_1' LUCAdup_allcb_outFigTree > LUCAdup_allcb_calibnames.tree
sed -i 's/..*\[\&R\] //' LUCAdup_allcb_calibnames.tree
sed -i 's/\:[0-9]*\.[0-9]*//g' LUCAdup_allcb_calibnames.tree
sed -i 's/&label=//g' LUCAdup_allcb_calibnames.tree
sed -i 's/\"//g' LUCAdup_allcb_calibnames.tree
# Delete quotation marks!
sed -i "s/'//g" LUCAdup_allcb_calibnames.tree
# Add PHYLIP header -- we have now 246 taxa!
sed -i '1s/^/246 1\n/' LUCAdup_allcb_calibnames.tree
```

We have created a [mapping file](scripts/Calib_converter_allcb.txt) in which the flags/labels included in the tree have been assigned the corresponding calibrations in `MCMCtree` format. If you run [our R in-house script `Include_calibrations_allcb.R`](scripts/Include_calibrations_allcb.R), you will generate a tree file inside `01_inp_data` called `LUCAdup_allcb_calib_MCMCtree.tree`, where now the node labels will have the corresponding fossil calibration in `MCMCtree` format. You will also generate a tree file that can be visualised with graphical viewers such as `FigTree` called `LUCAdup_allcb_outR_fordisplay_calib_MCMCtree.tree`. As soon as we obtained this last output file, we visualised it in `FigTree` and saved the project as `LUCAdup_allcb_calibs_visFigTree`. You can access [this `FigTree` project](00_raw_data/trees/02_round3/LUCAdup_allcb_calibs_visFigTree) in case you want to colour specific branches, modify the labels, generate other output files, etc. In particular, we used this project to output a PDF file with which we could easily visualise the calibrated nodes: `LUCAdup_allcb_calibs.pdf`.

#### Cross-bracing only those duplicated nodes that match the same speciation event and for which there is a fossil calibration

Next, we generated a copy of the file `LUCAdup_allcb_calibnames.tree` (i.e., `cp LUCAdup_allcb_calibnames.tree LUCAdup_allnoncb_calibnames.tree`) and manually removed the flags/labels that identified those nodes that we had previously cross-braced for which there is no fossil constraints. The resulting file is [`LUCAdup_allnoncb_calibnames.tree`](00_raw_data/trees/02_round3/LUCAdup_allnoncb_calibnames.tree). In that way, we will only cross-brace those nodes for which fossil calibrations are used to constrain such node ages.

The corresponding mapping file can be found in the [`scripts` directory](scripts/Calib_converted_fosscb.txt). We then ran [our R in-house script](scripts/Include_calibrations_fosscb.R) to generate the output calibrated file in the [`01_inp_data` directory](01_inp_data): `LUCAdup_246sp_fosscb_calib_MCMCtree`. Note that this output file can be visualised directly in `FigTree`. Due to the way that `FigTree` reads numbers (and being this the format to identify mirrored nodes), only one of the various mirrored nodes will have the fossil information in `MCMCtree` format. The rest of the mirrored nodes, will have a `javascript` tag instead, but it has nothing to do with what `MCMCtree` will be using -- just a formatting issue. For visual purposes, you can always refer back to the `LUCAdup_allcb_outR_fordisplay_calib_MCMCtree` file previously generated.

#### Without cross-bracing

In addition, we used an updated mapping file in which the mirrored nodes for which fossil calibrations are available are not cross-braced. Instead, the traditional `MCMCtree` format is used to constrain them with the same fossil calibrations. You can find the [mapping file inside the `scripts` directory](scripts/Calib_converter_allnoncb.txt).

We then ran [our R in-house script](scripts/Include_calibrations_allnoncb.R) to generate the output calibrated file in the [`01_inp_data` directory](01_inp_data): `LUCAdup_246sp_allnoncb_calib_MCMCtree`. Note that this output file can be visualised directly in `FigTree`.

#### Uncalibrated tree

Last, we generated a file with un uncalibrated tree (i.e., Newick format with only the tree topology) as the input for `CODEML` to infer the branch lengths, gradient, and Hessian required to approximate the likelihood calculation in `MCMCtree` ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)):

```sh
# Run from `01_inp_data`
cp ../00_raw_data/trees/01_calibrations/LUCAdup_allnoncb_calibnames.tree LUCAdup_246sp_uncalib.tree
sed -i 's/\[[A-Z]*-[A-Z]*\]//g' LUCAdup_246sp_uncalib.tree
sed -i 's/\[[A-Z]*\]//g' LUCAdup_246sp_uncalib.tree
sed -i 's/\[[A-Z]*-[A-Z]*-[A-Z]*\]//g' LUCAdup_246sp_uncalib.tree
sed -i 's/\[[A-Z]*-[A-Z]*-[A-Z]*-[A-Z]*\]//g' LUCAdup_246sp_uncalib.tree
```

---

Now, we have all the input files correctly formatted to be analysed with PAML programs! We will subsequently run `CODEML` to calculate the branch lengths, the gradient, and the Hessian for both the concatenated and partitioned alignments so that `MCMCtree` can use them to approximate the likelihood calculation during timetree inference! You can move to the [`01_PAML` directory](../01_PAML) now and access directories [00_CODEML](../01_PAML/00_CODEML/README.md) and [01_CODEML](../01_PAML/01_CODEML/README.md) to follow the step-by-step guidelines for the next analyses!
