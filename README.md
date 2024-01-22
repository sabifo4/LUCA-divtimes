# The nature of the Last Universal Common Ancestor and its impact on the early Earth system -- timetree inference analyses

In this repository, you will find a very detailed tutorial with all the steps you need to follow from data formatting to timetree inference to reproduce the results that we obtained as part of our study: **The nature of the Last Universal Common Ancestor and its impact on the early Earth system**.

To get started, you can clone this repository on your PC and follow all the guidelines given in the various `README.md` files that you shall find inside each directory to follow all the steps we carried out for timetree inference. The summary of the workflow we followed to carry out our analyses is the following:

* Parsing and formatting the input data required to run `PAML` software: sequence files (concatenated and partitioned), calibrated tree files (three calibration strategies, more details below), control files.
* Inferring the mean evolutionary rate to specify a sensible rate prior.
* Running `PAML` software for timetree inference:
  * Using various in-house pipelines to set up the working environment, the file structure, and the control files required to run `PAML` software.
  * Running `CODEML` to calculate the branch lengths, the gradient, and the Hessian, which will be required by `MCMCtree` to enable the approximate likelihood calculation.
  * Running `MCMCtree` with the approximate likelihood calculation enabled for timetree inference. Our main analysis consisted of inferring the timetree for the **partitioned dataset** under both the Geometric Brownian motion (**GBM**, [Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal (**ILN**, [Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) **relaxed-clock models**, where nodes that correspond to the same divergences are cross-braced (i.e., hereby referred to as calibration strategy “cross-bracing A”). In addition, we ran 10 additional inference analyses to benchmark the effect that partitioning, cross-bracing, and relaxed-clock models can have on species divergence time estimation:
    * **GBM + concatenated alignment + cross-bracing A**.
    * **GBM + concatenated alignment + cross-bracing B**, (only nodes that correspond to the same divergences for which there are fossil constraints are cross-braced.
    * **GBM + concatenated alignment + without cross-bracing**.
    * **GBM + partitioned alignment + cross-bracing B**.
    * **GBM + partitioned alignment + without cross-bracing**.
    * **ILN + concatenated alignment + cross-bracing A**.
    * **ILN + concatenated alignment + cross-bracing B**.
    * **ILN + concatenated alignment + without cross-bracing**.
    * **ILN + partitioned alignment + cross-bracing B**.
    * **ILN + partitioned alignment + without cross-bracing**.
* Running MCMC diagnostics for all the chains under each analysis.
* General discussion.

To make it easier for you to navigate this GitHub repository, below you can find a summary of the content you shall find inside each directory.

## [Data formatting](00_data_formatting/README.md)

* **Molecular alignment**: we analysed a molecular alignment consisting of 5 pre-LUCA paralogues (i.e., F- and V- type subunits from ATPases; Elongation Factor Tu and G; Signal Recognition Protein and Signal Recognition Particle Receptor; Tyrosyl-tRNA and Tryptophanyl-tRNA synthetases; and Leucyl- and Valyl- tRNA synthetases) for 246 taxa. You can read the `README.md` file under [`00_data_formatting`](00_data_formatting/README.md) to reproduce all the steps that we carried out to parse the raw data.
  * [`Concatenated molecular alignment (raw)`](00_data_formatting/00_raw_data/alignment/concat5.fas): the raw molecular alignment consisted of 1,080 characters and 247 taxa. Gene families were identified using `BLASTp` ([Altschul 1990](https://pubmed.ncbi.nlm.nih.gov/2231712/)) and downloaded from the NCBI ([Sayers et al., 2010](https://pubmed.ncbi.nlm.nih.gov/21097890/)). Amino acid sequences were aligned with MUSCLE ([Edgar 2004](https://academic.oup.com/nar/article/32/5/1792/2380623)) and trimmed with `TrimAl` ([Capella-Gutiérrez 2009](https://academic.oup.com/bioinformatics/article/25/15/1972/213148)) (`-strict``). Individual gene trees were inferred under the LG+C20+F+G substitution model implemented in IQ-TREE 2 ([Minh et al. 2020](https://academic.oup.com/mbe/article/37/5/1530/5721363?login=false)). These trees were manually inspected and curated to remove non-homologous sequences or recent horizontal gene transfers. According to the filters applied, the resulting gene alignments included 246 species, and multiple copies of each taxa were represented at least twice (for some eukaryotes, they may be represented by plastid, mitochondrial, and nuclear sequences).
  * [`Individual gene alignments for each paralog (raw)`](00_data_formatting/00_raw_data/alignment/partitioned/): you will find the five gene alignments (i.e., one for each paralog) under the [`partitioned` directory](00_data_formatting/00_raw_data/alignment/partitioned/).
  * [`Input concatenated alignment`](00_data_formatting/01_inp_data/LUCAdup_246sp_aln.phy): after parsing the raw concatenated alignment aforementioned, this is the input concatenated alignment that we used for all timetree inference analyses.
  * [`Input partiitoned alignment`](00_data_formatting/01_inp_data/LUCAdup_246sp_5parts_aln.phy): after parsing the raw individual gene alignments aforementioned, this is the input partitioned alignment (i.e., five alignment blocks, one for each gene alignment) that we used for all timetree inference analyses.
* **Phylogeny**: given that a fixed tree topology is required for timetree inference with MCMCtree, we decided to infer the best-scoring ML tree with IQ-TREE 2 under the LG+C20+F+G4 model ([Hoang 2018](https://pubmed.ncbi.nlm.nih.gov/29077904/)) following our previous phylogenetic analyses. We then modified the resulting inferred tree topology to implement the available fossil calibrations (see below) and following consensus views of species-level relationships ([Coleman et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33958449/), [Aouad et al., 202](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-021-01952-0), [Burki et al., 2019](https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347(19)30257-5)).
* **Calibrations**: we used previously established dates for four fossil calibrations and included nine additional calibrations with justifications provided in the supplementary material of the manuscript. Note that five of the newly defined calibrations (i.e., total group Eukarya for the archaeal lineage, Embryophyta, Archaeplastida, and Eumetazoa) have been revised to follow the latest geochronological updates and literature as of December 2023, and so differ from those used in previous studies1–3. Note that the root of a phylogeny requires constraining in molecular clock-dating analyses. We constrained the root (representing the duplication of pre-LUCA paralogues) with a hard upper bound equal to the age of the moon forming impact. We generated 3 control files that we used as input fileS for our R in-house scripts to semi-automatically calibrate the tree topologies with the `MCMCtree` notation:
  * ["Cross-bracing A": `Calib_converter_allcb.txt`](00_data_formatting/scripts/Calib_converter_allcb.txt): text file that is used a converter file by our [R in-house script `Include_calibrations_allcb.R`](00_data_formatting/scripts/Include_calibrations_allcb.R) to calibrate the fixed tree topology aforementioned using `MCMCtree` notation. Note that the calibration format used in this file does not follow the conventional `MCMCtree` notation. Specifically, the format you shall find is used to cross-brace those nodes with the same speciation event, and so the posterior density for such "mirrored" nodes will be the same.
  * ["Without cross-bracing": `Calib_converter_allnoncb.txt`](00_data_formatting/scripts/Calib_converter_allnoncb.txt): text file that is used a converter file by our [R in-house script `Include_calibrations_allnoncb.R`](00_data_formatting/scripts/Include_calibrations_allnoncb.R) to calibrate the fixed tree topology aforementioned using `MCMCtree` notation. Note that the format used in this file is the conventional `MCMCtree` format in which the age of each node is constrained by a specific probability density. In that way, cross-bracing is not enabled and each node shall have the corresponding posterior density (i.e., no mirrored nodes).
  * ["Cross-bracing B": `Calib_converter_fosscb.txt`](00_data_formatting/scripts/Calib_converter_fosscb.txt): text file that is used a converter file by our [R in-house script `Include_calibrations_fosscb.R`](00_data_formatting/scripts/Include_calibrations_fosscb.R) to calibrate the fixed tree topology aforementioned using `MCMCtree` notation. There is one difference between this calibration strategy and "cross-bracing A": only those nodes for which there is a fossil calibration have been cross-braced.

## [PAML analyses](01_PAML)

Note that all the directories that you shall find inside the [`01_PAML` directory](01_PAML) have their own `README.md` file. In that way, you can navigate to each of these directories to read more about how we ran `CODEML` and `MCMCtree` when using either the concatenated or the partitioned alignments and when running the inference analyses under different settings (e.g., data type, calibration strategy, relaxed-clock model).

Note that each directory has all the input data, intermediate files, and scripts that you will need to follow the step-by-step guidelines if you clone this repository on your PC. Nevertheless, note that the `MCMCtree` output files with all the collected samples during the MCMC that we use to summarise our results are too large to be stored here. When you want to run the tutorials under the `01_PAML/*MCMCtree*` directories to reproduce the summary results we report in our manuscript, please make sure that you have done the following:

* Download the compressed file with all the results obtained during the timetree inference analyses from [our archive in figshare](https://figshare.com/s/4f4ab42f6a2e0886a161).
* Once downloaded, please decompress the file and find the `00_prior` and `01_posterior` directories inside each `01_PAML/*MCMCtree*/sum_analyses` directory depending on which analysis you want to reproduce.
* You will see that the `00_prior` and `01_posterior` directories are already inside the `01_PAML/*MCMCtree*` directories that we provide you with in this GitHub repository. Nevertheless, they are all missing the `mcmc.txt` files that you shall find inside `00_prior/CLK/*/` and `01_posterior/[GBM|ILN]/*/` directories, which are the output files generated by `MCMCtree` with all the samples collected for each parameter of interest. Our in-house script used to summarise the results use such files. We recommend that you do the following:
  * Move the content of the `sum_analyses` directories elsewhere so that you can have a copy of our original results, but please leave the `sum_analyses` directory empty.
  * Find the `00_prior` and `01_posterior` directories that you will find in the decompressed file for each `01_PAML/*MCMCtree*` analysis and save them in the corresponding `sum_analysis` directory of this cloned repository. Make sure that you save them in the same `01_PAML/*MCMCtre*` directory following this structure:
    * Directories `00_prior` and `01_posterior` that you shall find when you decompress `MCMCtree_conc_cbA.zip` need to go in `02_MCMCtree`.
    * Directories `00_prior` and `01_posterior` that you shall find when you decompress `MCMCtree_conc_notcb.zip` need to go in `03_MCMCtree_nonbraced`.
    * Directories `00_prior` and `01_posterior` that you shall find when you decompress `MCMCtree_conc_cbB.zip` need to go in `04_MCMCtree_fossbraced`.
    * Directories `00_prior` and `01_posterior` that you shall find when you decompress `MCMCtree_part_cbA.zip` need to go in `05_MCMCtree_part`.
    * Directories `00_prior` and `01_posterior` that you shall find when you decompress `MCMCtree_part_notcb.zip` need to go in `06_MCMCtree_part_nonbraced`.
    * Directories `00_prior` and `01_posterior` that you shall find when you decompress `MCMCtree_part_cbB.zip` need to go in `07_MCMCtree_part_fossbraced`.
  * If you follow the instructions and run the code snippets and in-house scripts as detailed in each `README.md` file under directories `01_PAML/*MCMCtree`, you will be able to reproduce all our results! You can always compare them to our original output files that should have already cloned and saved elsewhere :)

Once you have done this, you are ready to get started with the analyses! Below, you can find a short description of what you can find in each directory:

* [`00_CODEML`](01_PAML/00_CODEML/README.md): guidelines to infer the branch lengths, the gradient, and the Hessian for [the concatenated alignment](00_data_formatting/01_inp_data/LUCAdup_246sp_aln.phy) with `CODEML`.
* [`01_CODEML`](01_PAML/01_CODEML/README.md): guidelines to infer the branch lengths, the gradient, and the Hessian for each gene alignment with `CODEML`. Remember that each gene alignment is part of [the partitioned alignment](00_data_formatting/01_inp_data/LUCAdup_246sp_5parts_aln.phy) as individual alignment blocks. with `CODEML`.
* [`02_MCMCtree`](01_PAML/02_MCMCtree/README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the concatenated dataset and under the Geometric Brownian motion ([Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal ([Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) relaxed-clock models. Note that all nodes that correspond to the same speciation event were cross-braced during this timetree inference analysis.
* [`03_MCMCtree_nonbraced`](01_PAML/03_MCMCtree_nonbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) under the Geometric Brownian motion ([Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal ([Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) relaxed-clock models with the concatenated dataset. Note that cross-bracing was not enabled during this timetree inference analysis.
* [`04_MCMCtree_fossbraced`](01_PAML/04_MCMCtree_fossbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the concatenated dataset and under the Geometric Brownian motion ([Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal ([Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) relaxed-clock models. Note that only those nodes for which a fossil calibration was available were cross-braced during this timetree inference analysis.
* [`05_MCMCtree_part`](01_PAML/05_MCMCtree_part//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the partitioned dataset and under the Geometric Brownian motion ([Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal ([Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) relaxed-clock models. Note that all nodes that correspond to the same speciation event were cross-braced during this timetree inference analysis.
* [`06_MCMCtree_part_nonbraced`](01_PAML/06_MCMCtree_part_nonbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) under the Geometric Brownian motion ([Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal ([Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) relaxed-clock models with the partitioned dataset. Note that cross-bracing was not enabled during this timetree inference analysis.
* [`07_MCMCtree_part_fossbraced`](01_PAML/07_MCMCtree_part_fossbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the partitioned dataset and under the Geometric Brownian motion ([Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and independent-rate log-normal ([Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) relaxed-clock models. Note that only those nodes for which a fossil calibration was available were cross-braced during this timetree inference analysis.

## Software

Before you start this tutorial, please make sure you have the following software installed on your PCs/HPCs:

* **`PAML`**: we used `PAML` v4.10.7, available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml). Please download [the latest release on the PAML GitHub repository which, as of January 2024, corresponds to the PAML version we used: v4.10.7](https://github.com/abacus-gene/paml/releases/tag/4.10.7). You can use either the pre-compiled binaries for your OS or compile the source code. Note that we used the Linux version to run these analyses on an HPC. Numerical differences may occur depending on the OS where you run this software.

* **`R`** and **`RStudio`**: please download [R](https://cran.r-project.org/) and [RStudio](https://posit.co/download/rstudio-desktop/) as we used them throughout the practical. The packages we should work with R versions that are either newer than or equal to v4.1.2 (we used v4.1.2). If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another R version installed on your PC, please check whether you need to install `RTools 4.2` or `RTools 4.3`. For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2
    version$version.string
    # Now, install the packages we will be using
    # Note that it may take a while if you have not 
    # installed all these software before
    install.packages( c('rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan', 'colorBlindness'), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that you have installed the `RTools`
    ## aforementioned.
    ## The versions we used for each of the packages aforementioned are the following:
    ##   rstudioapi: v0.14
    ##   ape: v5.7.1
    ##   phytools: v1.5.1
    ##   sn: v2.1.1
    ##   stringr: v1.5.0
    ##   rstan: v2.21.7
    ##   colorBlindness: v0.1.9
    ```

* **`FigTree`**: you can use this graphical interface to display the timetrees that we provide in this repository and that you can reproduce if you follow our guidelines. You can then decide what you want to be displayed by selecting the buttons and options that you require for that to happen. You can [download the latest pre-compiled binaries, `FigTree v1.4.4`, from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest pre-compiled binaries, `Tracer v1.7.2`, from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

## Data analysis

If you have gone through the previous sections and have a clear understanding of the dataset we used, the workflow we followed and you shall follow to reproduce our analyses, and have installed the required software to do so... Then you are ready to go!

You can start by taking a look at how we formatted the raw dataset with which we started the timetree inference analyses [by following this link](00_data_formatting/README.md).

Happy timetree inference! :)
