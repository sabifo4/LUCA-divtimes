[![DOI](https://zenodo.org/badge/745536810.svg)](https://zenodo.org/doi/10.5281/zenodo.11260522)

# The nature of the Last Universal Common Ancestor and its impact on the early Earth system -- timetree inference analyses

In this repository, you will find a very detailed tutorial with all the steps you need to follow to reproduce the results for the timetree inference analyses we carried out as part of our study: **The nature of the Last Universal Common Ancestor and its impact on the early Earth system**.

To get started, you can clone this repository on your PC and follow all the guidelines given in the various `README.md` files that you shall find inside each directory so that you can go through all the steps we carried out for timetree inference. A summary of this workflow is given below:

* Parsing and formatting the input data required to run `PAML` programs `BASEML` and `MCMCtree` ([Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17483113/)): sequence files (concatenated and partitioned), calibrated tree files (three calibration strategies, more details below), and control files.
* Inferring the mean evolutionary rate to specify a sensible rate prior.
* Running `PAML` programs for timetree inference:
  * Using various in-house pipelines to set up the working environment, the file structure, and the control files required to run `PAML` programs.
  * Running `CODEML` to calculate the branch lengths, the gradient, and the Hessian; required by `MCMCtree` to enable the approximate likelihood calculation [dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613).
  * Running `MCMCtree` with the approximate likelihood calculation enabled for timetree inference. Our main analysis consisted of inferring the timetree for the **partitioned dataset** under both the Geometric Brownian motion (**GBM**, [Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and the independent-rates log-normal (**ILN**, [Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) **relaxed-clock models**, where nodes that correspond to the same divergences are cross-braced (i.e., hereby referred to as calibration strategy “cross-bracing A”). In addition, we ran 10 additional inference analyses to benchmark the effect that partitioning, cross-bracing, and relaxed-clock models can have on species divergence times estimation:
    * **GBM + concatenated alignment + cross-bracing A**.
    * **GBM + concatenated alignment + cross-bracing B** (i.e., only nodes that correspond to the same divergences for which there are fossil constraints are cross-braced).
    * **GBM + concatenated alignment + without cross-bracing**.
    * **GBM + partitioned alignment + cross-bracing B**.
    * **GBM + partitioned alignment + without cross-bracing**.
    * **ILN + concatenated alignment + cross-bracing A**.
    * **ILN + concatenated alignment + cross-bracing B**.
    * **ILN + concatenated alignment + without cross-bracing**.
    * **ILN + partitioned alignment + cross-bracing B**.
    * **ILN + partitioned alignment + without cross-bracing**.
* Running MCMC diagnostics for all the chains under each analysis.

To make it easier for you to navigate this GitHub repository, below you can find a summary of the content you shall find inside each directory.

## [Data formatting](00_data_formatting/README.md)

* **Molecular alignment**: we analysed a molecular alignment with 5 pre-LUCA paralogues (i.e., F- and V- type subunits from ATPases; Elongation Factor Tu and G; Signal Recognition Protein and Signal Recognition Particle Receptor; Tyrosyl-tRNA and Tryptophanyl-tRNA synthetases; and Leucyl- and Valyl- tRNA synthetases) for 246 taxa. You can read the `README.md` file under [`00_data_formatting`](00_data_formatting/README.md) to reproduce all the steps that we carried out to parse the raw data.
  * [`Concatenated molecular alignment (raw)`](00_data_formatting/00_raw_data/alignment/concat5.fas): the raw molecular alignment consisted of 1,080 characters and 247 taxa. Gene families were identified using `BLASTp` ([Altschul 1990](https://pubmed.ncbi.nlm.nih.gov/2231712/)) and downloaded from the NCBI ([Sayers et al., 2010](https://pubmed.ncbi.nlm.nih.gov/21097890/)). Amino acid sequences were aligned with MUSCLE ([Edgar 2004](https://academic.oup.com/nar/article/32/5/1792/2380623)) and trimmed with `TrimAl` ([Capella-Gutiérrez 2009](https://academic.oup.com/bioinformatics/article/25/15/1972/213148)) (`-strict`). Individual gene trees were inferred under the LG+C20+F+G substitution model implemented in `IQ-TREE 2` ([Minh et al. 2020](https://academic.oup.com/mbe/article/37/5/1530/5721363?login=false)). These trees were manually inspected and curated to remove non-homologous sequences or recent horizontal gene transfers. According to the filters applied, the resulting gene alignments included 246 species, and multiple copies of each taxa were represented at least twice (for some eukaryotes, they may be represented by plastid, mitochondrial, and nuclear sequences).
  * [`Individual gene alignments for each paralog (raw)`](00_data_formatting/00_raw_data/alignment/partitioned/): you will find the five gene alignments (i.e., one for each paralog) under the [`partitioned` directory](00_data_formatting/00_raw_data/alignment/partitioned/).
  * [`Input concatenated alignment`](00_data_formatting/01_inp_data/LUCAdup_246sp_aln.phy): after parsing the raw concatenated alignment aforementioned, this is the input concatenated alignment that we used for all timetree inference analyses.
  * [`Input partiitoned alignment`](00_data_formatting/01_inp_data/LUCAdup_246sp_5parts_aln.phy): after parsing the raw individual gene alignments aforementioned, this is the input partitioned alignment (i.e., five alignment blocks, one for each gene alignment) that we used for all timetree inference analyses.
* **Phylogeny**: given that a fixed tree topology is required for timetree inference with `MCMCtree`, we decided to infer the best-scoring ML tree with `IQ-TREE 2` under the LG+C20+F+G4 model ([Hoang 2018](https://pubmed.ncbi.nlm.nih.gov/29077904/)) following our previous phylogenetic analyses. We then modified the resulting inferred tree topology to implement the available fossil calibrations (see below) and following consensus views of species-level relationships ([Coleman et al., 2021](https://pubmed.ncbi.nlm.nih.gov/33958449/), [Aouad et al., 202](https://bmcecolevol.biomedcentral.com/articles/10.1186/s12862-021-01952-0), [Burki et al., 2019](https://www.cell.com/trends/ecology-evolution/fulltext/S0169-5347(19)30257-5)).
* **Calibrations**: we used previously established dates for four fossil calibrations and included nine additional calibrations with justifications provided in the supplementary material of the manuscript. In addition, five of the newly defined calibrations (i.e., total group Eukarya for the archaeal lineage, Embryophyta, Archaeplastida, and Eumetazoa) have been revised to follow the latest geochronological updates and literature as of December 2023. Note that the root of a phylogeny requires constraining in molecular clock-dating analyses. Therefore, we constrained the root (representing the duplication of pre-LUCA paralogues) with a hard upper bound equal to the age of the moon forming impact. We generated three control files that we used as input fileS for our R in-house scripts to semi-automatically calibrate the tree topologies with the `MCMCtree` notation:
  * [`Calib_converter_allcb.txt` file, strategy "cross-bracing A"](00_data_formatting/scripts/Calib_converter_allcb.txt): text file that our [R in-house script `Include_calibrations_allcb.R`](00_data_formatting/scripts/Include_calibrations_allcb.R) uses as a "converter" file to calibrate the fixed tree topology aforementioned using `MCMCtree` notation. Note that the calibration format used in this file does not follow the conventional `MCMCtree` notation. Specifically, the format you shall find is used to cross-brace those nodes with the same speciation event, and so the posterior density for such "mirrored" nodes will be the same.
  * [`Calib_converter_allnoncb.txt` file, strategy "no cross-bracing"](00_data_formatting/scripts/Calib_converter_allnoncb.txt): text file that our [R in-house script `Include_calibrations_allnoncb.R`](00_data_formatting/scripts/Include_calibrations_allnoncb.R) uses as a "converter" file to calibrate the fixed tree topology aforementioned using `MCMCtree` notation. Note that the format used in this file is the conventional `MCMCtree` format in which the age of each node is constrained by a specific probability density. In that way, cross-bracing is not enabled and each node shall have the corresponding posterior density (i.e., no mirrored nodes).
  * [`Calib_converter_fosscb.txt` file, strategy "cross-bracing B"](00_data_formatting/scripts/Calib_converter_fosscb.txt): text file that our [R in-house script `Include_calibrations_fosscb.R`](00_data_formatting/scripts/Include_calibrations_fosscb.R) uses as a "converter" file to calibrate the fixed tree topology aforementioned using `MCMCtree` notation. There is one difference between this calibration strategy and "cross-bracing A": only those nodes for which there is a fossil calibration have been cross-braced.

## [PAML analyses](01_PAML)

Note that all the directories that you shall find inside the [`01_PAML` directory](01_PAML) have their own `README.md` file. In that way, you can navigate to each of these directories to read more about how we ran `CODEML` and `MCMCtree` when using either the concatenated or the partitioned alignments and when running the inference analyses under different settings (e.g., data type, calibration strategy, relaxed-clock model).

Below, you can find a short description of what you can find in each directory:

* [`00_CODEML`](01_PAML/00_CODEML/README.md): guidelines to infer the branch lengths, the gradient, and the Hessian for [the concatenated alignment](00_data_formatting/01_inp_data/LUCAdup_246sp_aln.phy) with `CODEML`.
* [`01_CODEML`](01_PAML/01_CODEML/README.md): guidelines to infer the branch lengths, the gradient, and the Hessian for each gene alignment with `CODEML`. Remember that each gene alignment is part of [the partitioned alignment](00_data_formatting/01_inp_data/LUCAdup_246sp_5parts_aln.phy) as individual alignment blocks.
* [`02_MCMCtree`](01_PAML/02_MCMCtree/README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and **the newly implemented cross-bracing approach** with the concatenated dataset under the GBM and ILN relaxed-clock models. Note that all nodes that correspond to the same speciation event were cross-braced during this timetree inference analysis (i.e., "cross-bracing A").
* [`03_MCMCtree_nonbraced`](01_PAML/03_MCMCtree_nonbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) under the GBM and ILN relaxed-clock models with the concatenated dataset. Note that cross-bracing was not enabled during this timetree inference analysis.
* [`04_MCMCtree_fossbraced`](01_PAML/04_MCMCtree_fossbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the concatenated dataset under the GBM and ILN relaxed-clock models. Note that only those nodes for which a fossil calibration was available were cross-braced during this timetree inference analysis (i.e., "cross-bracing B").
* [`05_MCMCtree_part`](01_PAML/05_MCMCtree_part//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the partitioned dataset under the GBM and ILN relaxed-clock models. Note that all nodes that correspond to the same speciation event were cross-braced during this timetree inference analysis (i.e., "cross-bracing A").
* [`06_MCMCtree_part_nonbraced`](01_PAML/06_MCMCtree_part_nonbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) under the GBM and ILN relaxed-clock models with the partitioned dataset. Note that cross-bracing was not enabled during this timetree inference analysis.
* [`07_MCMCtree_part_fossbraced`](01_PAML/07_MCMCtree_part_fossbraced//README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) and cross-bracing with the partitioned dataset under the GBM and ILN relaxed-clock models. Note that only those nodes for which a fossil calibration was available were cross-braced during this timetree inference analysis (i.e., "cross-bracing B").

Please note that each directory has all the input data, intermediate files, and scripts that you will need to follow the step-by-step guidelines if you clone this repository on your PC. Nevertheless, note that the `MCMCtree` output files with all the collected samples during the MCMC that we use to summarise our results are too large to be stored here. When you want to run the tutorials under the `01_PAML/*MCMCtree*` directories to reproduce the summary results we report in our manuscript, please make sure that you have (i) either generated these output files yourself so that the scripts can run or (ii) do the following to download our results and use the scripts with our own results:

* Download the compressed file with all the results we obtained during the timetree inference analyses from [our archive (~71Gb)](https://www.dropbox.com/scl/fo/evl05sptqkgodvbfsppze/AGI5VAXPqyfMK5Hp4pQwx6o?rlkey=512vrq492xvd120jcqywmkzd7&st=jpv4919a&dl=0). Please note that we may change our archive location, and so we ask that you **keep an eye on [our GitHub repository](https://github.com/sabifo4/LUCA-divtimes) to track all the changes/updates since the publication of our manuscript**.
* Once downloaded, please decompress the file and find the `00_prior` and `01_posterior` directories inside each `01_PAML/*MCMCtree*/sum_analyses` directory (depending on which analysis you want to reproduce).
* You will see that the `00_prior` and `01_posterior` directories are already inside the `01_PAML/*MCMCtree*` directories that we provide you with in this GitHub repository. Nevertheless, they are all missing the `mcmc.txt` files that you shall find inside `00_prior/CLK/*/` and `01_posterior/[GBM|ILN]/*/` directories, which are the output files generated by `MCMCtree` with all the samples collected for each parameter of interest. Our in-house scripts use such files to summarise the results. We recommend that you do the following:
  * If you have your own results but want to reproduce our output, please either recreate the file structure and save the content of our `sum_analyses` directories that you will have downloaded or save your results elsewhere and reuse the file structure in this repository. If you just want to reanalyse our results, then keep reading.
  * Find the `00_prior` and `01_posterior` directories in the decompressed file for each `01_PAML/*MCMCtree*` analysis and save them in the corresponding `sum_analyses` directory in this cloned repository. Make sure that you save them in the same `01_PAML/*MCMCtre*` directory following the file structure!
  * If you follow the instructions and run the code snippets and in-house scripts as detailed in each `README.md` file under directories `01_PAML/*MCMCtree`, you will be able to reproduce all our results! If you had re-run these analyses, then you will be able to compare your results to ours!

## [Figures and tables](02_ms_figs_tables)

The figures and tables that we included in our supplementary file can also be found in this directory. Please refer to the supplementary material document for figure and table captions (i.e., some numbers may have been updated during the final publication process when rearranging figures/tables, so the numbers in this repository may not correspond to the final ones. **Please keep an eye on [our GitHub repository](https://github.com/sabifo4/LUCA-divtimes) to track all the changes/updates since the publication of our manuscript**). All the plots an tables included in these documents are generated by the scripts you will find inside [the `figs_tables` directory](01_PAML/figs_tables/scripts/).

## [Additional analyses](03_additional_analyses)

We carried out three additional sensitivity analyses to address some questions asked during the revisions of our manuscript:

* [Single-gene alignments strategy](03_additional_analyses/00_indgenes_analyses/): we ran a timetree inference analysis for each paralog separately, that is, a total of 5 timetree inference analysis for each gene alignment: "ATP", "EF", "Leu", "SRP", and "Tyr". The link above will redirect you to the workflow we followed to carry out these analyses.
* ["Leave-one-out" strategy](03_additional_analyses/01_lou_analyses/): we inferred 5 gene alignments in which each paralog had been iteratively removed (e.g., remove "ATP" but keep genes "EF", "Leu", "SRP", and "Tyr" concatenated in a unique alignment block; follow the same procedure for each gene family). The link above will redirect you to the workflow we followed to carry out these analyses.
* ["bs_inBV" strategy](03_additional_analyses/02_blunc_analysis/): we used the newly [`bs_inBV` approach](https://github.com/evolbeginner/bs_inBV) to assess the effect that using branch lengths, gradient, and Hessian estimated under a complex substitution model for timetree inference could have on our time estimates. The link above will redirect you to the workflow we followed to carry out these analyses.

After running all these analyses, we then generated summary plots and tables to compare the estimates we obtained for these sensitivity tests with our core analyses. You can find these files in the [`figs_tables`](03_additional_analyses/figs_tables/) together with the scripts we used to generate them.

## [In-house R functions](src)

The main functions that you call when you execute the various in-house scripts that are mentioned in this tutorial can be found in these three R scripts:

* [Functions.R](src/Functions.R)
* [Functions_plots_MCMCtreeR.R](src/Functions_plots_MCMCtreeR.R)
* [Functions_plots.R](src/Functions_plots.R)

If you open the scripts aforementioned, you will find all the details regarding the arguments they need and how to run them (which is also referred to when they are called when summarising our results and/or generating output files such as plots or tables).

## What do you need to run our analyses?

### Software

Before you go through this step-by-step tutorial to reproduce our results, please make sure you have the following software installed on your PCs/HPCs:

* **`PAML`**: we used `PAML` v4.10.7, available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml). Please download [the latest release on the PAML GitHub repository which, as of January 2024, corresponds to the PAML version we used: v4.10.7](https://github.com/abacus-gene/paml/releases/tag/4.10.7). You can use either the pre-compiled binaries for your OS or compile the source code. Note that we used the Linux version to run these analyses on an HPC. Numerical differences may occur depending on the OS where you run this software.

  > **NOTE for Linux users**: you may need to install the `intel` compilers before you run `make -f Makefile`. Please visit this link to download the [Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#dpcpp-cpp). Thank you, [@ERRMoody](https://github.com/ERRMoody), for pointing this out!

* **`R`** and **`RStudio`**: please download [R](https://cran.r-project.org/) (i.e., select either `Download R for Linux`, `Download R for macOS`, or `Download R for Windows` depending on your OS; then follow the installation instructions) and [RStudio](https://posit.co/download/rstudio-desktop/) as you will need them to run various of our in-house R scripts. The packages you will need to install work with R versions that are either newer than or equal to v4.1.2 (we used v4.1.2). If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another R version installed on your PC, please check whether you need to install `RTools 4.2` or `RTools 4.3`. For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2
    version$version.string
    # Now, install the packages you will need to run
    # our in-house R scripts
    # Note that it may take a while to install them all
    install.packages( c('rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan', 'colorBlindness'), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that you have installed the `RTools`
    ## aforementioned. Otherwise, you will get errors during the installation.
    ## The versions we used for each of the packages aforementioned are the following:
    ##   rstudioapi: v0.14
    ##   ape: v5.7.1
    ##   phytools: v1.5.1
    ##   sn: v2.1.1
    ##   stringr: v1.5.0
    ##   rstan: v2.21.7
    ##   colorBlindness: v0.1.9
    ```
  
  > **NOTE for Linux users**: you may experience some problems when you try to execute the `install.packages()` command that you see in the code snippet above. If that is the case, please install each package separately (e.g., `install.packages( 'rstudioapi' )` to install `rstudioapi`, etc.). If you experience problems with `stringr`, please follow the [suggestions given in this StackOverflow page](https://stackoverflow.com/questions/38987157/libicu-and-stringi-on-fedora-24-causing-r-headaches/39411793#39411793) -- despite the question being addressed for Fedora, the solution suggested also works for Ubuntu. Thank you, [@ERRMoody](https://github.com/ERRMoody), for pointing this out!

* **`FigTree`**: you can use this graphical interface to display the timetrees that we provide in this repository and that you can reproduce if you follow our guidelines. You can then decide what you want to be displayed by selecting the buttons and options that you require for that to happen. You can [download the latest pre-compiled binaries, `FigTree v1.4.4`, from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest pre-compiled binaries, `Tracer v1.7.2`, from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

### Additional note for Mac users when using the `sed` command

If you are running a UNIX-based system (e.g., Mac users), you will experience some problems with the Linux-based `sed` command that you shall see (i) in the various code snippets that you need to run as part of the tutorial described in this repository and (ii) used in most of our in-house bash scripts. By default, this command is different from Linux-based systems, and hence you will have problems to execute it as intended unless you follow one of these two approaches:

* (EASY): Instead of running the `sed` commands using the format `sed -i 's/PATTERN/REPLACEMENT/'` (i.e., what you shall see in the code snippets), you should include `''` between `-i` and `'s/PATTERN/REPLACEMENT/'`: `sed -i '' 's/PATTERN/REPLACEMENT/'`. Please remember to modify the commands in this tutorial accordingly before you copy and paste them on your terminal!
* (MORE DIFFICULT): You should install `GNU sed` and establish it as the "standard" `sed` command instead of the one you will have by default. At the time of writing, [this post](https://medium.com/@bramblexu/install-gnu-sed-on-mac-os-and-set-it-as-default-7c17ef1b8f64) is available and has a detailed explanation that you could follow to install `GNU sed`. Nevertheless, there are many other tutorials out there that you can follow to achieve the same goal. If you decide to follow this approach, you will not need to change the commands in this tutorial!

Once the `sed` incompatibility is sorted out, there should not be any problems with regards to following this tutorial and/or running our in-house bash scripts!

<br>

----

<br>

If you have gone through the previous sections and have a clear understanding of the data we used, the workflow we followed (which you shall follow to reproduce our analyses), and have installed the required software aforementioned... Then you are ready to go!

You can start by taking a look at how we formatted the raw dataset with which we started the timetree inference analyses [by following this link](00_data_formatting/README.md).

Happy timetree inference! :)

----

## Contact

This repository and its content was created by **Dr Sandra Álvarez-Carretero** ([`@sabifo4`](https://github.com/sabifo4/)), who actively maintains this repository too. If you have any queries with regards to the analyses detailed in this repository, please do not hesitate to <a href="mailto:sandra.ac93@gmail.com">reach me via email</a>!
