paper repository for Vazquez-Garcia & Obermayer et al.

contents:

- `data`: contains input data for paper figures (prepared using `prepare_data.R`

- `paper_figures.Rmd`: R code to produce paper figures (uses data in `data`, nothing else needed)

<details>
<summary>sessionInfo</summary>
**R version 4.3.2 (2023-10-31)**

**Platform:** x86_64-pc-linux-gnu (64-bit) 

**locale:**
_LC_CTYPE=en_US.UTF-8_, _LC_NUMERIC=C_, _LC_TIME=en_US.UTF-8_, _LC_COLLATE=en_US.UTF-8_, _LC_MONETARY=en_US.UTF-8_, _LC_MESSAGES=en_US.UTF-8_, _LC_PAPER=en_US.UTF-8_, _LC_NAME=C_, _LC_ADDRESS=C_, _LC_TELEPHONE=C_, _LC_MEASUREMENT=en_US.UTF-8_ and _LC_IDENTIFICATION=C_

**attached base packages:** 
_grid_, _stats_, _graphics_, _grDevices_, _utils_, _datasets_, _methods_ and _base_

**other attached packages:** 
_dendextend(v.1.17.1)_, _lme4(v.1.1-35.1)_, _gtools(v.3.9.5)_, _RColorBrewer(v.1.1-3)_, _variancePartition(v.1.32.5)_, _BiocParallel(v.1.36.0)_, _limma(v.3.58.1)_, _readxl(v.1.4.3)_, _pROC(v.1.18.5)_, _glmnet(v.4.1-8)_, _Matrix(v.1.6-5)_, _car(v.3.1-2)_, _carData(v.3.0-5)_, _ggrepel(v.0.9.5)_, _circlize(v.0.4.16)_, _ComplexHeatmap(v.2.18.0)_, _cowplot(v.1.1.3)_, _scales(v.1.3.0)_, _caret(v.6.0-94)_, _lattice(v.0.21-9)_, _lubridate(v.1.9.3)_, _forcats(v.1.0.0)_, _stringr(v.1.5.1)_, _purrr(v.1.0.2)_, _readr(v.2.1.5)_, _tidyr(v.1.3.1)_, _tibble(v.3.2.1)_, _tidyverse(v.2.0.0)_, _dplyr(v.1.1.4)_, _ggpubr(v.0.6.0)_ and _ggplot2(v.3.5.1)_

**loaded via a namespace (and not attached):** 
_bitops(v.1.0-7)_, _Rdpack(v.2.6)_, _gridExtra(v.2.3)_, _rlang(v.1.1.3)_, _magrittr(v.2.0.3)_, _clue(v.0.3-65)_, _GetoptLong(v.1.0.5)_, _matrixStats(v.1.2.0)_, _compiler(v.4.3.2)_, _png(v.0.1-8)_, _vctrs(v.0.6.5)_, _reshape2(v.1.4.4)_, _pkgconfig(v.2.0.3)_, _shape(v.1.4.6.1)_, _crayon(v.1.5.2)_, _backports(v.1.4.1)_, _pander(v.0.6.5)_, _caTools(v.1.18.2)_, _utf8(v.1.2.4)_, _prodlim(v.2023.08.28)_, _tzdb(v.0.4.0)_, _nloptr(v.2.0.3)_, _xfun(v.0.42)_, _EnvStats(v.2.8.1)_, _recipes(v.1.0.10)_, _remaCor(v.0.0.18)_, _broom(v.1.0.5)_, _parallel(v.4.3.2)_, _cluster(v.2.1.4)_, _R6(v.2.5.1)_, _stringi(v.1.8.3)_, _boot(v.1.3-28.1)_, _parallelly(v.1.37.0)_, _rpart(v.4.1.21)_, _numDeriv(v.2016.8-1.1)_, _cellranger(v.1.1.0)_, _Rcpp(v.1.0.12)_, _iterators(v.1.0.14)_, _knitr(v.1.45)_, _future.apply(v.1.11.1)_, _IRanges(v.2.36.0)_, _splines(v.4.3.2)_, _nnet(v.7.3-19)_, _timechange(v.0.3.0)_, _tidyselect(v.1.2.0)_, _viridis(v.0.6.5)_, _rstudioapi(v.0.15.0)_, _abind(v.1.4-5)_, _timeDate(v.4032.109)_, _gplots(v.3.1.3.1)_, _doParallel(v.1.0.17)_, _codetools(v.0.2-19)_, _listenv(v.0.9.1)_, _lmerTest(v.3.1-3)_, _plyr(v.1.8.9)_, _Biobase(v.2.62.0)_, _withr(v.3.0.0)_, _future(v.1.33.1)_, _survival(v.3.5-7)_, _pillar(v.1.9.0)_, _KernSmooth(v.2.23-22)_, _foreach(v.1.5.2)_, _stats4(v.4.3.2)_, _generics(v.0.1.3)_, _S4Vectors(v.0.40.2)_, _hms(v.1.1.3)_, _aod(v.1.3.3)_, _munsell(v.0.5.0)_, _minqa(v.1.2.6)_, _globals(v.0.16.2)_, _RhpcBLASctl(v.0.23-42)_, _class(v.7.3-22)_, _glue(v.1.7.0)_, _tools(v.4.3.2)_, _fANCOVA(v.0.6-1)_, _data.table(v.1.15.0)_, _ModelMetrics(v.1.2.2.2)_, _gower(v.1.0.1)_, _ggsignif(v.0.6.4)_, _mvtnorm(v.1.2-4)_, _rbibutils(v.2.2.16)_, _ipred(v.0.9-14)_, _colorspace(v.2.1-0)_, _nlme(v.3.1-163)_, _cli(v.3.6.2)_, _fansi(v.1.0.6)_, _viridisLite(v.0.4.2)_, _lava(v.1.8.0)_, _corpcor(v.1.6.10)_, _gtable(v.0.3.4)_, _rstatix(v.0.7.2)_, _digest(v.0.6.34)_, _BiocGenerics(v.0.48.1)_, _pbkrtest(v.0.5.2)_, _rjson(v.0.2.21)_, _lifecycle(v.1.0.4)_, _hardhat(v.1.3.1)_, _GlobalOptions(v.0.1.2)_, _statmod(v.1.5.0)_ and _MASS(v.7.3-60)_
</details>

- `swibrid_runs`: contains config files for various SWIBRID runs on human or mouse data, or the simulations

    - `benchmarks`:  config files for the benchmarks

        - `dense`: using dense MSA, can be run as is using `swibrid test`
        - `sparse`: using sparse MSA. for this, the [`sparsecluster`](https://github.com/bobermayer/sparsecluster) package needs to be installed

    - `mouse`: config files for mouse data

        - download raw fastq files from SRA (accession PRJNA1190672) into `raw_data` and run `demultiplex_dataset.sh`; this will put fastq and `info.csv` files for individual samples into `input` and make it possible to run all samples in one go
        - download mm10 genome from UCSC or elsewhere
        - download gencode M12 reference and use `swibrid prepare_annotation`
        - use `config.yaml` for running all mouse data
        - use `config_noSg.yaml` for running everything only on Sm + Sa (potentially restrict info files in `input` to reads with Sa primer)

    - `human`: config files for human data

        **raw sequencing data for human donors cannot be shared due to patient privacy legislation**

        - `demultiplex_dataset.sh` is used to demultiplex input for each run, demultiplexed fastq and `info.csv` files would be expected in `input`
        - get hg38 genome and gencode v33 reference, create LAST index
        - `config.yaml` for "regular" runs
        - `config_reads_averaging.yaml` to use averaging of features over reads not clusters
        - `combine_replicates.sh` to pool reads from technical replicates
        - `plot_bars.sh` and `plot_bars.py` to plot isotype fractions as in Fig. 1
        - `plot_circles.sh` and `plot_circles.py` to create bubble plots of Fig. 1
        - `plot_clustering.sh` to create read plots for Fig. 1 and S2
        - `plot_breakpoints.sh` and `plot_breakpoint_stats.py` to create breakpoint matrix plot of Fig. 2A

     - `external`: config files for public datasets (Vincendeau et al. and Panchakshari et al.)

        - for Vincendeau et al., download data from SRA (PRJNA831666) into the Vincendeau subfolder and run `make_info.py` on every sample to create dummy files with primer locations
        - for Panchakshari et al., use `get_data.sh` in the `HTGTS` folder to download data, collapse read mates with `bbmerge` and create info files

- `supplementary_note.ipynb`: python code to make plots for supplementary note
