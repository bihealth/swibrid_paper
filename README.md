paper repository for Vazquez-Garcia & Obermayer et al.

contents:

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

- `data`: contains input data for paper figures (prepared using `prepare_data.R`

- `paper_figures.Rmd`: R code to produce paper figures

- `supplementary_note.ipynb`: python code to make plots for supplementary note
