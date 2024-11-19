paper repository for Vazquez-Garcia & Obermayer et al.

contents:

- `swibrid_runs`: contains config files for various SWIBRID runs on human or mouse data, or the simulations

    - `benchmarks`:  config files for the benchmarks

        - `dense`: using dense MSA, can be run as is using `swibrid test`
        - `sparse`: using sparse MSA. for this, the [`sparsecluster`](github.com/bobermayer/sparsecluster) package needs to be installed

    - `mouse`: config files for mouse data

        - download raw fastq files from SRA and run `demultiplex_dataset.sh` for all batches (20220425, 20221118, 20221212, 20230120, 20230213, 20240830); this will put fastq and info files for individual samples into `input`
        - download mm10 genome from UCSC or elsewhere
        - download gencode M12 reference and use `swibrid prepare_annotation`
        - use `config.yaml` for running all mouse data
        - use `config_noSg.yaml` for running everything only on Sm + Sa (potentially use only reads with Sa primer from the info files in `input`)

    - `human`: config files for human data

        **raw sequencing data for human donors cannot be shared due to patient privacy legislation**

- `data`: contains input data for paper figures (prepared using `prepare_data.R`

- `paper_figures.Rmd`: R code to produce paper figures
