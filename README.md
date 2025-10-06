# VSG Processing Pipeline

This repository houses the tooling that takes the raw count tables produced by the [myRNA-seq pipeline](https://github.com/mtinti/myRna-seq) and transforms them into curated, quality-controlled inputs for the [VSGs web server](https://vsgs-web-server.pages.dev/). The workflow focuses on extracting differential expression signatures for Variant Surface Glycoproteins (VSGs), summarising silent VSG activation, monitoring the dominant VSGs across experiments, and compiling the quality metrics that power the web leaderboard and plots.

## Why this project exists

Researchers who run myRNA-seq on *Trypanosoma brucei* datasets receive comprehensive alignment and counting outputs, but preparing those results for comparative VSG analyses requires several extra steps. This project bridges that gap by:

- Normalising raw VSG count data and calculating fold-changes between control and treatment conditions.
- Flagging the top-expressed VSGs per sample to support experiment curation.
- Aggregating silent VSG activity to highlight derepression trends across experiments.
- Computing experiment-level quality control metrics that surface potential anomalies before visualisation.
- Emitting the CSV files consumed directly by the VSGs leaderboard and plotting interface.

## Repository structure

```
vsg_processing/
├── vsg_analysis_pipeline.py      # End-to-end orchestration of the processing workflow
├── experiment_table.txt          # Metadata table describing samples and experimental design
├── vsg_dic.txt                   # Mapping of VSG identifiers to gene families
├── de_analysis/                  # Differential expression summaries per experiment (generated)
├── plot_data.csv                 # Example aggregated output for plotting
├── vsgs_web_server/              # Static assets mirrored by the public VSGs site
├── TriTrypDB-68_TbruceiTREU927.gff# Genome annotation used to add gene descriptions
└── prepare_data_and_plots.ipynb   # Notebook for exploratory analysis and figure generation
```

> **Note**: The `de_analysis` directory is populated when the pipeline runs.

The `prepare_data_and_plots.ipynb` notebook can be run end-to-end to execute the same processing script that powers the pipeline, regenerate the curated datasets, and build the Plotly scatter plot visualisations that are embedded within the VSGs web server.

## End-to-end workflow

1. **Collect metadata** – `experiment_table.txt` captures SRA accessions, experimental groups, and run-level details. This file underpins the automation of control/treatment pairing and sample exclusion.
2. **Load VSG annotations** – `vsg_dic.txt` and the TriTrypDB GFF supply gene family labels and descriptions that enrich downstream summaries.
3. **Validate myRNA-seq outputs** – For every SRA entry, the pipeline checks the presence of count and QC files inside `myRna-seq/results/result_vsgs/<SRA>/`.
4. **Normalise and compare** – Counts are converted to relative abundances, controls and treatments are averaged across replicates, and differential expression statistics (fold change and log₂FC) are computed for every VSG.
5. **Summarise by experiment** – The script produces experiment configuration metadata (`exp_config.json`), silent VSG activity tables (`silentCsvData.csv`), main VSG dynamics (`mainCsvData.csv`), and QC summaries (`QC.csv`). Individual differential expression tables are written to `de_analysis/` for auditability.
6. **Publish-ready outputs** – The generated CSVs feed directly into the VSGs web server leaderboards and plots, ensuring the website reflects the most recent analyses while keeping provenance with the raw counts.

## Running the pipeline

1. **Install dependencies** (Python ≥3.8):
   ```bash
   pip install pandas numpy tqdm
   ```
2. **Fetch myRNA-seq results** by running the upstream pipeline or copying the relevant `result_vsgs` directory into this repository.
3. **Execute the processor**:
   ```bash
   python vsg_analysis_pipeline.py \
     --experiment-table experiment_table.txt \
     --vsg-dict vsg_dic.txt \
     --base-path myRna-seq/results/result_vsgs \
     --output-dir output
   ```

The script logs progress and writes all derived files into the specified `output` directory. Regenerating the summary tables allows updating the public-facing dashboards in a single command once fresh counts are available.

## Quality assurance focus

- **Automated exclusions**: Known problematic accessions are filtered to prevent contaminated comparisons.
- **QC ratios**: The pipeline records the dominance of the main VSG and mitochondrial read fractions, helping identify mixed populations or sequencing artefacts.
- **Traceable provenance**: All experiment-level differential expression tables remain accessible for manual review (`de_analysis/*.csv`).

## Acknowledgements
- TriTrypDB for providing the genomic annotations powering gene descriptions.
