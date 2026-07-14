# MAP — Metabarcoding Analysis Pipeline

MAP is a fully containerized metabarcoding pipeline that supports both **Illumina paired‑end** and **Oxford Nanopore (ONT) / long‑read single‑end** data. All tools (incl. R and Python packages, etc. ) and reference databases are baked into the container or automatically created upon running MAP. MAP is customizable using command line flags and a customizable parameters file designed for customizing parameters for different barcode markers, UMIs, and primers, if desired. 

---

## Contents
- [What MAP does](#what-map-does)
- [How MAP works](#how-map-works)
- [Requirements](#requirements)
- [Quick start (run the included demo after Docker setup)](#quick-start-run-the-included-demo)
- [Run on your own data](#run-on-your-own-data)
- [Interactive container](#interactive-container)
- [Parameters](#parameters)
- [Outputs](#outputs)
- [Build the image yourself](#build-the-image-yourself)
- [Reference library](#reference-library)
- [Demo data](#demo-data)
- [Troubleshooting](#troubleshooting)
- [Citation](#citation)

---

## What MAP does
- MAP designates OTUs, taxonomy and BIN classifications for use with BOLD (Barcode of Life Database) from one to several non-demultiplexed fastq.gz files. MAP handles demultiplexing, primer/UMI trimming, denoising/clustering, reference‑frame sequence correction (COI), chimera screening, contaminant removal, taxonomic assignment (SINTAX), and BOLD BIN matching.

- MAP was specifically developed for use with COI barcodes and for BIN assignment for direct output to BOLD (Barcode of Life Database).
---

## How MAP works

1. **Demultiplex & trim** — assign reads to samples; remove primers/UMIs (`cutadapt`, `seqkit`). Illumina reads are merged (`PEAR`).
2. **Per‑sample clustering** — denoise/cluster each sample into OTUs (`vsearch`); filter by *Min reads per OTU*.
3. **Run‑wide clustering** — cluster per‑sample OTUs across the run into run‑wide OTUs.
4. **Sequence correction** (COI) — reference‑frame indel/homopolymer correction to clean coding sequences.
5. **Chimera screening** — `vsearch --uchime_denovo`.
6. **Taxonomy** — SINTAX assignment against the bundled reference library.
7. **Contaminant removal** — negative‑control‑aware filtering + replicate‑support checks.
8. **BIN matching** (COI) — match OTUs to BOLD BINs.
9. **Report** — Quarto renders the interactive HTML report; final tables written to Excel/TSV.
---


## Requirements
- Supported architectures: **linux/amd64** and **linux/arm64** (Apple Silicon).
- Disk space and memory. Large intermediate files may be created during runtime.
- [Docker](https://docs.docker.com/get-docker/) (Desktop or Engine). 

### Instructions for easy Docker setup
- There are many ways to get Docker up and running, including:

**Easiest Docker setup (Linux):**

``` bash
curl -fsSL https://get.docker.com | sh      # installs Docker Engine + compose plugin
sudo usermod -aG docker "$USER"             # so you don't need sudo
newgrp docker                               # apply group now (or just log out/in)
docker run hello-world                      # verify, no sudo
```

**Easiest Docker setup (Mac/Windows with WSL):**

- Download Docker Desktop from their website and follow their guidelines.
https://www.docker.com/products/docker-desktop/
---

## Quick start (run the included demo after Docker setup)
- **Runs the demo out of the box.** With no arguments, MAP runs the bundled `PHAUS_1K` test dataset, so you can confirm the whole pipeline works before touching your own data.

**Note:** The instructions are meant to be run from the '<PATH>/workdir' directory, but can be run from any directory with the fastq.gz, parameters.xlsx, and compose.yaml files. These 3 files are required for MAP to run. 

**Pull the image:**
```bash
docker pull ghcr.io/cbg-innov/map:latest
```
### Demo data

Bundled inside the image (under `/MAP/Metabarcoding/`) and easily viewable from 'workdir' directory in Github repo:
- `PHAUS_1K_RawReads.fastq.gz` — small COI test read set.
- `parameters.xlsx` — matching parameters spreadsheet (use as your template).
- `compose.yaml` - compose file to aid command line 

Running MAP with no `--fastq`/`--params` runs this demo end‑to‑end.

**Run the bundled demo and write results to the directory** — this single command runs the entire pipeline on the included `PHAUS_1K` test data.
**Note:** The '/data/' directory name is used internally within the Docker container, but 'workdir' is how it will apppear on your device

```bash
cd <PATH>/workdir #or rename
MAP_DATA="$(pwd)" \
  docker compose -f compose.yaml run --rm map \
  bash /MAP/SCRIPTS/MAP.sh \
  --wd /data
```
#### Here, the fastq.gz and parameters files are internal and do not need to be named directly. 

When it finishes, your results are in `./MAP_results/` (Excel + TSV tables, an interactive HTML report, and demultiplexing summaries). See [Outputs](#outputs).

---

## Run on your own data

Put your reads (`*.fastq.gz`) and a filled‑in parameters spreadsheet in a folder, mount it, and point MAP at them. Use `--wd` so all outputs land in your mounted folder:

```bash
# host folder containing: reads.fastq.gz,  parameters.xlsx, and compose.yaml
MAP_DATA="$(pwd)" docker compose -f compose.yaml \
  run --rm map \
  bash /MAP/SCRIPTS/MAP.sh \
  --fastq /data/*.fastq.gz \
  --params /data/my_parameters.xlsx \
  --wd /data
```

Results appear in `./my_run/output/`.

- **Illumina paired‑end:** Set *Paired‑end Reads = Yes* in the parameters spreadsheet. Use the common prefix and/or suffix for both files, with '*' where names diverge (R1/R2 handling and merging are driven by the parameters file.)
- **ONT / long‑read:** set *Paired‑end Reads = No* in the parameters spreadsheet.
- Start from the bundled spreadsheet as a template (in 'workdir' directory in Github repo) (see [Demo data](#demo-data)).

---

## Interactive container

For exploring intermediate files or running repeatedly, keep a container alive. `compose.yaml` also persists the reference library in a named volume so it survives container recreation.

```bash
docker compose up -d            # start a long‑running 'map' container
docker compose exec -it map bash        # drop into a shell (the 'map' env auto‑activates)

# inside the container:
bash SCRIPTS/MAP.sh        # runs the demo, or add --fastq/--params/--wd

# copy results out to your Desktop:
docker cp map:/MAP/Metabarcoding/output ~/Desktop/MAP_results

docker compose down             # stop & remove the container (volume persists)
```

To make your own data visible inside the interactive container, add a bind mount under the `map` service in `compose.yaml`, e.g.:
```yaml
    volumes:
      - reflib:/MAP/REFS
      - ~/Desktop/my_run:/data        # <-- your data here
```

---

## Parameters
**Note:** UMI map and primers are provided by the user in an Excel (.xlsx) file, using the template provided. 
The 'UMIs and Primers', 'Bulk Sample Metadata', (except Instructions) must be filled out.

### The parameters spreadsheet
Important parameters are set in the `.xlsx` parameters file. Start from the bundled `parameters.xlsx` in the workdir directory on Github to understand formatting. Key fields:

| Tab | Field | Meaning |
|---|---|---|
| *UMIs and Primers* | *Many* | *Primer and UMI information* |
| | **Run Name** | Name used for outputs |
| | **Plate** | Indicate plate name or number (e.g, Plate06)
| | **Well** | Indicate well number (e.g., B09)
| | **Sample** | Indicate sample name corresponding to well (e.g., GMP-58226_Rep1). Any replicate should be denoted with '_Rep#', but the first does not require this suffix.
| | **Forward / Reverse UMI** | UMI sequences corresponding to sample (e.g., ACAGATTTTTA) |
| | **Forward / Reverse Primer Name** | Use the names that match the sequences as defined in the Dictionary Update tab (see below) (e.g., PHAUS_F2, PHAUS_R3) |
| | **Negative Control** | indicate whether it is a *Negative* control (e.g., yes). Otherwise, leave blank. |
| | |
| *Bulk Sample Metadata* | *Sample collection information* | *Per‑sample metadata for the report* |
| | **Sample** | Should match Sample in 'UMIs and Primers' tab, minus any replicate denotion (e.g., should be GMP-58226, not GMP-58226_Rep2)|
| | **Collection Site** | Any name (e.g., Australia) |
| | **Latitude/Longitude** | Please enter separately in respective fields in decimal format (e.g., -134.5544) |
| | **Collection Start Date/Collection End Date** | Please enter separately in respective fields in this format: DD-Mon-YY (e.g., 20-Dec-24) |
| | |
| *Dictionary Update* | *Many* | *Parameter informaton to be used with each primer combination* |
| | **Forward / Reverse Primer Name** | This is where you indicate the name you will use for the corresponding Forward / Reverse Primer Sequences that are used in the 'UMIs and Primers' tab (e.g., PHAUS_F2, PHAUS_R3) |
| | **Forward / Reverse Primer Sequence** | This is where you indicate the sequence that corresponds with the primer name you provided  (e.g., AYATRGCHTTYCCHCG) |
| | **Marker** | Locus name (e.g. `COI-5P`) |
| | **Reference Library** | SINTAX‑formatted reference DB name (e.g., BOLDistilled_COI_Apr2026). Note: BOLDistilled_COI_ can be used to allow for updated/different libraries to be used without errors. |
| | **Min / Max Amplicon Length** | Length filter  WITHOUT UMIs or primers or single reads (e.g., ONT) WITH UMIs and primers attached. |
| | **Target amplicon length** | Expected amplicon length (no primers/UMIs) |
| | **Paired‑end Reads** | `Yes` (Illumina paired-end) or `No` (ONT/long‑read) |
| | **Min reads per OTU** | OTUs below this read count are discarded (e.g., 2) |
| | **Replicates per sample** | Number of replicates per sample (e.g., 8) |
| | **Intra-OTU Clustering Threshold** | Threshold for clustering OTUs within samples (i.e., across replicates of the same sample) (e.g., 2.5). Somewhat akin to denoising ASVs. |
| | **Inter-OTU Clustering Threshold** | Threshold for final OTU clustering across samples (e.g., 2.3). |


### Command‑line flags
Override paths and advanced parameters at run time (defaults shown):

| Flag | Default | Description |
|---|---|---|
| `--fastq` | `/MAP/Metabarcoding/PHAUS_1K_RawReads.fastq.gz` | Input reads. You may use a wildcard to refer to multiple files, but assign a similar prefix (e.g., `PHAUS_Illumina_*.fastq.gz`). |
| `--params` | `/MAP/Metabarcoding/parameters.xlsx` | Parameters spreadsheet |
| `--refs` | `/MAP/REFS` | Reference library directory |
| `--wd` | `/MAP/Metabarcoding` | Working directory (outputs go to `<wd>/output`) |
| `--scripts` | `/MAP/SCRIPTS` | Pipeline scripts directory |
| `--sintax_cutoff` | `0.6` | SINTAX confidence cutoff (0–1) |
| `--componentreads` | `0` | Save per‑OTU component reads (`1` yes / `0` no) |
| `--cores_to_leave` | `2` | How many cores to leave free. MAP will use the rest. |
| `--ref_seq_corr` | `/MAP/REFS/reference_seqs_327K.fasta` | File used for sequence correction |
| `--umi_overlap_min` | `0.75` | Multiplier applied to UMI lengths during demultiplexing. We recommend 0.75 for UMIs ≥12 nucleotides, and 1.0 for UMIs <12 nucleotides. |
| `--primer_overlap_min` | `0.75` | Multiplier applied to primer lengths. We recommend 0.75. |
| `--error_umi1` | `0.125` | Max error rate (Cutadapt) for the forward UMI; could be lower for shorter UMI sequences (e.g., <10bp) |
| `--error_umi2` | `0.125` | Max error rate (Cutadapt) for the reverse UMI; could be lower for shorter UMI sequences (e.g., <10bp) |
| `--error_primer1` | `0.2` | Max error rate (Cutadapt) for the forward primer |
| `--error_primer2` | `0.2` | Max error rate (Cutadapt) for the reverse primer |
| `--maxee` | `1` | (Illumina paired‑end only) Max expected error value permitted |
| `--minqual` | `10` | Min quality score permitted with Chopper (ONT/long‑read). We recommend 10 for long‑read sequences with lower overall predicted quality for better retention; a score closer to 20 is likely more useful for higher‑quality sequences (e.g., HiFi PacBio). |
| `--min_read_and_primer_length` | `100` | Min length of sequence+primers kept by Chopper. Kept loose to account for variable markers with variable polymorphic lengths. |
| `--max_read_and_primer_length` | `1000` | Max length of sequence+primers kept by Chopper. Kept loose to account for variable markers with variable polymorphic lengths. |
| `--Ill_abskew` | `10` | (Paired‑end/Illumina only) VSEARCH `uchime_denovo` abskew parameter |
| `--Ill_mindiv` | `0.0005` | (Paired‑end/Illumina only) VSEARCH `uchime_denovo` mindiv parameter |
| `--LR_abskew` | `10` | (Long‑read/ONT only) VSEARCH `uchime_denovo` abskew parameter |
| `--LR_mindiv` | `0.0005` | (Long‑read/ONT only) VSEARCH `uchime_denovo` mindiv parameter |
| `--minsize_unoise` | `2` | Min cluster size (within samples) for paired‑end reads using VSEARCH `cluster_unoise` |
| `--min_rep_prop` | `0.65` | Min proportion of replicate samples an OTU must be present in to be retained |
| `--high_dens_prop` | `0.65` | Min proportion of replicates (where the OTU is present) whose sequence count falls within the high kernel‑density area, for the OTU to be retained. For Illumina, we recommend 0.15 |
| `--alpha_default` | `0.001` | Default contamination rate; only used if `alpha_quantile` is not defined. |
| `--alpha_quantile` | `0.9` | Quantile of the calculated alpha (contaminated‑reads floor estimate) to use. For Illumina, we recommend 0.8 |
| `--k_multiplier` | `1` | Contaminant‑floor multiplier (`<1` permissive, `>1` strict). Note: we do not recommend changing this unless widespread contamination occurred. Values over 1 are likely to remove a lot of good data. |
| `--BIN_percent_ID` | `0.85` | We recommend 0.85 for short‑read data. Threshold to keep matches to BINs in the final database, regardless of 'formal' BIN assignment. Matches below this threshold are labeled "NO MATCH" in the final data. |
| `--BIN_maxaccepts` | `3` | Parameter fed to VSEARCH's `usearch_global` command. |
| `--BIN_maxhits` | `3` | Parameter fed to VSEARCH's `usearch_global` command. |

> **Tip (ONT):** raising **Min reads per OTU** is an effective lever for trimming low‑read long‑read error variants and tightening per‑sample richness.

---

## Outputs

All results are written to `<working_dir>/output/` (e.g. `./MAP_results/` or `./my_run/output/`):

- **`Metabarcoding Results - <run>.xlsx`** — main results workbook (*By Sample*, *By Replicate*, *Sample Metadata* sheets).
- **`2-TSV Versions of Results/`** — the same tables as plain TSV (`…_BySample.tsv`, `…_ByReplicate.tsv`).
- **`1-Results and Report/MAP Report - <run>.html`** — interactive report (richness, maps, treemaps, BIN matches).
- **`Demultiplexing_Results_<run>.pdf`** — per‑sample/plate read‑count summaries.
- **`<run>_<marker>_<len>bp_NegativeControlOTUs.tsv`** — OTUs detected in negative controls.
- **OTU component reads** (zipped) when `--componentreads 1`.

For COI markers ≥ ~300 bp, results include **BOLD BIN** matches (BIN hit, % identity, BIN taxonomy).

---

## Build the image yourself

**Locally (single-arch, matches your own machine):**

From the `Docker/` directory (the build context):

```bash
cd Docker
docker build -t map .
```

This installs the full environment via `micromamba`, lays out `SCRIPTS/`, `Metabarcoding/`, and `REFS/`, installs Quarto and `iNEXT`, and downloads + unpacks teh reference library used for sequence correction of COI.

Run a locally‑built image by replacing `ghcr.io/cbg-innov/map:latest` with `map:latest` in any command above.

**Multi‑arch (published to ghcr.io):** the published `ghcr.io/cbg-innov/map:latest` image is a single multi‑arch manifest (linux/amd64 + linux/arm64), built by the `.github/workflows/docker-publish.yml` GitHub Actions workflow. amd64 builds natively on the GitHub runner and arm64 is built via QEMU emulation under `docker buildx`, so Docker automatically pulls the right layer whether you're on Apple Silicon or an x86 Linux/CI box.

---

## Reference library

The latest **BOLDdistilled** COI SINTAX reference set is downloaded and unpacked into `/MAP/REFS` **at first run time** and needs an internet connection, and may prolong the MAP demo runtime slightly. The COI correction reference set is provided **at first build time** — no manual download needed. With `docker compose`, `REFS` is mounted as a named volume (`reflib`) so it persists across container recreations and can be shared between containers. If you wish to change the reference library, make sure that the parameters.xlsx sheet reflects the new name and copy your vsearch reference file into your working directory. Make sure to also change the --refs flag while running via command line to include the file name, minus '.fasta', e.g.:
```bash
MAP_DATA="$(pwd)" docker compose -f compose.yaml \
  run --rm map \
  bash /MAP/SCRIPTS/MAP.sh \
  --fastq /data/*.fastq.gz \
  --params /data/my_parameters.xlsx \
  # Calls the reference file name 'my_refs' that is in the same working directory as the parameters and fastq files.
  --refs /data/my_refs \
  --wd /data
```

---


## Troubleshooting

- **No output on the host?** Make sure you mounted a volume to the output location (`-v "$(pwd)/MAP_results:/MAP/Metabarcoding/output"`), or use `--wd /data` with `-v "$(pwd)/my_run:/data"`. Alternatively retrieve results with `docker cp map:/MAP/Metabarcoding/output ./MAP_results`.
- **Disk space.** Large/deep runs (especially ONT) can generate many intermediate files; ensure adequate free disk on the Docker host.
- **`map` environment not active.** Use `bash -lc "…"`, or run inside `micromamba run -n map bash -c "…"`.

---

## Citation 

> Sean WJ Prosser, Nicholas W Bard, Ken A Thompson, Robin M Floyd, Sameer Padhye, Emine Ozsahin, Saeideh Jafarpour, and Paul DN Hebert. The Metabarcoding Analysis Pipeline (MAP): Simple, accurate, and flexible metabarcoding. <i> In prep. </i>


