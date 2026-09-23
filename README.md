# Snakemake pipeline for the calling and genotyping of genomic Structural Variants based long-read sequencing

**Authors: Thomas Brazier<sup>1</sup>, Lune Angevin<sup>1</sup>, Claire Lemaitre<sup>2</sup> and Claire Mérot<sup>1</sup>**

*Institutions: (1) UMR 6553 ECOBIO, University of Rennes (2) GenScale Team, IRISA-INRIA lab, University of Rennes*

This pipeline performs ensemble calling of Structural Variants (SV) from PacBio HiFi or Oxford Nanopore (ONT) long-read sequencing of one or more individuals. Each individual is processed on its own. SV calling is performed by a combination of two aligners (minimap2 + ngmlr) and four different tools: SVIM, Sniffles2, CuteSV2 and Debreak ([poster](images/Poster_PopGroup58_BRAZIER.pdf)). The eight independent callsets are then merged with JasmineSV and SVs are genotyped with SVJediGraph.

![The complete workflow of EvolSV.](images/workflow.png)

The combination of an ensemble of independent callsets and genotype likelihoods allows us to estimate the relative performance of each tool on a given dataset, and to weight these tools in order to get the best calls with a reliable proxy of their uncertainty.

![The ensemble method to estimate calls uncertainty. How the relative performance scores fo each tool are calculated.](images/sv_validation.png)


## Installing the pipeline

The following dependencies must be installed before launching the pipeline:
* conda
* snakemake >= 8.9.0
* mamba
* python
* pandas

First you have to clone the github repository where you wish to run analyses.

```
git clone https://github.com/thomasbrazier/evolsv.git
cd evolsv
```

Then, you can install dependencies in a `conda` environment.

```
conda env create -f workflow/envs/snakemake.yaml
conda activate snakemake
```


## New features coming


The current version is fully functional, yet I plan to implement in a near future new features to address more types of data (e.g., ONT) and improve computation times:
* Simplify the rule 'final_report' to generalize better to different options and datasets.
* Better cleanup and compression for temporary and output files (optimize storage).
* Improved documentation and tests.
* CI/CD.



All these changes are being implemented in the `caterpillar` branch.


## Using the pipeline

There are three config files to set up your analysis:
* `config/config.yaml`, where you specify the working directory and all the settings for the different tools.
* `config/samples.tsv`, a three columns data frame to specify the sample name, the SRA accessions and the reference genome accession. Use one row per SRA run. The `sample_name` column identifies the individual: rows with the same `sample_name` are the runs of one individual, and their reads are merged. The sheet can contain many individuals, but all rows must use the same genome accession. A `sample_name` can contain only letters, digits, `.`, `_` and `-`, because it becomes a directory name.
* a profile under `profiles/`, specifying resources for each rule (number of cpus, memory, runtime) and how jobs reach your job scheduler. One is shipped per scheduler:

| Profile | Scheduler | Executor |
| --- | --- | --- |
| `profiles/slurm/` | SLURM | `slurm` |
| `profiles/pbs/` | PBS Pro / OpenPBS (and Torque, with one line swapped) | `cluster-generic` |
| `profiles/sge/` | SGE / SoGE / UGE | `cluster-generic` |
| `profiles/ci/` | none -- dry runs and DAG checks only | -- |

Edit the queue/account settings at the top of the profile you use before the first run. Thread counts live in that profile's `set-threads` block, which lists every rule in the workflow; `.test/test_profiles.py` fails the build if a rule is missing from it or if it names a rule that no longer exists.



After setting up the required config files, you can launch the pipeline in two ways: either by running the `launcher.sh` script (where you can configure your settings if you are using a SLURM job scheduler), or simply by executing the following command in the terminal:

```
snakemake --snakefile ./workflow/Snakefile --configfile ./config/config.yaml --profile ./profiles/slurm --use-conda
```

Do not add `--cores N` to that command. It caps the threads of *every* rule at N even when jobs run on a cluster, so a `--cores 1` silently overrides the whole `set-threads` block. The profile sets `cores:` to a sensible upper bound instead.

Before running the analysis, you can build Conda environments:

```
snakemake --snakefile ./workflow/Snakefile --cores 1 --use-conda --conda-frontend conda --conda-create-envs-only
```

An example dataset with a single sample of the `Vanessa cardui` species can be run with the command:

```
snakemake --snakefile ./workflow/Snakefile --configfile ./config/config_test.yaml --profile ./profiles/slurm --use-conda
```

Alternatively, if you wish to run a custom analysis, with config files in a subdirectory:

```
species=Vanessa_cardui

snakemake -s workflow/Snakefile --configfile data/config/config_$species.yaml \
--use-conda --conda-frontend conda --profile ./profiles/slurm \
--config samples="data/config/samples_$species.tsv"
```


### Running on another job scheduler

The pipeline ships a profile per scheduler. Pick one, edit the site-specific settings at the top of its `config.yaml`, and point `--profile` at it.

**SLURM.** Set `slurm_partition` and `slurm_account` in `profiles/slurm/config.yaml`. The executor plugin is in `workflow/envs/snakemake.yaml`.

**PBS Pro / OpenPBS.** Snakemake has no dedicated PBS executor plugin, so `profiles/pbs/` submits through `cluster-generic`: the `qsub` line in the profile is the only place resources are translated, and a resource not named there is not requested.

```
conda env create -f workflow/envs/snakemake_pbs.yaml
mkdir -p logs/pbs                       # qsub fails if -o/-e do not exist
chmod +x profiles/pbs/pbs-status.py

snakemake --snakefile ./workflow/Snakefile --configfile ./config/config.yaml \
--profile ./profiles/pbs --use-conda
```

Set `pbs_queue` in the profile, and add `-A <project>` to `cluster-generic-submit-cmd` if your site bills projects. On **Torque**, swap the two `-l` lines for the commented alternative in the same file: Torque spells the request `-l nodes=1:ppn=N,mem=Nmb` where PBS Pro uses `-l select=1:ncpus=N:mem=Nmb`.

**SGE / SoGE / UGE.** Same mechanism, via `profiles/sge/`.

```
conda env create -f workflow/envs/snakemake_sge.yaml
mkdir -p logs/sge
chmod +x profiles/sge/sge-submit.sh profiles/sge/sge-status.py

snakemake --snakefile ./workflow/Snakefile --configfile ./config/config.yaml \
--profile ./profiles/sge --use-conda
```

Set `sge_queue`, and set `sge_pe` to a parallel environment that exists on your cluster (`qconf -spl`) -- a missing or wrong PE makes multi-threaded jobs fail at submission on strict clusters. Submission goes through `profiles/sge/sge-submit.sh` rather than a bare `qsub` line for two reasons: SGE enforces `h_vmem` per *slot*, so the wrapper divides the profile's total `mem_mb` by the thread count, and `qsub -terse` is required because `cluster-generic` reads the job id from the first line of stdout.

`profiles/sge/` deliberately does not use `snakemake-executor-plugin-sge`. That plugin is PyPI-only (absent from bioconda and conda-forge), and in 0.6.24 it reads the CPU count from a resource named `threads` rather than from the job's threads -- so every job is submitted single-slot -- and appends `-l tmem=`, a site-specific complex, to every submission.

> **Not yet run against real PBS or SGE hardware.** Both profiles build the DAG correctly and are derived from the executor plugin's own source, but no job has been submitted to a live PBS or SGE cluster. Run the `config_test.yaml` dataset end to end before trusting them for a real analysis, and expect the status scripts (`pbs-status.py`, `sge-status.py`) to be the parts most likely to need tuning against your site's job-retention settings.

### Sequencing technology

The pipeline supports PacBio HiFi and Oxford Nanopore reads. A single key in `config/config.yaml` selects which:

```yaml
sequencing_technology: hifi # hifi (PacBio HiFi/CCS) | ont (Oxford Nanopore)
```

Any other value is rejected before the run starts. The key fills in the technology-dependent parameters of every tool that has them (`TECH_PRESETS` in `workflow/rules/common.smk`):

| Parameter | `hifi` | `ont` |
| --- | --- | --- |
| `minimap_ax` (minimap2 `-x`) | `map-hifi` | `map-ont` |
| `ngmlr_preset` (NGMLR `--presets`) | `pacbio` | `ont` |
| `read_group_platform` (`@RG PL`) | `PACBIO` | `ONT` |
| `max_cluster_bias_INS` (cuteSV) | 1000 | 100 |
| `diff_ratio_merging_INS` (cuteSV) | 0.9 | 0.3 |
| `max_cluster_bias_DEL` (cuteSV) | 1000 | 100 |
| `diff_ratio_merging_DEL` (cuteSV) | 0.5 | 0.3 |
| `longqc_preset` ([LongQC](https://github.com/yfukasawa/LongQC) `-x`) | `pb-hifi` | `ont-ligation` |

The cuteSV values are the ones [recommended by its authors](https://github.com/tjiangHIT/cuteSV#recommendation-parameters) for each technology.

These are defaults, not overrides: setting any of those keys explicitly in your config file wins over the preset. They are shipped commented out in `config/config.yaml` so the preset applies by default. For instance, to run ONT with a non-default minimap2 preset:

```yaml
sequencing_technology: ont
minimap_ax: lr:hq # overrides the map-ont from the preset; everything else stays ONT
```

Two things worth knowing:

* **Adapter removal depends on the technology.** Each SRA run goes through one adapter tool before the runs are merged:
  * `hifi`: [HiFiAdapterFilt](https://github.com/sheinasim/HiFiAdapterFilt) (Sim et al. 2022, *BMC Genomics* 23:157). It removes whole reads that match PacBio adapters (`hifiadapterfilt_min_length`, `hifiadapterfilt_min_match`). The counts of removed reads are in `{wdir}/{sample}/hifiadapterfilt/{run}.stats`.
  * `ont`: [Porechop_ABI](https://github.com/bonsai-team/Porechop_ABI) ([*Bioinformatics Advances*](https://doi.org/10.1093/bioadv/vbac085)). It trims adapters from the read ends and splits reads with an adapter in the middle. With `porechop_ab_initio: true` (default) it also infers the adapters from the reads. This step samples reads at random and has no seed, so the inferred adapters can differ between two runs. They are recorded in `{wdir}/{sample}/logs/{run}.porechop_abi.log`.
* **Read filtering is not technology-dependent.** `chopper_quality: 10` is used for both, and suits HiFi as well as modern ONT chemistries (R10.4+). Read length thresholds (`chopper_minlength`, `chopper_maxlength`) are likewise unchanged by the technology.
* **Only the aligners and cuteSV have technology presets.** Sniffles2, SVIM, DeBreak and SVJedi-graph publish none upstream — Sniffles2 derives its thresholds from coverage, DeBreak accepts HiFi/CLR/ONT/mixed BAMs without a flag, and SVJedi-graph maps onto its variation graph with minigraph, which has no per-technology presets. An ONT run therefore uses the same settings as a HiFi run in those four tools. `min_sv_size`, `mapq` and `min_mapq` may deserve a second look on noisier data.


### Choosing the aligners

By default the reads are aligned twice, with minimap2 and with NGMLR, and the four SV callers run on both — the eight callsets the ensemble consensus is built from. A single aligner can be used instead, which roughly halves the cost of a run:

```yaml
aligners:
  - minimap2 # one or both of minimap2, ngmlr; the order does not matter
```

`--config aligners=ngmlr` works too. An unknown name, or an empty list, stops the run before it starts. Every per-alignment stage — the alignment, the four callers, the QC, the callability BED, the diagnostic plots — then runs once instead of twice, and the aligner-specific keys (`minimap_ax`, `ngmlr_preset`, `min-identity`) only matter for a selected aligner.

**With one aligner the ensemble loses the cross-aligner agreement it is built on.** Jasmine merges 4 callsets instead of 8, so `SUPP` is out of 4 and does not mean what a two-aligner `SUPP` means: the four callers share every systematic error of the alignment they all read, so their agreement overstates confidence. `callability/{genome}_callable.bed` likewise stops being the intersection of two alignments and becomes that aligner's own callable regions. The final report labels its tables from the callsets actually present, so it reports what was run — on a different scale from a two-aligner run.


### Starting from pre-aligned BAM files

If your reads are already aligned — for instance because you mapped them to call SNPs — you can skip the SRA download, the read QC and the alignments, which are by far the most expensive stages of the workflow. Set in `config/config.yaml`:

```yaml
start_from_bam: true
reference_fasta: "/path/to/the/reference/used/for/the/alignment.fna" # optional
```

and declare the input files in three extra columns of the sample sheet (see `config/samples_bam.tsv` for a template). The `sra` column is not used in this mode and can be left empty:

```
sample_name	sra	genome	bam_minimap2	bam_ngmlr	fastq
SAMEA8724893		GCA_947247005.1	/path/mm2.bam	/path/ngmlr.bam	/path/reads.fastq.gz
```

**One BAM per selected aligner is required.** With the default two, the ensemble method relies on eight independent callsets produced by four callers on two different alignments, and the merging step (JasmineSV/IRIS) is given both BAM files. Supplying the same alignment twice would make the same evidence count as two independent observations and would inflate both the consensus and the per-tool performance scores — declare one aligner in `aligners` instead (see [Choosing the aligners](#choosing-the-aligners)), which is the honest way to express having a single alignment. Only the selected aligners' columns are read, so a sheet may keep a `bam_ngmlr` path that an `aligners: [minimap2]` run ignores.

**The reads are still needed, but the `fastq` column is optional.** Genotyping with SVJedi-graph maps reads onto a variation graph, so it cannot work from a linear BAM. Give the read file(s) in the `fastq` column when you have them; several rows of the same individual are concatenated, as in the SRA mode. Declare exactly one `bam_minimap2` and one `bam_ngmlr` file per individual.

Left blank, the column makes rule `bam_to_fastq` recover the reads from that individual's BAM with `samtools fastq -F 0x900` — the first selected aligner's, so minimap2 whenever minimap2 is selected. The decision is per individual, so one sheet can mix both forms. Prefer the original FASTQ when it is available: extracted reads are **not** the reads that were sequenced, and the difference is recorded in each `bam_check.txt` report.

* Secondary (`0x100`) and supplementary (`0x800`) records are excluded. They must be: a long read whose alignment is split would otherwise re-enter the callers as several reads, inflating both the coverage they see and the read support SVJedi-graph counts.
* Unmapped records (`0x4`) are kept, but reads the aligner never wrote cannot be recovered. A BAM produced with minimap2 `--sam-hit-only` — which this pipeline's own `minimap2` rule uses — or filtered to mapped reads holds fewer reads than the original FASTQ.
* Reverse-strand reads are restored to their original orientation, but bases removed by *hard* clipping on a primary alignment are gone.
* The extracted FASTQ is a temporary file: it is deleted once every rule that consumes it has run, and re-extracted if you later rerun one of them.

Requirements on the BAM files, all checked before the run proceeds (see `workflow/scripts/check_bam_reference.py`, which writes a report to `{wdir}/{sample}/bam/{genome}_{aligner}_bam_check.txt`):

* coordinate-sorted, with a `.bai` index next to the BAM (a `.csi` index is not accepted);
* an `@RG` line whose `SM` tag equals `sample_name` in the sample sheet, because Samplot selects reads by sample id;
* contig names *and* lengths matching the reference. This is the check most likely to fire: the pipeline downloads the GenBank assembly from NCBI, so a BAM aligned against a RefSeq or UCSC copy of the same assembly will be rejected. Point `reference_fasta` at the exact FASTA you aligned to. The assembly metadata is then still downloaded from NCBI, unless you also give it locally (see [Local reference genome and metadata](#local-reference-genome-and-metadata)).

Caveats to be aware of when interpreting the results:

* **The `chopper` read filters are not applied.** In the SRA mode, `chopper_quality`, `chopper_minlength`, `chopper_maxlength`, `chopper_headcrop` and `chopper_tailcrop` decide which reads reach every caller and genotyper. In BAM mode the alignment is used as supplied and the reads passed to SVJedi-graph are unfiltered, so those config keys have no effect. Filter your reads before aligning if you need the equivalent behaviour. HiFiAdapterFilt and Porechop_ABI are not applied either. This holds for reads extracted from a BAM too: the `_filtered` in their filename only keeps the downstream rules identical between the two entry points.
* **Read-level QC (FastQC, NanoPlot, LongQC) is skipped.** Alignment QC is still produced in `mapping_QC/` and `callability/`, and the final report is unaffected.
* **`sequencing_technology` still matters.** The aligner presets and the `@RG PL` tag are unused in this mode, since the alignments are supplied, but the key still drives the cuteSV clustering parameters. Set it to the technology the BAM files were produced from.
* **All the BAM files of one individual are assumed to come from the same read set.** This is not enforced: minimap2 (run with `--sam-hit-only`) and ngmlr legitimately retain different numbers of records, so comparing read counts would raise false alarms. Aligning two different read sets would bias the relative performance scores of the tools.
* **Changing `aligners` between runs on the same `datadir`** re-triggers `jasmine` and everything downstream of it, because its input list changed. The dropped aligner's own per-aligner files are left on disk, unused; delete them if the disk matters, but nothing reads them.

The BAM files are symlinked, not copied, so no extra storage is used.


### Local reference genome and metadata

By default the reference genome (the `genome` column of the sample sheet) is downloaded from NCBI with its metadata. Two metadata files are used:

* `sequence_report.jsonl`: rule `autosomes_sexchromosomes` uses it to split autosomes and sex chromosomes in the final VCF files, and the final report shows it;
* `assembly_data_report.jsonl`: only the final report reads it.

Three config keys give local files instead. Only these combinations are accepted; any other stops the run before it starts:

| `reference_fasta` | `sequence_report` | `assembly_data_report` | Downloaded from NCBI |
| --- | --- | --- | --- |
| – | – | – | FASTA and metadata |
| set | – | – | metadata only |
| set | set | optional | nothing |

```yaml
reference_fasta: "/path/to/assembly.fna"
sequence_report: "/path/to/sequence_report.jsonl"
assembly_data_report: "/path/to/assembly_data_report.jsonl" # optional
```

This lets you run offline, or on an assembly that is not in NCBI. Without `assembly_data_report`, the final report skips its assembly section.

`sequence_report.jsonl` uses the NCBI JSON Lines format: one JSON object per sequence. The workflow reads these fields:

* `assemblyAccession`: one value for the whole file;
* `role`: only `assembled-molecule` sequences go into the final VCF files;
* `chrName`: `X`, `Y`, `Z` and `W` are sex chromosomes; the names in `scaffolds_to_exclude` are removed;
* `length`;
* the contig name, as in the FASTA: `genbankAccession` if `assemblyAccession` contains `GCA_`, else `refseqAccession`;
* `assignedMoleculeLocationType` (final report only).

**The report must match the FASTA.** Rule `check_reference_names` compares the report with the FASTA index in every mode, before the autosome split. The run stops if a field is missing, if no assembled molecule is found in the FASTA, or if a contig length differs from the FASTA. Without this check, a name mismatch (e.g. a RefSeq FASTA with a GenBank report) gives empty final VCF files and no error. Assembled molecules partly missing from the FASTA, and assemblies with no assembled molecule (scaffold level), only give a warning in `{wdir}/genome/{genome}_reference_names_check.txt`.


## Data directory setup

Project data can be stored in the current `evolsv` git directory. The place where is the `data/` directory must be specified in the parameter `workingdir` in the `config.yaml`. The default is `workingdir: data/` which assumes `data/` to be in te current `evolsv/` directory (see below). `data/` will not be tracked by `git`.

```
.
├── evolsv/
│   ├── config/
│   │   ├── config.yaml
├   |   |── samples.tsv
│   ├── data/
├   ├── profiles/
├   ├   ├── slurm/
│   ├   │   ├── config.yaml
├   ├   ├── pbs/
│   ├   │   ├── config.yaml
│   ├   │   ├── pbs-status.py
├   ├   ├── sge/
│   ├   │   ├── config.yaml
│   ├   │   ├── sge-submit.sh
│   ├   │   ├── sge-status.py
├   ├   ├── ci/
│   ├   │   ├── config.yaml
│   ├── workflow/
```



## Output files

Results are written to one directory per individual, `{wdir}/{sample}/`, where `{wdir}` is `datadir` followed by the genome accession and `{sample}` is the `sample_name` of the individual. The reference genome and the files computed from it (`genome/`, `genmap/`, `mappability/`) are shared by all individuals and computed once:

```
{wdir}/
├── genome/            shared reference
├── genmap/            shared mappability
├── mappability/       shared mappability
├── SAMEA8724893/      one individual
│   ├── bam/ calling/ genotype/ ...
│   ├── {genome}_final.vcf.gz
│   └── {genome}_finalQC.html
└── SAMEA0000002/      another individual
    └── ...
```

The main output file is a VCF file containing the list of SVs of one individual, named `{wdir}/{sample}/{genome}_final.vcf.gz`. It is the result of merging the eigth SV catalogues generated for that individual. Additionnally, a `{wdir}/{sample}/{genome}_final.tsv` and a `{wdir}/{sample}/{genome}_final_light.vcf.gz` files are produced. They contain the same set of SV calls, but they are designed to be processed more easily than the full vcf. `{wdir}/{sample}/{genome}_final.tsv` is a data frame without sequences for an easy import in R for data analysis. `{wdir}/{sample}/{genome}_final_light.vcf.gz` is a lighter vcf without DNA sequences in REF/ALT and INFO fields (DNA sequences can be very large with structural variation).


We also produce an automatic report to assess the quality and empirical performance of the workflow for each individual. Please check `{wdir}/{sample}/{genome}_finalQC.html` for details.

**The results of different individuals are not merged.** Each individual is called, merged and genotyped alone. SV IDs, breakpoints and alleles are therefore not harmonised across individuals: the same SV can have different IDs and slightly different positions in two individuals. A population-level analysis needs a separate step that merges the per-individual VCFs (for example with JasmineSV) and genotypes all individuals at the same sites.




## Known issues


* **high impact/low probability** issue with **SVjdedi-graph**. The current version fails with some SVs close to chromosome boundaries. The fix is to filter these variants with the `filter_variant_positions` option in `config/config.yaml`, until a new corrected version of SVjedi-graph is released.