# NF_REFVAR_INDEL

## Description
A Nextflow pipeline for viral reference-based variant reporting and genome assembly.

## Workflow
- Trim reads with fastp (default minimum length 100bp)
- Consensus genome assembly using BBMap and iVar consensus
    - minimum coverage of 10
    - minimum base quality of 20
    - minimum frequency threshold of 0.6 
- Variant reporting using iVar variants
    - minimum coverage of 10
    - minimum base quality of 20
    - minimum frequency threshold of 0.01
- Generate summary and QC stats

## Requirements
Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

Install [`Docker`](https://docs.docker.com/engine/installation/)

## Usage (ex. RSVA)

### Run locally with Docker:
```bash
nextflow run main.nf \
    --input example_samplesheet.csv \
    --output example_output \
    --ref $(pwd)/assets/NC_038235.fa \
    --ref_index $(pwd)/assets/NC_038235.fa.fai \
    --ref_dict $(pwd)/assets/NC_038235.dict \
    --gff $(pwd)/assets/NC_038235.gff \
    --genomic_region "NC_038235.1:4688-5584" \
    --genomic_region_len 897 \
    -profile docker \
    -c your_nextflow_aws.config
```

### Run GitHub version with Docker:
```bash
nextflow run uwvirology-ngs/nf_refvar_indel -r realign_count_gaps -latest \
    --input example_samplesheet.csv \
    --output example_output \
    --ref $(pwd)/assets/NC_038235.fa \
    --ref_index $(pwd)/assets/NC_038235.fa.fai \
    --ref_dict $(pwd)/assets/NC_038235.dict \
    --gff $(pwd)/assets/NC_038235.gff \
    --genomic_region "NC_038235.1:4688-5584" \
    --genomic_region_len 897 \
    -profile docker
```

## Options

### Required Parameters
|Parameter|Explanation| Example Value |
|------|-----------|------|
| `--input` | samplesheet in csv format with fastq information | example_samplesheet.csv |
| `--output` | output directory (default: nf_output) | example_output |
| `--ref` | reference genome | assets/NC_045512.fa |
| `--ref_index` | corresponding index file | assets/NC_045512.fa.fai |
| `--gff` | general feature formal file | assets/NC_045512.gff |
| `--genomic_region` | genomic region of interest for variant calling | "NC_045512.2:21563-25384" | |
| `--genomic_region_len` | length of region of interest | 3822 |

### Optional Parameters
|Option|Explanation|
|------|-----------|
| `--run_name` | name for the summary tsv file (default: 'run') |
| `--skip_fastqc` | skip quality control using FastQC (default: false) |
| `--skip_fastp` | skip adapters and reads trimming using fastp (default: false) |
| `--min_trim_reads` | mininum number of trimmed reads required for downstream processes (default: 0) |
| `--trim_len` | minimum read length to keep (default:100) |
| `--save_trimmed_reads` | save trimmed fastq (default: false) |
| `--sample` | downsample fastq to a certain fraction or number of reads |
| `--min_mapped_reads` | minimum number of mapped reads for variants and consensus calling (default: 1000) |
| `--ivar_variants_t` | minimum frequency threshold to call consensus (default: 0.01) |
| `--ivar_variants_q` | minimum quality score threshold to call consensus (default: 20) |
| `--ivar_variants_m` | minimum depth to call consensus (default: 10) |
| `--save_mpileup` | save samtools mpileup used for iVar variants |
| `--ivar_consensus_t` | minimum frequency threshold to call consensus (default: 0.6) |
| `--ivar_consensus_q` | minimum quality score threshold to call consensus (default: 20) |
| `--ivar_consensus_m` | minimum depth to call consensus (default: 10) |
| `--trim_primers` | trim primers, requires bed file (default: false) |
| `--bed_file` | bed file containing primer coordinates (required if --trim_primers flag is set) |

## Notes
- Samplesheet example: `assets/samplesheet.csv`
- You can create a samplesheet using the bundled python script: `python bin/fastq_dir_samplesheet.py fastq_dir samplesheet_name.csv`
- Memory and CPU usage for pipeline processes can be adjusted in `conf/base.config`
- Process arguments can be adjusted in `conf/modules.config`
- If you are using Docker on Linux, check out these [post-installation steps](https://docs.docker.com/engine/install/linux-postinstall/) (especially cgroup swap limit capabilities support) for configuring Linux to work better with Docker.- By default, Docker has full access to full RAM and CPU resources of the host, but if you are using MacOS, go to Settings -> Resources in Docker Desktop to make sure enough resources are allocated to docker containers. 

## Contact
For bug reports please email aidants@uw.edu or raise an issue on Github.