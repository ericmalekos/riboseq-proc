# riboseq-proc

Ribo-seq reads to uniquely mapped, ncRNA-filtered BAMs with a 3-nt periodicity check. One container image, one
bash entrypoint, one sample sheet. The output BAMs are the input for ORF callers (Ribo-TISH, PRICE, RiboCode,
RiboTIE); ORF calling itself is not part of this tool.

```
source FASTQ -> [md5] -> [UMI extract] -> adapter trim -> [reverse complement] -> length window
             -> contaminant depletion (bowtie2)
             -> STAR (unique mappers, end to end) -> canonical chromosomes -> ncRNA locus filter -> [UMI dedup]
             -> ribotish quality + periodicity gate -> BAM, QC tables, done marker
```

## Image

`ghcr.io/ericmalekos/riboseq-proc:<VERSION>` (the tag is the content of `VERSION`; there is no `latest`, and a
published tag is never rebuilt). Every conda package is pinned by URL and md5 in `locks/`.

| tool | version | used for |
|---|---|---|
| STAR | 2.7.11b | alignment |
| cutadapt | 5.2 | adapter and length trimming |
| bowtie2 | 2.5.5 | contaminant depletion |
| samtools | 1.23.1 | filtering, indexing |
| UMI-tools | 1.1.6 | UMI extraction and deduplication |
| Ribo-TISH | 0.2.8 | read-length periodicity and P-site offsets |
| bedtools | 2.31.1 | ncRNA locus filter, reference building |
| rclone | 1.75.2 | reading from and writing to object storage |

```bash
docker run --rm ghcr.io/ericmalekos/riboseq-proc:0.1.2 riboseq-proc version
singularity exec docker://ghcr.io/ericmalekos/riboseq-proc:0.1.2 riboseq-proc version
```

The image runs as a non-root user and has no entrypoint: pass the command.

## Commands

```
riboseq-proc run     --sheet SHEET --run RUN [--stage prep|align|post|all]
riboseq-proc dataset --sheet SHEET --dataset DATASET [--par N]
riboseq-proc version
```

`run` processes one row of the sheet. `dataset` processes every row of one dataset, `PAR` runs at a time, each
from prep to post with `THREADS / PAR` threads. The STAR index is loaded into shared memory once
(`STAR --genomeLoad LoadAndExit`), every STAR run attaches to it (`LoadAndKeep`), and it is removed when the
command ends, also after a failure or a SIGTERM. With `SHARE_GENOME=0` each STAR run loads the index itself and
only one STAR runs at a time. A run that has a done marker is skipped, so a restarted job continues where the
last one stopped.

Exit codes of `run`: 0 done, 2 stopped by a data gate (a property of the reads, not a failure of the tool),
1 anything else. `dataset` exits 1 when at least one run failed for a reason other than a gate.

## Settings (environment variables)

| variable | default | meaning |
|---|---|---|
| `REF` | `/ref` | reference bundle (layout below); its parts can be set one by one with `STARIDX`, `CONTAM`, `NCRNA_BED`, `QC_GTF`, `CHROMS` |
| `WORK` | `/work` in the image | scratch; node-local disk is best |
| `OUT` | `$PWD/results` | results |
| `OUT_REMOTE` | empty | rclone remote path; when set, QC tables, BAMs and done markers are copied there and compared by checksum, and the local BAM is removed after it verified |
| `THREADS` | `$NRP_CPUS`, `$SLURM_CPUS_PER_TASK` or 4 | always set it inside containers: `nproc` reports the whole node |
| `PAR` | 1 | runs of a dataset processed at the same time (also `--par`) |
| `KEEP_CLEAN` | 0 | 1 keeps the depleted FASTQ in `OUT/<dataset>/fastq_clean/` |
| `BAMSORT_RAM` | 4000000000 | STAR `--limitBAMsortRAM`, per STAR run; a run that needs more stops with STAR's message, which names the value to set |
| `SHARE_GENOME` | 1 | `dataset` only. 1: the STAR index is held once in SysV shared memory for all runs; 0: no shared memory, one STAR at a time |
| `NCRNA_FILTER` | `blocks` | `blocks`: a read is removed when one of its aligned blocks overlaps an ncRNA locus; `span`: when its alignment span does, introns included (the behaviour of `samtools view -L`) |
| `MAX_READS` | 0 | test mode: process only the first N reads of each source; the md5 is then not checked, and the summaries and the done marker say so |

Memory: the STAR index (about 27 GB for a mouse or human genome) is held once, whatever `PAR` is. Each run adds its
BAM sort buffer (`BAMSORT_RAM`) and about 2 GB for cutadapt, bowtie2 and the post stage, so plan for the index plus
`PAR` x (`BAMSORT_RAM` + 2 GB). With `SHARE_GENOME=1` the index stays in memory during the prep and post stages
too, which keeps the memory use of a job steady. Shared memory needs `kernel.shmmax` and `kernel.shmall` at least
as large as the index (the Linux default is unlimited), and in a container it counts toward the memory limit.
A SIGKILL leaves the segment behind: `ipcs -m` lists it, `STAR --genomeLoad Remove --genomeDir DIR` or
`ipcrm -m ID` frees it. Only one `dataset` command per index and host should use `SHARE_GENOME=1` at a time, since
the first to finish removes the segment.

## Sample sheet

Tab separated, one row per run, this exact header, `-` for "none" (no empty fields):

| column | content |
|---|---|
| `dataset` | study accession; results are grouped by it |
| `run` | run accession, unique in the sheet |
| `sample`, `condition` | free labels, copied to the summaries |
| `source_type` | `fastq` (a path the container can read) or `rclone` (a remote path read with `rclone cat`) |
| `source` | the FASTQ, gzipped or plain; for paired-end libraries give read 1 |
| `md5` | md5 of the source file; checked on the same stream that is processed, `-` skips the check |
| `adapter` | 3' adapter sequence, or `-` for reads deposited already trimmed |
| `adapter_mode` | `require` drops reads without the adapter (a footprint is shorter than the read, so it must reach the adapter); `optional` keeps them |
| `orient` | `fwd` when the read is the footprint, `rev` when it is its reverse complement (for example read 2 of a pair when read 1 is unusable); `rev` reads are reverse-complemented after adapter trimming |
| `umi` | `-`, `raw:<pattern>` (extract before adapter trimming) or `trimmed:<pattern>` (extract from the insert after adapter trimming). A pattern with named groups such as `^(?P<umi_1>.{2}).+(?P<umi_2>.{5})$` is a UMI-tools regex, any other a 5' string pattern such as `NNNN` |
| `pass2` | extra cutadapt options for the second pass, for example `-u 3` or `-q 20`; `-` for none |
| `minlen`, `maxlen` | footprint length window after trimming (20 and 40 for standard RNase I footprints) |
| `expect` | range the modal length of the depleted reads should fall in, for example `28-32` (gate 2) |
| `umi_dedup` | `yes` or `no` (`umi_tools dedup` on the final BAM) |
| `rtlens` | read lengths tested by `ribotish quality`, for example `25,35` |
| `minpct` | minimum share of reads in periodic lengths (gate 3), 50 for RNase I libraries |
| `notes` | free text |

`examples/sheet.example.tsv` has four typical rows.

## Stages and gates

**prep.** The first 1 M raw reads are checked for the 3' adapter named in the sheet, for the Illumina TruSeq
adapter `AGATCGGAAGAGC`, for the modal read length and for the SRA Lite quality signature. The run stops when
more than half of the reads carry one quality character (`?` or `$`) throughout, or when the sheet names no adapter
but more than 5% of reads contain the TruSeq adapter. Then one stream runs from the source to the depleted reads:
md5 of the source bytes, optional UMI extraction, `cutadapt -a ADAPTER -e 0.1 -O 3 -m 2`
(with `--discard-untrimmed` in `require` mode), reverse complement for `orient rev`, a second cutadapt pass with
the options of `pass2` and the length window, and
`bowtie2 --sensitive-local` against the contaminant index, keeping the reads that do not align
(0.1.1 and earlier used `--very-fast-local`, which lets rRNA fragments with one to three extra 5' bases through;
some of them then align uniquely to unannotated genomic rRNA copies and reach the final BAM).
Gate 2 compares the modal length of the depleted reads with `expect`: PASS inside, WARN outside, FAIL (run stopped)
when the mode sits at the upper length cap or no read is left. The reads before depletion are not used for the
gate: in an rRNA-rich library their mode is often one rRNA fragment.

**align.** `STAR --alignEndsType EndToEnd --outFilterMultimapNmax 1 --outFilterMismatchNmax 1
--outFilterIntronMotifs RemoveNoncanonicalUnannotated`, sorted BAM. Only uniquely mapped reads are kept at the
aligner: the common ORF callers do not remove multi-mappers themselves.

**post.** Reads on the chromosomes of `chroms.txt` are kept, reads that overlap `ncrna.bed` are removed, UMIs
are deduplicated when asked, the reads on annotated CDSs are counted by strand (`strand` in the summary: `sense`
at 80% or more, `antisense` at 20% or less, which means the row needs `orient rev`, else `mixed`), and
`ribotish quality` measures the frame distribution per read length over annotated CDSs. Gate 3 (`bin/periodicity_gate.py`) gives PASS, WARN or FAIL:
FAIL when no read length is periodic or when the periodic lengths hold less than `minpct` of the reads; WARN
below 5 M reads in periodic lengths or for implausible P-site offsets. The BAM is delivered whatever the verdict;
the verdict is recorded in the summaries and in the done marker.

The ncRNA filter has two modes (`NCRNA_FILTER`). `blocks`, the default, is `bedtools intersect -v -split`: only
the aligned blocks of a read are compared with the loci (with `-nonamecheck`, since the BED may name sequences
the BAM header lacks), so a read spliced across an intron that holds a snoRNA, a miRNA or a tRNA is kept. `span` is `samtools view -L ... -U`, which compares the whole alignment span
and removes such reads although none of their bases lies on the locus.

## Outputs

```
OUT/<dataset>/
  bam/<run>.bam, .bam.bai          final BAM (removed locally after a verified upload when OUT_REMOTE is set)
  qc/<run>.raw_checks.tsv          adapter, length and quality checks on the first 1 M reads
  qc/<run>.cutadapt_*.log, .bowtie2_contam.log, .umi_extract.log, .umi_dedup.log
  qc/<run>.trim.lengths.tsv, .clean.lengths.tsv
  qc/<run>.prep_summary.tsv        read counts through trimming and depletion, gate 2
  qc/<run>.STAR.Log.final.out, .SJ.out.tab (.STAR.Log.out when STAR failed)
  qc/<run>.qual.txt, .qual.pdf, .para.py   Ribo-TISH quality output and P-site offsets
  qc/<run>.gate3.txt, .gate3.tsv   periodicity verdict
  qc/<run>.align_summary.tsv       read counts through alignment and filtering, gate 3, CDS strand
  status/<run>.tsv                 one line per stage: stage, state, detail, time
  done/<run>.json                  written last; holds the md5 of the BAM (used to verify the upload)
  logs/<run>.log                   (dataset command)
  prep_summary.tsv, align_summary.tsv, run_states.tsv   (dataset command) all runs of the dataset
```

Two runs of the same input give BAMs with the same records but not the same md5: STAR with more than one thread
does not fix the order of reads that start at the same position, and the header holds the command lines. Compare
`samtools view FILE | sort | md5sum`, or the read counts of `align_summary.tsv`.

## References

`refs/build_refs.sh` builds the contaminant index, the ncRNA BED and the chromosome list from one definition, so
that the depletion before STAR and the locus filter after STAR agree:

- contaminant index: exons of `gene_type` rRNA, Mt_rRNA, rRNA_pseudogene, Mt_tRNA and miRNA in the annotation;
  the records of the GENCODE tRNA GTF (cytoplasmic tRNAs are not in the main GENCODE GTF; `--trna-skip` leaves
  out records by `gene_type`, for example `Pseudo_tRNA`); whole rDNA sequences given as FASTA (the 45S repeat is
  not part of a primary assembly and is the largest contaminant of a Ribo-seq library);
- `ncrna.bed`: those loci plus exons of snoRNA, snRNA, scaRNA, misc_RNA, vault_RNA, sRNA, scRNA and ribozyme
  genes;
- `chroms.txt`: sequence names of the genome that match `--chrom-regex` (default `^chr([0-9]+|X|Y)$`);
- optionally a STAR index (`--star-gtf GTF --overhang N`).

```bash
build_refs.sh --gtf gencode.vM38.annotation.gtf --genome GRCm39.primary_assembly.genome.fa \
    --trna-gtf gencode.vM38.tRNAs.gtf.gz --rdna BK000964.3.fa --rdna NR_030686.1.fa --out ref
```

Bundle layout expected under `REF`:

```
ref/star/              STAR index (built with the annotation used for ORF calling; --sjdbOverhang about read length minus 1)
ref/contam/contam.*    bowtie2 index
ref/ncrna.bed
ref/chroms.txt
ref/qc.gtf             annotation with CDS records for the periodicity check (for example the GENCODE GTF)
```

## Running a dataset

```bash
docker run --rm -v /data/ref:/ref:ro -v /data/fastq:/data/fastq:ro -v /scratch/out:/out \
    -e THREADS=16 -e OUT=/out ghcr.io/ericmalekos/riboseq-proc:0.1.2 \
    riboseq-proc dataset --sheet /out/sheet.tsv --dataset GSE00001 --par 2
```

With object storage, give the sources as `rclone` rows, set `OUT_REMOTE` and pass the rclone configuration through
the environment (`RCLONE_CONFIG_<REMOTE>_...` variables from a secret, never written into the job specification).

## License

The scripts of this repository are under the MIT license (`LICENSE`). The tools installed in the image keep their own
licenses (among them STAR and Ribo-TISH: GPL-3.0; cutadapt and UMI-tools: MIT).
