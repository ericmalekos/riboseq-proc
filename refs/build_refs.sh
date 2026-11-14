#!/usr/bin/env bash
# build_refs.sh: the reference files riboseq-proc needs, built from one definition of "ncRNA" so that the read
# depletion before STAR and the locus filter after STAR agree.
#
#   contaminant bowtie2 index (reads that align to it are removed before STAR):
#     exons of gene_type rRNA, Mt_rRNA, rRNA_pseudogene, Mt_tRNA and miRNA in the annotation GTF;
#     every record of the GENCODE tRNA GTF (cytoplasmic tRNAs are not in the main GENCODE annotation), except
#     records whose gene_type matches --trna-skip (for example Pseudo_tRNA: in mouse most tRNAscan-SE predictions
#     are SINE-derived pseudo-tRNAs, 22,983 of the 26,277 records of GENCODE vM38);
#     whole rDNA sequences given as FASTA (the 45S repeat and 5S rRNA are not annotated loci of a primary assembly,
#     and they are the largest contaminant of a Ribo-seq library).
#   ncrna.bed (reads over these loci are removed after STAR): the genomic loci above plus exons of snoRNA, snRNA,
#     scaRNA, misc_RNA, vault_RNA, sRNA, scRNA and ribozyme genes.
#   chroms.txt: the chromosomes kept after STAR (sequence names of the genome that match --chrom-regex).
#   star/ (only with --star-gtf): STAR genome index.
#
#   build_refs.sh --gtf ANNOTATION.gtf[.gz] --genome GENOME.fa --out REFDIR
#                 [--trna-gtf tRNAs.gtf[.gz]] [--trna-skip REGEX] [--rdna RDNA.fa]...
#                 [--chrom-regex '^chr([0-9]+|X|Y)$']
#                 [--contam-prefix PREFIX] [--ncrna-bed FILE] [--star-gtf GTF --overhang N] [--threads N]
#
# REFDIR receives contam/contam.fa, contam/contam.*.bt2, ncrna_deplete.bed, ncrna_postfilter.bed, ncrna.bed,
# chroms.txt and README.txt. --contam-prefix and --ncrna-bed write the index and the union BED elsewhere (for a
# shared index directory). Tools from PATH, or BEDTOOLS, BOWTIE2_BUILD, SAMTOOLS, STAR_BIN.
set -euo pipefail
BEDTOOLS=${BEDTOOLS:-bedtools}; BOWTIE2_BUILD=${BOWTIE2_BUILD:-bowtie2-build}; SAMTOOLS=${SAMTOOLS:-samtools}
STAR_BIN=${STAR_BIN:-STAR}
DEPLETE='rRNA|Mt_rRNA|rRNA_pseudogene|Mt_tRNA|miRNA'
POSTFILT='snoRNA|snRNA|scaRNA|misc_RNA|vault_RNA|sRNA|scRNA|ribozyme'
GTF=""; GENOME=""; OUT=""; TRNA=""; TRNA_SKIP=""; RDNA=(); CHROM_RE='^chr([0-9]+|X|Y)$'; PREFIX=""; NCBED=""; STAR_GTF=""; OH=""
T=${THREADS:-4}
die() { echo "ABORT: $*" >&2; exit 1; }
while (( $# )); do
  case "$1" in
    --gtf) GTF=$2; shift 2 ;;
    --genome) GENOME=$2; shift 2 ;;
    --out) OUT=$2; shift 2 ;;
    --trna-gtf) TRNA=$2; shift 2 ;;
    --trna-skip) TRNA_SKIP=$2; shift 2 ;;
    --rdna) RDNA+=("$2"); shift 2 ;;
    --chrom-regex) CHROM_RE=$2; shift 2 ;;
    --contam-prefix) PREFIX=$2; shift 2 ;;
    --ncrna-bed) NCBED=$2; shift 2 ;;
    --star-gtf) STAR_GTF=$2; shift 2 ;;
    --overhang) OH=$2; shift 2 ;;
    --threads) T=$2; shift 2 ;;
    *) die "unknown option $1" ;;
  esac
done
[[ -n "$GTF" && -n "$GENOME" && -n "$OUT" ]] || die "--gtf, --genome and --out are required"
for f in "$GTF" "$GENOME" ${TRNA:+"$TRNA"} ${RDNA[@]+"${RDNA[@]}"} ${STAR_GTF:+"$STAR_GTF"}; do [[ -s "$f" ]] || die "missing $f"; done
[[ "$GENOME" != *.gz ]] || die "the genome FASTA must be uncompressed"
[[ -z "$STAR_GTF" || "$OH" =~ ^[0-9]+$ ]] || die "--star-gtf needs --overhang N (read length minus 1)"
mkdir -p "$OUT"
PREFIX=${PREFIX:-$OUT/contam/contam}; NCBED=${NCBED:-$OUT/ncrna.bed}
mkdir -p "$(dirname "$PREFIX")" "$(dirname "$NCBED")"
readf() { case "$1" in *.gz) gzip -dc -- "$1" ;; *) cat -- "$1" ;; esac; }
exons() {   # gene_type pattern -> BED4 of the exons of those genes (name = gene_id)
  readf "$GTF" | awk -F'\t' -v pat="$1" '$3 == "exon" && $9 ~ ("gene_type \"(" pat ")\"") {
    match($9, /gene_id "[^"]+"/); printf "%s\t%d\t%d\t%s\n", $1, $4 - 1, $5, substr($9, RSTART + 9, RLENGTH - 10) }'
}
merged() { sort -k1,1 -k2,2n | "$BEDTOOLS" merge -c 4 -o distinct -i -; }
echo "[$(date '+%F %T')] build_refs on $(hostname): gtf=$GTF genome=$GENOME out=$OUT"

# 1. intervals
{ exons "$DEPLETE"
  if [[ -n "$TRNA" ]]; then
    readf "$TRNA" | awk -F'\t' -v skip="$TRNA_SKIP" '!/^#/ && ($3 == "exon" || $3 == "tRNA") {
      if (skip != "" && match($9, /gene_type "[^"]+"/) && substr($9, RSTART + 11, RLENGTH - 12) ~ ("^(" skip ")$")) { ns++; next }
      g = "tRNA"; if (match($9, /gene_id "[^"]+"/)) g = substr($9, RSTART + 9, RLENGTH - 10)
      nk++; printf "%s\t%d\t%d\t%s\n", $1, $4 - 1, $5, g }
      END { printf "tRNA GTF records used %d, skipped %d\n", nk, ns > "/dev/stderr" }'
  fi
} | merged > "$OUT/ncrna_deplete.bed"
exons "$POSTFILT" | merged > "$OUT/ncrna_postfilter.bed"
cat "$OUT/ncrna_deplete.bed" "$OUT/ncrna_postfilter.bed" | cut -f1-3 | sort -k1,1 -k2,2n | "$BEDTOOLS" merge -i - > "$NCBED.tmp"
mv "$NCBED.tmp" "$NCBED"
n_dep=$(wc -l < "$OUT/ncrna_deplete.bed"); n_post=$(wc -l < "$OUT/ncrna_postfilter.bed"); n_all=$(wc -l < "$NCBED")
bp_all=$(awk '{s += $3 - $2} END {print s + 0}' "$NCBED")
echo "deplete intervals $n_dep, postfilter intervals $n_post, union $n_all intervals ($bp_all bp)"
(( n_dep > 100 )) || die "only $n_dep depletion intervals: is gene_type present in $GTF?"
[[ -n "$TRNA" ]] || echo "WARNING: no --trna-gtf: cytoplasmic tRNAs are not depleted"
(( ${#RDNA[@]} > 0 )) || echo "WARNING: no --rdna: the 45S rDNA repeat is not depleted"

# 2. contaminant FASTA and bowtie2 index
"$BEDTOOLS" getfasta -fi "$GENOME" -bed "$OUT/ncrna_deplete.bed" -fo "$PREFIX.loci.fa"
cat "$PREFIX.loci.fa" ${RDNA[@]+"${RDNA[@]}"} > "$PREFIX.fa"
rm -f "$PREFIX.loci.fa"
n_seq=$(grep -c '^>' "$PREFIX.fa"); bp_seq=$(grep -v '^>' "$PREFIX.fa" | tr -d '\n' | wc -c)
[[ "$(grep '^>' "$PREFIX.fa" | cut -d' ' -f1 | sort | uniq -d | wc -l)" -eq 0 ]] || die "duplicate sequence names in $PREFIX.fa"
echo "contaminant FASTA: $n_seq sequences, $bp_seq bp"
"$BOWTIE2_BUILD" --threads "$T" -q "$PREFIX.fa" "$PREFIX"
[[ -s "$PREFIX.1.bt2" ]] || die "bowtie2-build wrote no index"

# 3. chromosomes kept after STAR
[[ -s "$GENOME.fai" ]] || "$SAMTOOLS" faidx "$GENOME"
cut -f1 "$GENOME.fai" | grep -E "$CHROM_RE" > "$OUT/chroms.txt" || die "no sequence name of $GENOME matches $CHROM_RE"
echo "chromosomes kept: $(tr '\n' ' ' < "$OUT/chroms.txt")"

# 4. STAR index (optional)
if [[ -n "$STAR_GTF" ]]; then
  mkdir -p "$OUT/star"
  "$STAR_BIN" --runMode genomeGenerate --runThreadN "$T" --genomeDir "$OUT/star" --genomeFastaFiles "$GENOME" \
    --sjdbGTFfile "$STAR_GTF" --sjdbOverhang "$OH" --outTmpDir "$OUT/star_tmp" --outFileNamePrefix "$OUT/star/"
  [[ -s "$OUT/star/SA" ]] || die "STAR wrote no index"
fi

readme() {
  cat <<TXT
riboseq-proc references, built $(date -Iseconds) by refs/build_refs.sh
annotation  $GTF (md5 $(md5sum "$GTF" | cut -d' ' -f1))
tRNA GTF    ${TRNA:-none}${TRNA:+ (md5 $(md5sum "$TRNA" | cut -d' ' -f1))}${TRNA_SKIP:+; records with gene_type $TRNA_SKIP skipped}
rDNA FASTA  $(if (( ${#RDNA[@]} )); then for f in "${RDNA[@]}"; do printf '%s (md5 %s, %s) ' "$f" "$(md5sum "$f" | cut -d' ' -f1)" "$(head -n 1 "$f" | cut -c2-60)"; done; else echo none; fi)
genome      $GENOME (.fai md5 $(md5sum "$GENOME.fai" | cut -d' ' -f1))
contaminant index $PREFIX: $n_seq sequences, $bp_seq bp = merged exons of gene_type $DEPLETE${TRNA:+, the tRNA GTF records}${RDNA[@]+, the rDNA FASTA}
  use: bowtie2 --very-fast-local -x $PREFIX -U reads.fq --un-gz unaligned.fq.gz -S /dev/null
ncRNA loci  $NCBED: $n_all merged intervals, $bp_all bp = the $n_dep depletion intervals plus $n_post intervals of gene_type $POSTFILT
  use: bedtools intersect -v -split -abam in.bam -b $NCBED > kept.bam   (aligned blocks only; samtools view -L
       $NCBED -U kept.bam also removes reads spliced across a locus)
chromosomes $OUT/chroms.txt: $(wc -l < "$OUT/chroms.txt") names matching $CHROM_RE
$( [[ -z "$STAR_GTF" ]] || echo "STAR index  $OUT/star: $STAR_GTF, --sjdbOverhang $OH, $("$STAR_BIN" --version)" )
tools       $("$BEDTOOLS" --version), bowtie2-build $("$BOWTIE2_BUILD" --version | head -n 1 | sed 's/.*version //')
TXT
}
readme > "$OUT/README.txt"
[[ "$PREFIX" == "$OUT/contam/contam" ]] || cp "$OUT/README.txt" "$PREFIX.README"
cat "$OUT/README.txt"
echo "[$(date '+%F %T')] done"
