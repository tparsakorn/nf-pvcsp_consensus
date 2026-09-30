#!/usr/bin/env bash
# Usage: collect_consensus.sh <barcode_id> <ngspeciesid_outdir> <filtered.fastq>
# Renames each medaka consensus to <barcode>_consensus_cluster<N>.fasta (N = 1 is the most
# abundant cluster) and writes read-count tables:
#   <barcode>_cluster_read_counts.tsv : clusters that received a consensus
#   <barcode>_cluster_sizes.tsv       : ALL clusters found (incl. those under --abundance_ratio)
set -euo pipefail
bc="$1"; outdir="$2"; fastq="$3"

total=$(( $(wc -l < "$fastq") / 4 ))

# all clusters (used to spot low-abundance haplotypes that missed the cutoff)
{
    printf 'ngspeciesid_cl_id\treads\tpct_of_filtered_reads\n'
    cut -f1 "$outdir/final_clusters.tsv" | sort | uniq -c | sort -k1,1nr \
        | awk -v t="$total" '{printf "%s\t%d\t%.2f\n", $2, $1, 100*$1/t}'
} > "${bc}_cluster_sizes.tsv"

# clusters with a consensus, ordered by read support
shopt -s nullglob
files=( "$outdir"/medaka_cl_id_*/consensus.fasta )
if [ ${#files[@]} -eq 0 ]; then
    echo "ERROR: no consensus produced for $bc (no cluster passed --abundance_ratio?)" >&2
    exit 1
fi

for f in "${files[@]}"; do
    hdr=$(head -1 "$f")
    reads=$(printf '%s' "$hdr" | sed -E 's/.*total_supporting_reads_([0-9]+).*/\1/')
    clid=$(printf '%s' "$hdr"  | sed -E 's/.*consensus_cl_id_([0-9]+)_.*/\1/')
    printf '%s\t%s\t%s\n' "$reads" "$clid" "$f"
done | sort -k1,1nr > .ranked.tmp

printf 'cluster\tngspeciesid_cl_id\treads\tpct_of_filtered_reads\n' > "${bc}_cluster_read_counts.tsv"
n=0
while IFS=$'\t' read -r reads clid f; do
    n=$((n+1))
    {
        printf '>%s_cluster%d cl_id=%s reads=%s\n' "$bc" "$n" "$clid" "$reads"
        tail -n +2 "$f"
    } > "${bc}_consensus_cluster${n}.fasta"
    awk -v c="cluster${n}" -v i="$clid" -v r="$reads" -v t="$total" \
        'BEGIN{printf "%s\t%s\t%d\t%.2f\n", c, i, r, 100*r/t}' >> "${bc}_cluster_read_counts.tsv"
done < .ranked.tmp
rm -f .ranked.tmp
