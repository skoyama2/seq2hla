version 1.0

task extract_mhc_with_mates {
  input {
    File cram
    File crai

    # Reference Files
    File ref_fasta
    File ref_fasta_fai

    # Docker with samtools
    String samtools_docker = "quay.io/biocontainers/samtools:v1.19.2-1-deb_cv1"
  }

  command <<<
    set -euo pipefail

    SAMPLE=$(basename "~{cram}" | sed 's/\.cram$//')

    # Extract contig names in MHC region + alt contigs (chr6_*_alt / 6_*_alt / HLA*) with mapped read

    # Get MHC alt contig names from CRAM header (no index required)
    samtools view -H "~{cram}" \
      | awk '/^@SQ/ && ($2 ~ /SN:chr6_.*_alt$/ || $2 ~ /SN:6_.*_alt$/ || $2 ~ /SN:HLA/) {
          sub(/SN:/, "", $2); print $2
        }' \
      > "${SAMPLE}.mhc_alt_contigs.txt"

    # Extract read names in MHC region (GRCh38 chr6:28510120-33480577) + alt contigs

    samtools view \
      -T "~{ref_fasta}" \
      -X "~{crai}" \
      -F 0x900 \
      "~{cram}" \
      chr6:28510120-33480577 \
      $(cat "${SAMPLE}.mhc_alt_contigs.txt" || true) \
        | cut -f 1 | sort -u > "${SAMPLE}.mhc_qnames.txt"

    # Extract unmapped read names

    samtools view \
      -T "~{ref_fasta}" \
      -X "~{crai}" \
      -F 0x900 \
      -f 0xC \
      "~{cram}" \
        | cut -f 1 | sort -u > "${SAMPLE}.unmap_qnames.txt"

    # Merge read names in a file

    cat "${SAMPLE}.mhc_qnames.txt" "${SAMPLE}.unmap_qnames.txt" > "${SAMPLE}.qnames.txt"

    # Extract reads (+mates) and create .bam files

    samtools view \
      -T "~{ref_fasta}" \
      -X "~{crai}" \
      -b -F 0x900 -N "${SAMPLE}.qnames.txt" \
      "~{cram}" \
        > "${SAMPLE}.mhc_with_mates.bam"

    # Indexing

    samtools index "${SAMPLE}.mhc_with_mates.bam"

  >>>

  output {
    File bam = basename(cram, ".cram") + ".mhc_with_mates.bam"
    File bai = basename(cram, ".cram") + ".mhc_with_mates.bam.bai"
  }

  runtime {
    docker: samtools_docker
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 200 HDD"
  }
}
