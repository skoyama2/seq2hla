version 1.0

task extract_mhc_with_mates {
  input {
    File cram
    File crai

    # Reference Files
    File ref_fasta
    File ref_fasta_fai
    File ref_dict

    # Docker with samtools
    String samtools_docker = "biocontainers/samtools:v1.19.2-1-deb_cv1"
  }

  command <<<
    set -euo pipefail

    ln -sf "~{cram}"  sample.cram
    ln -sf "~{crai}"  sample.cram.crai

    ln -sf "~{ref_fasta}"     ref.fa
    ln -sf "~{ref_fasta_fai}" ref.fa.fai
    ln -sf "~{ref_dict}"      ref.dict

    SAMPLE=$(basename "sample.cram" | sed 's/\.cram$//')

    # Extract contig names in MHC region + alt contigs (chr6_*_alt / 6_*_alt / HLA*) with mapped read

    samtools idxstats sample.cram \
      | awk '$3>0 && ($1 ~ /^chr6_.*_alt$/ || $1 ~ /^6_.*_alt$/ || $1 ~ /^HLA/){print $1}' \
      > "${SAMPLE}.mhc_alt_contigs.txt"

    # Extract read names in MHC region (GRCh38 chr6:28510120-33480577) + alt contigs

    samtools view \
      -T ref.fa \
      -F 0x900 \
      sample.cram \
      chr6:28510120-33480577 \
      $(cat "${SAMPLE}.mhc_alt_contigs.txt" || true) \
        | cut -f 1 | sort -u > "${SAMPLE}.mhc_qnames.txt"

    # Extract unmapped read names

    samtools view \
      -T ref.fa \
      -F 0x900 \
      -f 12 \
      sample.cram \
        | cut -f 1 | sort -u > "${SAMPLE}.unmap_qnames.txt"

    # Merge read names in a file

    cat "${SAMPLE}.mhc_qnames.txt" "${SAMPLE}.unmap_qnames.txt" > "${SAMPLE}.qnames.txt"

    # Extract reads (+mates) and create .bam files

    samtools view \
      -T ref.fa \
      -b -F 0x900 -N "${SAMPLE}.qnames.txt" \
      sample.cram \
        > "${SAMPLE}.mhc_with_mates.bam"

    # Indexing

    samtools index "${SAMPLE}.mhc_with_mates.bam"

  >>>

  output {
    File bam = "${basename(cram, '.cram')}.mhc_with_mates.bam"
    File bai = "${basename(cram, '.cram')}.mhc_with_mates.bam.bai"
  }

  runtime {
    docker: samtools_docker
    cpu: 2
    memory: "8 GB"
    disks: "local-disk 200 HDD"
  }
}
