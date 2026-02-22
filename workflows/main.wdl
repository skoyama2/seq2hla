version 1.0

import "../tasks/cram_to_mhcbam.wdl" as mhc

workflow main {
  input {
    File cram
    File crai
    File ref_fasta
    File ref_fasta_fai
    String samtools_docker = "quay.io/biocontainers/samtools:v1.19.2-1-deb_cv1"
  }

  call mhc.extract_mhc_with_mates {
    input:
      cram = cram,
      crai = crai,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      samtools_docker = samtools_docker
  }

  output {
    File mhc_bam = extract_mhc_with_mates.bam
    File mhc_bai = extract_mhc_with_mates.bai
  }
}
