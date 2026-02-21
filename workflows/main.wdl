version 1.0

import "https://raw.githubusercontent.com/skoyama2/seq2hla/main/tasks/cram_to_mhcbam.wdl" as mhc

workflow main {
  input {
    File cram
    File crai
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
    String samtools_docker = "biocontainers/samtools:v1.7.0_cv4"
  }

  call mhc.extract_mhc_with_mates {
    input:
      cram = cram,
      crai = crai,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict,
      samtools_docker = samtools_docker
  }

  output {
    File mhc_bam = extract_mhc_with_mates.bam
    File mhc_bai = extract_mhc_with_mates.bai
  }
}
