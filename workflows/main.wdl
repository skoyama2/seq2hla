version 1.0

import "https://raw.githubusercontent.com/skoyama2/seq2hla/main/tasks/cram_to_mhcbam.wdl" as mhc

workflow main {
  input {
    File cram
    File crai
    File ref_fasta
    File ref_fasta_fai
    File ref_dict
  }

  call mhc.extract_mhc_with_mates {
    input:
      cram = cram,
      crai = crai,
      ref_fasta = ref_fasta,
      ref_fasta_fai = ref_fasta_fai,
      ref_dict = ref_dict
  }

  output {
    File mhc_bam = mhc.extract_mhc_with_mates.bam
    File mhc_bai = mhc.extract_mhc_with_mates.bai
  }
}
