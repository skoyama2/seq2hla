# seq2hla — Code Reference

## Overview

WDL-based bioinformatics pipeline that extracts HLA/MHC reads from whole-genome CRAM files.
Designed to run on Terra (Cromwell backend) with Google Cloud Storage inputs.

---

## Directory Structure

```
seq2hla/
├── workflows/
│   └── main.wdl              # Entry point workflow
├── tasks/
│   └── cram_to_mhcbam.wdl   # Core extraction task
├── inputs/
│   └── inputs.json           # Input template (GCS paths)
├── tools/                    # Empty
├── README.md
└── LICENSE
```

---

## File Descriptions

### `workflows/main.wdl`

Top-level WDL 1.0 workflow. Imports the task and wires inputs/outputs.

**Inputs:**

| Name | Type | Description |
|---|---|---|
| `cram` | File | Whole-genome CRAM file |
| `crai` | File | CRAM index |
| `ref_fasta` | File | Reference genome FASTA (GRCh38) |
| `ref_fasta_fai` | File | FASTA index |
| `ref_dict` | File | Sequence dictionary |

**Outputs:**

| Name | Type | Description |
|---|---|---|
| `mhc_bam` | File | BAM containing MHC/HLA reads + mates |
| `mhc_bai` | File | BAM index |

---

### `tasks/cram_to_mhcbam.wdl`

Task `extract_mhc_with_mates`. Performs the actual read extraction using samtools.

**Processing steps:**

1. Symlink input files to standardized names (`sample.cram`, `ref.fa`, etc.)
2. Run `samtools idxstats` to identify alt contigs with mapped reads matching:
   - `chr6_*_alt`
   - `6_*_alt`
   - `HLA*`
3. Extract read names from the MHC region (`chr6:28510120-33480577`) and identified alt contigs
4. Extract read names of fully unmapped pairs (`-f 0xC`)
5. Merge all read names into a single list
6. Use `samtools view -N` to extract all reads (and mates) matching the list
7. Index the output BAM with `samtools index`

**Runtime:**

| Parameter | Value |
|---|---|
| Docker | `biocontainers/samtools:v1.19.2-1-deb_cv1` |
| CPU | 2 |
| Memory | 8 GB |
| Disk | 200 GB HDD (local-disk) |

---

### `inputs/inputs.json`

Input template for Terra submission. All paths use Google Cloud Storage (`gs://`).

```json
{
  "main.cram":         "gs://PATH_TO/sample.cram",
  "main.crai":         "gs://PATH_TO/sample.cram.crai",
  "main.ref_fasta":    "gs://PATH_TO/Homo_sapiens_assembly38.fasta",
  "main.ref_fasta_fai":"gs://PATH_TO/Homo_sapiens_assembly38.fasta.fai",
  "main.ref_dict":     "gs://PATH_TO/Homo_sapiens_assembly38.dict"
}
```

---

## Pipeline Workflow

```
Input: CRAM + CRAI + Reference (GRCh38)
        |
        v
  [main.wdl]
        |
        v
  [extract_mhc_with_mates]
   |
   |-- samtools idxstats  →  identify MHC alt contigs
   |-- samtools view      →  read names from chr6:28510120-33480577 + alt contigs
   |-- samtools view      →  read names of unmapped pairs
   |-- merge read name lists
   |-- samtools view -N   →  extract reads + mates
   |-- samtools index     →  index output BAM
        |
        v
Output: mhc_with_mates.bam + .bam.bai
```

---

## Compatibility Notes

### Terra-Specific

**Relative imports require ZIP upload**

Relative imports (`../tasks/...`) are only supported when the workflow is uploaded to Terra as a ZIP archive preserving the directory structure. Uploading `main.wdl` alone will fail.

### Warnings

**`samtools idxstats` without `-T` — `tasks/cram_to_mhcbam.wdl:31`**

Reference is not passed to `samtools idxstats`. Works if the CRAM embeds a reference URI, but fragile for CRAMs that require explicit reference specification.

**`-f 12` notation — `tasks/cram_to_mhcbam.wdl:50`**

Functionally correct (`0x4 | 0x8 = 12`), but `-f 0xC` is the conventional hex notation for samtools flag filtering.
