# pySnpEffToMaf
Simple tool to convert SnpEff annotated VCFs to MAF files to be used with maftools. Inspired by [snpeffToMaf](https://github.com/tsy19900929/snpeffToMaf).

## Usage

```
# Simple
python snpeffToMaf.py input.vcf/input.vcf.gz

# Additional parameters:
-h, --help            show this help message and exit
-o OUTPUT, --output OUTPUT  Optional output MAF file
-minDP MINDP          Minimum depth (default: 0)
-minAF MINAF          Minimum allele frequency (default: 0)
--build BUILD         NCBI Build (default: GRCh38)
--filterPASS          Filter for PASS (default: True)
```

## Loading MAF output 

Loading the maf file might require reclassifying some of the annotations. Here's how you can do that:

```
process_maf <- function(file) {
  # Read the MAF file
  maf_df <- read.delim(file, comment.char = "#")

  # Standardize Variant_Classification
  maf_df$Variant_Classification <- recode(maf_df$Variant_Classification,
    "frameshift_variant" = "Frame_Shift_Ins",
    "frameshift_variant&splice_region_variant" = "Frame_Shift_Ins",
    "frameshift_variant&splice_acceptor_variant&splice_donor_variant&splice_region_variant&intron_variant" = "Splice_Site",
    "missense_variant" = "Missense_Mutation",
    "missense_variant&splice_region_variant" = "Splice_Site",
    "stop_gained" = "Nonsense_Mutation",
    "stop_gained&conservative_inframe_insertion" = "Nonsense_Mutation",
    "splice_region_variant" = "Splice_Site",
    "splice_region_variant&intron_variant" = "Splice_Site",
    "splice_region_variant&non_coding_transcript_exon_variant" = "Splice_Site",
    "splice_region_variant&synonymous_variant" = "Silent",
    "synonymous_variant" = "Silent"
  )

  return(maf_df)
}

"output.maf" |>
  process_maf() |>
  read.maf() |>
  plotmafSummary(, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE, top = 20)
```
