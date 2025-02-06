import os
import argparse
import gzip

def parse_vcf(snpeff_file, output_file=None, min_dp=0, min_af=0, build="GRCh38", filter_pass=True):
    if output_file is None:
        output_file = f"{os.path.splitext(snpeff_file)[0]}.maf"
    
    # Handle gzipped VCF files
    if snpeff_file.endswith('.gz'):
        infile = gzip.open(snpeff_file, 'rt', encoding='utf-8')
    else:
        infile = open(snpeff_file, 'r', encoding='utf-8')
    
    with infile, open(output_file, 'w', encoding='utf-8') as outfile:
        outfile.write("Hugo_Symbol\tEntrez_Gene_Id\tCenter\tNCBI_Build\tChromosome\tStart_Position\tEnd_Position\tStrand\tVariant_Classification\tVariant_Type\tReference_Allele\tTumor_Seq_Allele1\tTumor_Seq_Allele2\tTumor_Sample_Barcode\tProtein_Change\ti_TumorVAF_WU\ti_transcript_name\n")
        
        for line in infile:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if fields[4] == '.' or (filter_pass and fields[6] != "PASS"):
                continue
            
            chromosome, start_pos, ref_allele, alt_allele = fields[0], int(fields[1]), fields[3], fields[4]
            chromosome = chromosome.replace("chr", "")
            end_pos = start_pos
            
            if len(ref_allele) == 1 and len(alt_allele) == 1:
                var_type = "SNP"
            elif len(ref_allele) > len(alt_allele):
                var_type = "DEL"
            elif len(ref_allele) < len(alt_allele):
                var_type = "INS"
            else:
                var_type = "Complex"
            
            tumor_seq_allele1 = ref_allele
            tumor_sample_barcode = os.path.splitext(snpeff_file)[0]
            
            variant_classification, hugo_symbol, transcript_name, protein_change = "NA", "NA", "NA", "NA"
            annotation_info = next((f for f in fields[7].split(';') if f.startswith("ANN=")), None)
            if annotation_info:
                ann_fields = annotation_info.split('|')
                variant_classification = ann_fields[1] if len(ann_fields) > 1 else "NA"
                hugo_symbol = ann_fields[3] if len(ann_fields) > 3 else "NA"
                transcript_name = ann_fields[6] if len(ann_fields) > 6 else "NA"
                protein_change = ann_fields[10] if len(ann_fields) > 10 else "NA"
            
            format_fields = fields[8].split(':')
            sample_values = fields[-1].split(':')
            
            dp_index = format_fields.index('DP') if 'DP' in format_fields else -1
            ad_index = format_fields.index('AD') if 'AD' in format_fields else -1
            
            dp = int(sample_values[dp_index]) if dp_index >= 0 and sample_values[dp_index].isdigit() else "NA"
            ad = int(sample_values[ad_index].split(',')[-1]) if ad_index >= 0 and ',' in sample_values[ad_index] else "NA"
            
            tumor_vaf = f"{ad / dp:.5f}" if dp != "NA" and ad != "NA" and dp > 0 else "NA"
            
            if dp != "NA" and dp < min_dp:
                continue
            if tumor_vaf != "NA" and float(tumor_vaf) < min_af:
                continue
            
            outfile.write(f"{hugo_symbol}\tNA\tNA\t{build}\t{chromosome}\t{start_pos}\t{end_pos}\tNA\t{variant_classification}\t{var_type}\t{ref_allele}\t{tumor_seq_allele1}\t{alt_allele}\t{tumor_sample_barcode}\t{protein_change}\t{tumor_vaf}\t{transcript_name}\n")
    
    print(f"Output written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert SnpEff VCF to MAF format")
    parser.add_argument("snpeff_file", help="Input SnpEff-annotated VCF file")
    parser.add_argument("-o", "--output", type=str, help="Optional output MAF file")
    parser.add_argument("-minDP", type=int, default=0, help="Minimum depth (default: 0)")
    parser.add_argument("-minAF", type=float, default=0, help="Minimum allele frequency (default: 0)")
    parser.add_argument("--build", type=str, default="GRCh38", help="NCBI Build (default: GRCh38)")
    parser.add_argument("--filterPASS", action='store_true', default=True, help="Filter for PASS (default: True)")
    
    args = parser.parse_args()
    
    parse_vcf(args.snpeff_file, args.output, args.minDP, args.minAF, args.build, args.filterPASS)
