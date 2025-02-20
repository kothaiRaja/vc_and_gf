process EXTRACT_VCF {
    container params.container
	publishDir "${params.outdir}/annotations", mode: 'copy'

    input:
    path vcf_file
	path html

    output:
    path "extracted_vcf_data.csv"

    script:
    """
    python - <<EOF
import pandas as pd
import vcfpy

# File paths
vcf_file_path = "${vcf_file}"

# Initialize list to store parsed VCF data
vcf_rows = []

# Parse the VCF file
reader = vcfpy.Reader.from_path(vcf_file_path)

# Extract headers for samples
sample_names = reader.header.samples.names

# Extract data from VCF
for record in reader:
    chrom = record.CHROM
    pos = record.POS
    ref = record.REF
    alt = ",".join(str(a) for a in record.ALT)
    filt = ",".join(record.FILTER)
    dp = record.INFO.get('DP', 'NA')  # Extract Depth of Coverage
    ann_field = record.INFO.get('ANN', ['NA'])

    if ann_field != 'NA':
        ann_details = [a.split('|') for a in ann_field]
        for ann in ann_details:
            gene = ann[3] if len(ann) > 3 else "NA"
            impact = ann[2] if len(ann) > 2 else "NA"
            variant_type = ann[1] if len(ann) > 1 else "NA"

            # Iterate over samples to include sample-specific data
            for sample in sample_names:
                sample_data = record.call_for_sample[sample].data
                genotype = sample_data.get('GT', './.') if sample_data else './.'

                vcf_rows.append({
                    "CHROM": chrom,
                    "POS": pos,
                    "REF": ref,
                    "ALT": alt,
                    "FILTER": filt,
                    "DP": dp,
                    "Impact": impact,
                    "Gene": gene,
                    "Variant": variant_type,
                    "Sample_ID": sample,
                    "Genotype": genotype
                })

# Convert the VCF data to a DataFrame
vcf_df = pd.DataFrame(vcf_rows)

# Save VCF data to a file for reference
vcf_output_path = "extracted_vcf_data.csv"
vcf_df.to_csv(vcf_output_path, index=False)

print(f"Extracted VCF data saved to: {vcf_output_path}")
EOF
    """
}