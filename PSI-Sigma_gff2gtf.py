import argparse

def gff_to_gtf(gff_filename, gtf_filename):
    with open(gff_filename, 'r') as gff_file, open(gtf_filename, 'w') as gtf_file:
        for line in gff_file:
            if line.startswith('#') or line.strip() == '':
                continue

            fields = line.strip().split('\t')
            attributes = fields[8]

            # Extract ID, Name, and Parent if available
            attr_dict = {attr.split('=')[0]: attr.split('=')[1] for attr in attributes.split(';') if '=' in attr}
            gene_id = attr_dict.get('ID', '')
            gene_name = attr_dict.get('Name', '')
            transcript_id = attr_dict.get('Parent', '')

            # Change 'mRNA' feature type to 'transcript' for GTF
            feature_type = fields[2]
            if feature_type == 'mRNA':
                feature_type = 'transcript'

            # Construct GTF attribute field
            gtf_attributes = []
            if feature_type == 'gene':
                gtf_attributes.append(f'gene_id "{gene_id}"')
                gtf_attributes.append(f'gene_name "{gene_name}"')
            elif feature_type == 'transcript':
                gtf_attributes.append(f'gene_id "{transcript_id}"')
                gtf_attributes.append(f'transcript_id "{gene_id}"')
            else: # exon or CDS
                gtf_attributes.append(f'gene_id "{transcript_id.split(".")[0]}"')
                gtf_attributes.append(f'transcript_id "{transcript_id}"')
                exon_number = gene_id.split(".")[-1].lstrip('exon')
                gtf_attributes.append(f'exon_number "{exon_number}"')
                if feature_type == 'exon':
                    gtf_attributes.append('transcript_biotype "unknown"')

            gtf_attributes_str = '; '.join(gtf_attributes) + ';'
            gtf_line = '\t'.join(fields[:2] + [feature_type] + fields[3:8]) + '\t' + gtf_attributes_str + '\n'
            gtf_file.write(gtf_line)

def main():
    parser = argparse.ArgumentParser(description="Convert GFF file to GTF format.")
    parser.add_argument("--gff", required=True, help="Input GFF file")
    parser.add_argument("--gtf", required=True, help="Output GTF file")
    args = parser.parse_args()

    gff_to_gtf(args.gff, args.gtf)

if __name__ == "__main__":
    main()
