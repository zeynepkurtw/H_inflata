from ete3 import NCBITaxa
import pandas as pd
import argparse
import os

def check_lineage(taxon, outside_plancto):
    lineage = ncbi.get_lineage(taxon)
    names = ncbi.get_taxid_translator(lineage)
    
    if 'Planctomycetes' not in [names[taxid] for taxid in lineage]:
        outside_plancto=True
    return outside_plancto
    

def extract_exclusive_OG(path_blast_files, outfile_name):
    ncbi = NCBITaxa()
    outfile = open(outfile_name, 'w')
    blast_headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore', 'staxids', 'sscinames']
    
    # for each OG blast result file for each species of ineterest, see what OGs does not have hits outside planctomycetes
    for blast_file in os.listdir(path_blast_files):
        
        orthogroup = os.path.split(blast_file)[1]
        orthogroup = orthogroup.split('_')[0]
        blast_df = pd.read_csv(str(path_blast_folder) + str(blast_file), sep='\t', names=blast_headers)

        protein_results = blast_df.groupby('qseqid')

        for query, df in protein_results:
            outside_plancto = False
            taxa = df['staxids'].dropna().values.tolist()
            
            # check if the taxon/taxa in per hit is from planctomycetes
            for taxon in taxa:
                if ';' in str(taxon):
                    taxon = str(taxon).split(';')

                    for part in taxon:
                        outside_plancto = check_lineage(part, outside_plancto)
                else:
                    outside_plancto = check_lineage(taxon, outside_plancto)
                    
        # if all the hits to this OG proteins do not come from outside planctomycetes in the NR database save it to our exclsuive OG file
        if outside_plancto == False:
            outfile.write(orthogroup + '\n')
    
    outfile.close()
    return outfile

extract_exclusive_OG(snakemake.input, snakemake.output)


