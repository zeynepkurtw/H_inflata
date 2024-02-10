
import argparse
import pandas as pd
parser = argparse.ArgumentParser()
from collections import Counter 


def open_file(my_file, header, seperator, col_names):
    if col_names is not None:
        df = pd.read_csv(my_file, sep=seperator, header=header, names=col_names) 
    else:
        df = pd.read_csv(my_file, sep=seperator, header=header)
    df = df.fillna('Missing')

    return df

def open_interproscan_file(interproscan_file, analysis):#, selected_domains):
    """
    extract the Pfam domains of a interproscan file
    """
    interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession', 
                           'Signature description','Start location', 'Stop location', 'Score', 'Status', 
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description', 
                           'optional: GO terms', 'optional: Pathway annotation']

    domain_df = pd.read_csv(interproscan_file, sep='\t', header=None, names = interproscan_header)
    domain_df = domain_df.fillna('Missing')
    #pids = domain_df[domain_df['Signature accession'].isin(selected_domains)]['protein_id'].values.tolist()
    #domains_df = domain_df[domain_df['protein_id'].isin(pids)]
    domain_df = domain_df[domain_df['Analysis']==analysis]
    domain_df = domain_df.sort_values(by=['protein_id','Start location'])
    
    return domain_df

def open_file_tsv(in_file):
    df = pd.read_csv(in_file, sep='\t')
    return df

def open_file(in_file):
    df = pd.read_csv(in_file)
    return df

def reshape_interproscan_file(interpro_df):
    interpro_df['Start location'] = interpro_df['Start location'].astype(str) 
    interpro_df['Stop location'] = interpro_df['Stop location'].astype(str)

    interpro_df_pid = interpro_df.groupby('protein_id').agg({'Signature accession': lambda x: ';'.join(x)})
    interpro_df_pid = interpro_df_pid.reset_index()


    interpro_df_anno = interpro_df.groupby('protein_id').agg({'Signature description': lambda x: ';'.join(x)})
    interpro_df_anno = interpro_df_anno.reset_index()

    interpro_df = interpro_df_pid.join(interpro_df_anno.set_index('protein_id'), on='protein_id')
    #interpro_df = interpro_df[interpro_df['Signature accession'].str.contains('PF00069')]
    ## one protein per line and the correct order, sort on protein and start values before joining the rows
    return interpro_df #_pid

def get_domain_arc(interpro_df):
    domains_all = interpro_df['Signature accession']
    cnt_our = Counter(domains_all.values.tolist())
    domains = domains_all.drop_duplicates().values.tolist()
    return domains, cnt_our

def create_schematic_table(domains_list, cnt_our, interpro_df, interpro_original_df): #, cnt_homolog):
    schematic_df = []
    prot_id = len(domains_list)

    for domain_arc in domains_list:

        amount_our = cnt_our[domain_arc]
        #amount_homolog = cnt_homolog[domain_arc]

        protein_ids = interpro_df[interpro_df['Signature accession']==domain_arc]['protein_id'].sort_values().str.cat(sep=';').rstrip() #values.tolist().replace("'")

        if ';' in domain_arc:
            prot_domains = domain_arc.split(';')
            prot_len = len(prot_domains)*150
            start = 25
            stop = 125

            for domain in prot_domains:
                domain_anno = interpro_original_df[interpro_original_df['Signature accession']==domain]['Signature description'].drop_duplicates().values.tolist()[0]
                schematic_df.append(['Architecture ' + str(prot_id), prot_len, domain, start, stop, amount_our, protein_ids, domain_anno]) #, amount_homolog])
                start+=125
                stop+=125
        else:
            prot_len = 150
            start = 25
            stop = 125
            domain_anno = interpro_original_df[interpro_original_df['Signature accession']==domain_arc]['Signature description'].drop_duplicates().values.tolist()[0]
            schematic_df.append(['Architecture ' + str(prot_id), prot_len, domain_arc, start, stop, amount_our, protein_ids, domain_anno]) #, amount_homolog])
            
        prot_id-=1 # the protein id counter must descrese here to get the correct order in the figure, smallest number at the top of the figure
    schematic_df = pd.DataFrame(schematic_df)

    return schematic_df

def number_of_dom(domain):
    val = 1
    if ';' in domain:
        val = len(domain.split(';'))
    return val

def dom_loc(domain):
    val = 1
    if ';' in domain:
        domain = domain.split(';')
        for e in domain:
            if e in selected_domains:
                val +=1
                return val
            val+=1
    return val



def merge_similar_domains(domains_to_merge, interpro_df):
    
    domains_to_merge_df = pd.read_csv(domains_to_merge, sep='\t')
    #domain_id_group_df = domains_to_merge_df[['Domain id','Merged domain id']] # removed
    domain_id_group_df = domains_to_merge_df[['Domain arc id','Merged domain arc id', 'Domain annotation merged']]

    domain_id_group_df.columns = ['Signature accession','Merged domain id', 'Merged annotation']

    interpro_df_added = interpro_df.set_index('Signature accession').join(domain_id_group_df.set_index('Signature accession')).reset_index()
    ##interpro_df_added = interpro_df.merge(domain_id_group_df, how='left', on='Signature accession')

    interpro_df = interpro_df_added.drop(columns=['Signature accession'])
    interpro_df = interpro_df.drop(columns=['Signature description']) # added

    interpro_df = interpro_df.rename(columns={"Merged domain id": "Signature accession"})
    interpro_df = interpro_df.rename(columns={"Merged annotation": "Signature description"}) # added

    interpro_df['Signature accession'] = interpro_df['Signature accession'].astype(str)

    #interpro_df.reset_index()
    return interpro_df

def main(interproscan_file, analysis, out_file): #, merged_domains):#, domains_to_merge):

    interpro_original_df = open_interproscan_file(interproscan_file, analysis)#, selected_domains)

    interpro_df = reshape_interproscan_file(interpro_original_df)

    domain_list, cnt_our = get_domain_arc(interpro_df)

    #domain_list.sort(key=dom_loc, reverse=True)
    
    domain_list.sort(key=number_of_dom, reverse=True)
    schematic_df = create_schematic_table(domain_list, cnt_our, interpro_df, interpro_original_df)

    f = schematic_df.to_csv(out_file, header=None, index=False, sep='\t')
    return f


# get input and run the script
parser = argparse.ArgumentParser()
interpro_path = parser.add_argument("--interpro_path")
out_file_path = parser.add_argument("--out_file_path")
analysis = parser.add_argument("--analysis")
#og_orthomcl = parser.add_argument("--orthomcl")
#og_orthofinder = parser.add_argument("--orthofinder")

#pfam_domain_selection = parser.add_argument("--pfam_domain_selection")
#merged_domains = parser.add_argument("--merged_domains") #--merge_domains 'Domain_anno_color_v3.txt'


args = parser.parse_args()

#global selected_domains
#selected_domains = ['PF00069'] #['PF00069'], IPR011009 ['PF04542', 'PF08281', 'PF04539', 'PF04545']

main(args.interpro_path, args.analysis, args.out_file_path)#, args.merged_domains)



