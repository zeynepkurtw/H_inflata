import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import numpy as np
import matplotlib
from matplotlib import cm
import matplotlib.colors as c
import pandas as pd
import numpy as np
import argparse
parser = argparse.ArgumentParser()

############### Open files ############################################

def open_file(my_file, general_domains):
    """
    open the interproscan file as a dataframe
    """
    interproscan_header = ['protein_id','Seq MD5 digest','Sequence length','Analysis', 'Signature accession',
                           'Signature description','Start location', 'Stop location', 'Score', 'Status',
                           'Date of run', 'Interpro annotations - accession', 'Interpro annotation - description',
                           'optional: GO terms', 'optional: Pathway annotation']

    my_header = ['Architecture_id','Sequence length','Signature accession','Start location', 'Stop location', 'Abundance in interproscan file', 'Protein_id', 'Signature description']

    if general_domains == 'No':
        df = pd.read_csv(my_file, header = None, sep='\t', names = interproscan_header)
    else:
        df= pd.read_csv(my_file, header = None, sep='\t', names = my_header)
    
    df = df.fillna('Missing')
    return df

def open_TM(TM_file):
    """
    open the Transmembrane prediction file
    """
    TM_data = pd.read_csv(TM_file, header = None, sep='\t')
    return TM_data

############### Extract info from files ############################################

def extract_analysis(domain_file, analysis, general_domains):
    """
    extract the domains of the interproscan file, default the analysis is Pfam
    """
    
    interpro_df = open_file(domain_file, general_domains)
    
    if general_domains == 'No':
        domain_df = interpro_df[interpro_df['Analysis']==analysis]
    else:
        domain_df = interpro_df

    return domain_df


def extract_protein_ids(protein_data):
    """
    extract the protein ids
    """
    protein_ids = protein_data['Architecture_id']
    protein_ids = protein_ids.drop_duplicates()
    return protein_ids

def protein_length(protein_id, protein_data):
    """
    extract the protein length
    """
    p_length = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()].drop_duplicates()['Sequence length']
    p_length = int(p_length.values.tolist()[0])
    return p_length


def protein_amount(protein_id, protein_data):
    """
    extract the protein length
    """
    amount_s = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()].drop_duplicates()['Abundance in interproscan file'].values.tolist()[0]
    #amount_h = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()].drop_duplicates() ['Abundance in homolog'].values.tolist()[0]
    #amount = str(amount_s) + ' / ' + str(amount_h) 
    return amount_s

def extract_domains(protein_id, protein_data):
    """
    extract the domains with domain id, start and stop for each domain prediction per protein
    """
    p_info = protein_data[protein_data['Architecture_id'] == str(protein_id).rstrip()]
    domain_info = p_info[['Signature accession','Start location', 'Stop location']]
    return domain_info

def extract_TM(protein_id, TM_data):
    """
    extract the transmembrane info / signaling too?
    """
    TM = TM_data[TM_data['protein_id'] == str(protein_id).rstrip()]
    return TM

def get_groups(interpro_df, protein_id, groups):
    group_list = []
    group_temp = ''
    
    real_protein_ids = interpro_df[interpro_df['Architecture_id']==protein_id]['Protein_id'].values.tolist()[0]
    
    for protein in real_protein_ids.split(';'): 
        group_temp = groups[groups['protein_id'].str.contains(protein)]['specie'].drop_duplicates().values.tolist()[0]
        
        if group_temp not in group_list:
            group_list.append(group_temp)
        group_temp = ''
    return group_list

############### Figure text ############################################


def add_text_2_fig(general_domain, grid, section, protein_name, p_length, p_amount, space):
    """
    Add protein length and protein name as text to the images
    """
    if general_domains == 'No':
        plt.text(grid[section][0] + p_length + 100, grid[section][1],  str(p_length), fontsize = 10, horizontalalignment = 'center')  
    else:
        plt.text(grid[section][0] + p_length + 100, grid[section][1],  str(p_amount), fontsize = 10, horizontalalignment = 'center') 
    plt.text(grid[section][0] - space, grid[section][1], protein_name + '  ', fontsize = 10, horizontalalignment = 'right')
    return plt 

def label(xy,text):
    """
    Create the text label layout
    """
    plt.text(xy[0], xy[1] -15, text, ha='left', va='center',horizontalalignment= 'right',rotation =0, family='sans-serif', size = 8)
    return plt

def label_bold(xy,text):
    """
    Create the bold text label layout
    """
    plt.text(xy[0], xy[1] -15, text, ha='left', va='center', family='sans-serif', size = 8, weight='bold')
    return plt

def create_lable(grid, section, lable, patches, margin):
    """
    create lable text box of color and domain ids
    """
    plt.text(grid[section][0] + 2500, grid[section][1] - margin, str(lable), ha='left', va='center', family='sans-serif', size = 16) #, weight='bold')
    
    return plt

def create_lable_rect(grid, section, lable, color, patches, colors, margin):
    rect = mpatches.Rectangle([grid[section][0]+2400, grid[section][1]- margin - 25], 50, 50) #50 is the hight of the domain
    patches.append(rect)
    
    #circle = mpatches.Ellipse((int(grid[section][0] + 2400), int(grid[section][1] - margin)), 25, 100)
    #patches.append(circle)
    colors.append(color)
    return patches, colors

def create_lable_circle(grid, section, lable, color, patches, colors, margin):
    circle = mpatches.Ellipse((int(grid[section][0] + 2425), int(grid[section][1] - margin)), 25, 100)
    patches.append(circle)
    colors.append(color)
    return patches, colors
############# Add shapes; lines and rectangles representing protein length, domains and transmembrane regions

def create_line(grid, section, p_length, patches, colors):
    """
    create the line representing the protein length
    """
    rect = mpatches.Rectangle([grid[section][0],grid[section][1]+25],p_length,5) # create the black line that symbolises sequence length (grid[section][1]+25) 
    patches.append(rect)
    colors.append([0.0,0.0,0.0])

    return patches, colors

def add_domains_2_fig(grid, section, domain_info, patches, colors):
    """
    Add domains to the figure as orange rectangles
    """
    
    top_10_domain_id_list_keys = top_10
    global colors_for_top_10_values
    colors_for_top_10_values = ['#FFBF65', '#C05780', '#8DD78F', '#FF96C5', '#E77577', '#00B0BA', '#FFD872', '#4DD091', '#6C88C4', '#FF828B']
    
    color_dict = dict(zip(top_10_domain_id_list_keys, colors_for_top_10_values))
        
        
    domain = str(domain_info[0]) # name
    index = [int(domain_info[1]), int(domain_info[2])] # start and stop
    rect = mpatches.Rectangle([grid[section][0]+index[0], grid[section][1]], (index[1] - index[0]), 50) #50 is the hight of the domain
    label([grid[section][0]+index[0], grid[section][1]-44], domain)
    patches.append(rect)
    
    #if domain in ['PF00069', 'PF04542', 'PF08281', 'PF04539', 'PF04545']:
    #    colors.append('#8762aa')#[245/255.0,145/255.0,64/255.0]) # make the domains scilifelab-orange
    #else:
    #    colors.append('#6699cc')#[17/255.0,134/255.0,195/255.0])
    
    if domain not in color_dict.keys():
        colors.append('#00A5E3')
    else:
        colors.append(color_dict[domain])
    return patches, colors

def add_group_circle_2_fig(grid, section, group, patches, colors, p_length, cnt_group):   
    """
    add what group has this domain architecture
    """
    global color_dict
    color_dict = {'HI': 'mediumaquamarine',
                 'SS':'pink',
                 'GM':'crimson',
                 'GL': 'darkorange',
                 'TP': 'blue',
                  #'GI': 'green',
                  #'GC' : 'green',
                  #'KA': 'black',
                 }
   #'species?': 'gray'

    
    buffer = 50*cnt_group
    circle = mpatches.Ellipse((int(grid[section][0]+ p_length  + 100  + int(buffer)), int(grid[section][1])),25,100)
    patches.append(circle)
    colors.append(color_dict[group])
    
    return patches, colors

def add_TM_2_fig(TM, section, grid, patches, colors):
    """
    Add transmembrane regions to the figure as blue rectangles
    """
    rect = mpatches.Rectangle([grid[section][0] + int(TM[1]), grid[section][1]], (int(TM[2])-int(TM[1])), 50)
    patches.append(rect)
    colors.append([17/255.0,134/255.0,195/255.0])

    return patches, colors

def add_signal_2_fig(SIG, section, grid, patches, colors):
    """
    Add transmembrane regions to the figure as blue rectangles
    """
    rect = mpatches.Rectangle([grid[section][0] + int(SIG[1]), grid[section][1]], (int(SIG[2])-int(SIG[1])), 50)
    patches.append(rect)
    colors.append([17/255.0,134/255.0,195/255.0]) # change color to green

    return patches, colors

############## Save figure with final layout settings ##################

def save_fig(protein_ids, out_name, patches, colors, plt, ax, x_lim, y_lim):
    """
    Save the figure with layout options
    """
    collection = PatchCollection(patches)
    collection.set_color(colors)
    ax.add_collection(collection)
    ax.figure.set_size_inches(60,0.5*int(len(protein_ids)))# 25 for pvc
    ax.set_xlim([-200,x_lim]) #1500
    ax.set_ylim([-2,y_lim]) #30000
    #plt.axis('off')
    ax.axis('off')
    #ax.get_xaxis().set_visible(False)
    #ax.get_yaxis().set_visible(False)
    
    #plt.figure.frameon(False)
    plt.savefig(out_name, orientation = 'portrait', bbox_inches='tight')

    return plt



def create_protein_domain(sec, general_domain, patches, colors, section, protein_id, interpro_df, grid, groups):#, TM_data):
    """
    create a protein domain figure for one protein
    """
    domain_info = extract_domains(protein_id, interpro_df).values.tolist()
    #print('here')
    #print(protein_id)
    
    p_length = protein_length(protein_id, interpro_df)
    p_amount = 0

    if general_domain == 'Yes':
        p_amount = protein_amount(protein_id, interpro_df)
    
    ## add protein id and protein length as text in the figure
    add_text_2_fig(general_domain, grid, section, str(protein_id), p_length, p_amount, 90)
    
    ## add a black line representing the protein length
    patches, colors = create_line(grid, section, p_length, patches, colors)

    ## add domain region in the protein

    for domain in domain_info:
        patches, colors = add_domains_2_fig(grid, section, domain, patches, colors)
        
    ## add the group as a colored circle:
    cnt_group = 1
    for group in groups:
        patches, colors = add_group_circle_2_fig(grid, section, group, patches, colors, p_length, cnt_group)
        cnt_group+=1
    
    ## add TM region in blue in the protein
    #if TM_data != []:
        #TM_info = []
    #    TM_info = extract_TM(protein_id, TM_data).values.tolist()
    #    for TM in TM_info:
    #        patches, colors = add_TM_2_fig(TM, section, grid, patches, colors)

    ## add SIG region in green
    
    #section += sec

    return patches, colors, section, grid

def get_anno(top_10, interpro_df):
    anno_list = []
    for domain_id in top_10:
        anno_list.append(interpro_df[interpro_df['Signature accession']==domain_id]['Signature description'].values.tolist()[0])
    return anno_list

def generate_domain_fig(sec, interpro_df, protein_ids, general_domain, groups):#, TM_data):
    """
    Generate the domain figures for all proteins
    """
    ## figure variables
    fig, ax= plt.subplots()
    #ax = plt.axes()
    patches=[]
    colors=[]
    grid = np.mgrid[0:1,0:int(201 * len(protein_ids))].reshape(2,-1).T 
    #grid = np.mgrid[0:1:1j,0:4500:540j].reshape(2,-1).T # create grid to build the images on
    section = 0
    
    top_10_domain_colors = interpro_df
    top_10_domain_colors = top_10_domain_colors[top_10_domain_colors['Signature accession']!='Missing']
    top_10_domain_colors = top_10_domain_colors[['Signature accession']].value_counts().reset_index()
    top_10_domain_colors.columns = ['Signature accession', 'Count']
    global top_10
    top_10 = top_10_domain_colors.iloc[0:10,0].values.tolist()
    top_10_anno = get_anno(top_10, interpro_df)
                  

    # create domain figure for each protein in the interproscan file
    for protein_id in protein_ids:
        section+=sec
        
        #rint(protein_id)
        
        # if schematic table as input:
        group_list = get_groups(interpro_df, protein_id, groups)
        #print(group_list)
        

        patches, colors, section, grid = create_protein_domain(sec, general_domain, patches, colors, section, protein_id, interpro_df, grid, group_list)#, TM_data)
    
    margin = 0
    counter = 0
    

    for lable in top_10_anno:
        create_lable(grid, section, str(lable), patches,margin)
        patches, colors = create_lable_rect(grid, section, lable, colors_for_top_10_values[counter], patches, colors, margin)
        margin += 150
        counter+=1
        
    # create a lable for the domains not in the top 10, that all has the same color:
    create_lable(grid, section, 'Other domains', patches,margin)
    patches, colors = create_lable_rect(grid, section, 'Other domains', '#00A5E3', patches, colors, margin)
    margin += 300

    
    for group in groups['specie'].drop_duplicates().values.tolist():
        create_lable(grid, section, str(group), patches, margin)
        patches, colors = create_lable_circle(grid, section, str(group), color_dict[group], patches, colors, margin)
        margin += 150
                
    return patches, colors, fig, ax

def main(Interproscan, analysis, out_name, general_domains, groups):
    """
    Generates a domain architecture figure of proteins in an interproscan file. 
    Use the input flag analysis to change between domain prediction softwares, default is Pfam.
    """

    interpro_df = extract_analysis(Interproscan, analysis, general_domains)
    protein_ids = extract_protein_ids(interpro_df)
    
    x_lim = int(max(interpro_df['Sequence length'])) + 1500 #110
    y_lim = len(protein_ids)* 210#201#240
    sec = 200
    patches, colors, plt, ax = generate_domain_fig(sec, interpro_df, protein_ids, general_domains, groups)

    f = save_fig(protein_ids, out_name, patches, colors, plt, ax, x_lim, y_lim)

    return f

############## Run script ##################
parser.add_argument("--out_file")
parser.add_argument("--interproscan")
parser.add_argument("--analysis")
parser.add_argument("--general_domains")
parser.add_argument("--groups_protein_ids")

args = parser.parse_args()


Interproscan = args.interproscan
out_name = args.out_file
analysis = str(args.analysis)
general_domains = str(args.general_domains)
groups_protein_ids = pd.read_csv(args.groups_protein_ids, sep='\t')[['specie','protein_id']]

main(Interproscan, analysis, out_name, general_domains, groups_protein_ids)


