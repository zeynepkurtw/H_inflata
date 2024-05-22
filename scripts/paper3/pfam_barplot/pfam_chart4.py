import plotly.graph_objects as go
import pandas as pd
import matplotlib.colors as mcolors

# Load InterProScan output data
#interproscan_file = "scripts/paper3/pfam_barplot/interpro_HIN_LRR_pfam.csv"
interproscan_file = "jupyter/Annas_script/multifamily/output/interpro_HIN_LRR_pfam.csv"
df = pd.read_csv(interproscan_file, sep='\t', header="infer")

# Determine the frequency of each Pfam domain
pfam_counts = df['Signature accession'].value_counts()

# Identify the most and least frequent Pfam domains
most_frequent_pfam = pfam_counts.idxmax()
least_frequent_pfam = pfam_counts.idxmin()

# Create a color gradient based on frequency
color_map = {}
cmap = mcolors.LinearSegmentedColormap.from_list("frequency", ["#D3D3D3", "#000000"])

for domain, count in pfam_counts.items():
    if domain == most_frequent_pfam:
        color_map[domain] = '#FF0000'  # Red for the most frequent domain
    elif domain == least_frequent_pfam:
        color_map[domain] = '#0000FF'  # Blue for the least frequent domain
    else:
        norm_count = (count - pfam_counts.min()) / (pfam_counts.max() - pfam_counts.min())
        color_map[domain] = mcolors.to_hex(cmap(norm_count))

# Create a helper function to determine unique architectures
def get_architecture_string(protein_df):
    return ';'.join(protein_df.sort_values(by='Start location')['Signature accession'])

# Add an architecture column to the dataframe
df['Architecture'] = df.groupby('protein_id').apply(get_architecture_string).reset_index(level=0, drop=True)

# Get unique architectures and their counts
unique_architectures = df['Architecture'].value_counts().reset_index()
unique_architectures.columns = ['Architecture', 'Count']

# Sort unique architectures based on their length (number of domains)
unique_architectures['Length'] = unique_architectures['Architecture'].apply(lambda x: len(x.split(';')))
unique_architectures = unique_architectures.sort_values(by='Length', ascending=False)

# Create a bar plot for each unique architecture
fig = go.Figure()

# Track domains that have been added to the legend
legend_domains = set()

# Define the fixed segment length for each domain
fixed_segment_length = 100

# Iterate through each unique architecture
for idx, row in unique_architectures.iterrows():
    architecture = row['Architecture']
    count = row['Count']
    prot_id = f'Architecture {idx + 1}'

    domains = architecture.split(';')
    start_position = 0

    # Add a line to represent the full length of the protein (behind the bars)
    fig.add_trace(go.Scatter(
        x=[0, len(domains) * (fixed_segment_length + 10)],
        y=[prot_id, prot_id],
        mode='lines',
        line=dict(color='lightgrey', width=2),
        showlegend=False
    ))

    for domain in domains:
        domain_color = color_map.get(domain, 'skyblue')  # Default color if not in color_map
        domain_description = df[df['Signature accession'] == domain]['Signature description'].iloc[0]

        show_legend = domain not in legend_domains
        if show_legend:
            legend_domains.add(domain)

        fig.add_trace(go.Bar(
            x=[fixed_segment_length],
            y=[prot_id],
            orientation='h',
            base=start_position,
            text=f"{domain_description} ({domain})",
            textposition='none',  # Hide text on the bars themselves
            hoverinfo='text',
            marker=dict(color=domain_color),
            name=domain,
            legendgroup=domain,
            showlegend=show_legend
        ))

        start_position += fixed_segment_length + 10  # Increment start position for consistent spacing

    # Add the count at the end of the line
    fig.add_trace(go.Scatter(
        x=[start_position + 200],
        y=[prot_id],
        text=[str(count)],
        mode='text',
        textfont=dict(size=6),  # Smaller text size
        showlegend=False
    ))

# Update layout for better visualization
fig.update_layout(
    title='Pfam Domains of Proteins with Unique Architectures',
    xaxis=dict(title='Protein Sequence Position'),
    yaxis=dict(title='Protein Architecture'),
    barmode='stack',
    showlegend=True,
    height=800,  # Adjust the height as needed
    width=1200  # Adjust the width as needed
)

fig.write_html('scripts/paper3/pfam_barplot/LRR_domain_arch.html')
# Display the plot
fig.show()

# Print statistics from data
num_proteins = df['protein_id'].nunique()
num_unique_architectures = unique_architectures.shape[0]
avg_protein_size = df.groupby('protein_id')['Sequence length'].mean().mean()

print(f"Number of proteins: {num_proteins}")
print(f"Number of unique domains: {len(pfam_counts)}")
print(f"Number of unique architectures: {num_unique_architectures}")
print(f"Average protein size: {avg_protein_size:.2f}")
