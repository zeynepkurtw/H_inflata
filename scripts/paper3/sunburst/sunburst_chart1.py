import plotly.graph_objects as go

# Define the data based on the provided information
labels = ["Proteome",
          "Functional", "Hypothetical",
          "Conserved Hypothetical",
          "Top7",
          "LRR", "HB", "CP", "PK", "CRP", "WD40", "Ankyrin",
          "D1_LRR", "D2_LRR", "D3_LRR",
          "D1_HB", "D2_HB", "D3_HB"]
parents = ["",
           "Proteome", "Proteome",
           "Hypothetical",
           "Functional",
           "Top7", "Top7", "Top7", "Top7", "Top7", "Top7", "Top7",
           "LRR", "LRR", "LRR",
           "HB", "HB", "HB"]
values = [79341,
          29874, 49467,
          14945,
          10477,
          3702, 3568, 885, 789, 682, 432, 419,
          1473, 432, 281,
          2020, 43, 25]

# Create the sunburst chart
fig = go.Figure(go.Sunburst(
    labels=labels,
    parents=parents,
    values=values,
    branchvalues="total",
    hovertemplate='%{label}<br>%{value} (%{percentParent:.2%})<extra></extra>',
    marker=dict(colors=[
        "white",
        "#FF9999", "#D3D3D3",
        "#A9A9A9",
        "#FF9999",
        "#FF3333", "#FF4D4D", "#FF6666", "#FF8080", "#FF9999", "#FFB3B3", "#FFCCCC",
        "#FF3333", "#FF3333", "#FF3333",
        "#FF4D4D", "#FF4D4D", "#FF4D4D"],
        line=dict(width=2)
    )
))

# Update layout for better visualization
fig.update_layout(
    title="Proteome Composition of H. inflata",
    margin=dict(t=40, l=0, r=0, b=0)
)

fig.write_html('scripts/paper3/plot.html')
# Display the plot
fig.show()
