import pandas as pd
import matplotlib.pyplot as plt
from Bio import Phylo
from matplotlib.patches import Patch

# File paths
data_file = '../data/all_events_1_10_100.txt'
tree_file = '../data/treerooted.nwk'

# Read deletion data
df = pd.read_csv(data_file, sep='\t')
df.columns = df.columns.str.strip()

# Print column names to verify
print("Columns in df:", df.columns.tolist())

# Clean values in 'State' column
df['State'] = df['State'].str.strip().str.lower()

# Parse phylogenetic tree
tree = Phylo.read(tree_file, 'newick')

# Assign Y positions to all nodes (internal and terminal)
def assign_y_positions(clade, positions_dict, current_y=0):
    if clade.is_terminal():
        positions_dict[clade] = current_y
        return current_y + 1
    else:
        ys = []
        for subclade in clade.clades:
            current_y = assign_y_positions(subclade, positions_dict, current_y)
            ys.append(positions_dict[subclade])
        positions_dict[clade] = sum(ys) / len(ys)
        return current_y

node_positions = {}
assign_y_positions(tree.root, node_positions)

# Create a mapping of node names to Y positions
# Also convert "nod_" to "node_" in the mapping
name_to_y = {}
for clade in node_positions:
    if clade.name:
        name = clade.name
        # Change "nod_" to "node_" in the tree as well
        if isinstance(name, str) and name.startswith("nod_"):
            name = name.replace("nod_", "node_", 1)
        name_to_y[name] = node_positions[clade]

# Draw tree without node names
def draw_clade(clade, x_start, ax, positions_dict):
    x_here = x_start + 1
    y_here = positions_dict[clade]
    if clade.clades:
        for child in clade.clades:
            y_child = positions_dict[child]
            # Draw vertical line to child
            ax.plot([x_here, x_here], [y_here, y_child], color='k')
            # Draw horizontal line to child
            ax.plot([x_here, x_here + 1], [y_child, y_child], color='k')
            draw_clade(child, x_here, ax, positions_dict)
    # Optional: add node labels if desired

# Calculate total number of nodes to adjust figure size
total_nodes = max(node_positions.values()) + 1

# Adjust figure size
fig_height = total_nodes * 0.25
fig = plt.figure(figsize=(30, 70))

# Create axes for tree and deletions, adjusting positions to increase space between graphs
ax_tree = fig.add_axes([0.05, 0.25, 0.15, 0.65])

# Adjust ax_del and ax_del_full_genes to increase space between graphs
# We slightly increase the 'left' value and adjust the 'width' accordingly
ax_del = fig.add_axes([0.26, 0.25, 0.33, 0.65], sharey=ax_tree)
ax_del_full_genes = fig.add_axes([0.64, 0.25, 0.33, 0.65], sharey=ax_tree)

# Draw tree
draw_clade(tree.root, 0, ax_tree, node_positions)
ax_tree.axis('off')
ax_tree.set_ylim(-1, total_nodes)

# Gene data with specified colors
genes = pd.DataFrame({
    'Gene': ['5´UTR', 'ORF1ab', 'NS2', 'HE', 'ORF2', 'S', '4.9/4.8kDaNSP', 'NSP3/ORF3',
             'NS4/15kDaNSP', 'ORF4/5', 'E', 'M', 'NS6', 'ORFx', 'ORF5', 'ORF6', 'ORF7',
             'ORF8', 'ORF9', 'ORF10', 'N2', 'N', 'IP', 'p10', 'ORF7/3x-lp', 'ORF8',
             'ORF9', 'ORF10', 'ORF11', 'ORF', 'ORFx', 'ORFy', '3´UTR'],
    'Start': [0, 1017, 31842, 33115, 34466, 35971, 45085, 45469, 47014, 47932,
              49820, 50246, 52427, 52938, 53302, 54776, 55488, 56320, 56892,
              57353, 58097, 58790, 61550, 62905, 63198, 65975, 67783, 68236,
              68532, 69075, 69288, 69851, 70385],
    'End': [1016, 30926, 32687, 34450, 35821, 43557, 45389, 46925, 47758,
            49740, 50176, 51778, 52753, 53222, 54772, 55480, 56167, 56817,
            57350, 57985, 58789, 60911, 62218, 63180, 65642, 66637, 68080,
            68531, 69074, 69287, 69850, 70384, 71920],
    'Color': ['#FF0000', '#008000', '#FFA500', '#0000FF', '#FFD700', '#800080',
              '#00FA9A', '#FF69B4', '#808000', '#6495ED', '#D2691E', '#5F9EA0',
              '#FF00FF', '#A52A2A', '#7FFF00', '#8A2BE2', '#F08080', '#4B0082',
              '#FF7F50', '#008080', '#CD5C5C', '#ADFF2F', '#A6CEE3', '#33A02C',
              '#DEB887', '#4B0082', '#FF7F50', '#008080', '#000080', '#800000',
              '#A52A2A', '#1F78B4', '#FF0000']
})

# Create an interval for each gene
genes_intervals = []
for idx, row in genes.iterrows():
    genes_intervals.append({
        'Gene': row['Gene'],
        'Start': row['Start'],
        'End': row['End'],
        'Color': row['Color']
    })

# Plot deletions colored according to genes in ax_del
for idx, row in df.iterrows():
    node_name = row['Node']
    
    # Change "nod_" to "node_" in node names
    if isinstance(node_name, str) and node_name.startswith("nod_"):
        node_name = node_name.replace("nod_", "node_", 1)
    
    if node_name in name_to_y:
        y = name_to_y[node_name]
        del_start = row['Deletion_Start']
        del_end = row['Deletion_End']
        state = row['State']
        # Get transparency according to state
        if state == 'ancestral':
            alpha = 0.2  # More transparent
        elif state == 'new':
            alpha = 1.0  # Solid color
        else:
            alpha = 1.0  # Default value
        # Find intersections with genes
        for gene in genes_intervals:
            gene_start = gene['Start']
            gene_end = gene['End']
            # Find overlap
            overlap_start = max(del_start, gene_start)
            overlap_end = min(del_end, gene_end)
            if overlap_start < overlap_end:
                color = gene['Color']
                ax_del.hlines(y, overlap_start, overlap_end, colors=color, linewidth=2, alpha=alpha)
    else:
        print(f"Node '{node_name}' not found in tree.")

# Now, identify deletions that completely remove genes
full_gene_deletions = []

for idx, row in df.iterrows():
    del_start = row['Deletion_Start']
    del_end = row['Deletion_End']
    node_name = row['Node']
    state = row['State']
    # Find genes that are completely removed by this deletion
    completely_deleted_genes = genes[
        (genes['Start'] >= del_start) & (genes['End'] <= del_end)
    ]
    if not completely_deleted_genes.empty:
        for idx2, gene_row in completely_deleted_genes.iterrows():
            full_gene_deletions.append({
                'Node': node_name,
                'State': state,
                'Del_Start': gene_row['Start'],
                'Del_End': gene_row['End'],
                'Gene': gene_row['Gene'],
                'Color': gene_row['Color']
            })

# Create a DataFrame with complete gene deletions
df_full_gene_dels = pd.DataFrame(full_gene_deletions)

# Plot complete gene deletions in ax_del_full_genes
for idx, row in df_full_gene_dels.iterrows():
    node_name = row['Node']
    
    # Change "nod_" to "node_" in node names
    if isinstance(node_name, str) and node_name.startswith("nod_"):
        node_name = node_name.replace("nod_", "node_", 1)
    
    if node_name in name_to_y:
        y = name_to_y[node_name]
        del_start = row['Del_Start']
        del_end = row['Del_End']
        color = row['Color']
        state = row['State']
        # Get transparency according to state
        if state == 'ancestral':
            alpha = 0.2  # More transparent
        elif state == 'new':
            alpha = 1.0  # Solid color
        else:
            alpha = 1.0  # Default value
        ax_del_full_genes.hlines(y, del_start, del_end, colors=color, linewidth=2, alpha=alpha)
    else:
        print(f"Node '{node_name}' not found in tree.")

# Adjust Y axis labels to include internal and terminal nodes
yticks = []
yticklabels = []
label_colors = []  # Para guardar colores de fondo según género

# Definir colores de fondo para cada género
genus_colors = {
    'A': '#FFA500',  # Naranja - Alphacoronavirus
    'B': '#800080',  # Morado - Betacoronavirus
    'D': '#FF0000',  # Rojo - Deltacoronavirus
    'G': '#008000',  # Verde - Gammacoronavirus
}

for clade in node_positions:
    if clade.name:
        y = node_positions[clade]
        yticks.append(y)
        # Change "nod_" to "node_" in labels
        label = clade.name
        if isinstance(label, str) and label.startswith("nod_"):
            label = label.replace("nod_", "node_", 1)
        yticklabels.append(label)
        
        # Determine background color by genus (only for species, not nodes)
        bg_color = None
        if isinstance(label, str) and '-' in label:
            # It's a species (terminal taxon)
            genus = label.split('-')[0].strip()
            bg_color = genus_colors.get(genus, None)
        label_colors.append(bg_color)

# Assign Y axis labels to both graphs with custom colors
ax_del.set_yticks(yticks)
ax_del.set_yticklabels(yticklabels, fontsize=8)

ax_del_full_genes.set_yticks(yticks)
ax_del_full_genes.set_yticklabels(yticklabels, fontsize=8)

# Apply background colors to labels in both axes
for i, (tick_label, bg_color) in enumerate(zip(ax_del.get_yticklabels(), label_colors)):
    if bg_color:
        tick_label.set_bbox(dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.3, edgecolor='none'))

for i, (tick_label, bg_color) in enumerate(zip(ax_del_full_genes.get_yticklabels(), label_colors)):
    if bg_color:
        tick_label.set_bbox(dict(boxstyle='round,pad=0.3', facecolor=bg_color, alpha=0.3, edgecolor='none'))

# Invert Y axis so nodes appear from top to bottom
ax_del.invert_yaxis()
ax_del_full_genes.invert_yaxis()

# Adjust labels and title
ax_del.set_xlabel('Genomic Position')
ax_del.set_title('Deletions along Sequence')

ax_del_full_genes.set_xlabel('Genomic Position')
ax_del_full_genes.set_title('Complete Gene Deletions')

# Integrate gene bar within ax_del
gene_bar_y = total_nodes + 0.5  # Y position for gene bar
for idx, row in genes.iterrows():
    x_start = row['Start']
    x_end = row['End']
    color = row['Color']
    ax_del.hlines(gene_bar_y, x_start, x_end, colors=color, linewidth=5)

# Integrate gene bar within ax_del_full_genes
for idx, row in genes.iterrows():
    x_start = row['Start']
    x_end = row['End']
    color = row['Color']
    ax_del_full_genes.hlines(gene_bar_y, x_start, x_end, colors=color, linewidth=5)

# Adjust Y axis range of ax_del and ax_del_full_genes to include gene bar
ax_del.set_ylim(-1, gene_bar_y + 1)
ax_del_full_genes.set_ylim(-1, gene_bar_y + 1)

# Invert Y axis again after adjusting range
ax_del.invert_yaxis()
ax_del_full_genes.invert_yaxis()

# Create legend axis
ax_legends = fig.add_axes([0.97, 0.25, 0.02, 0.65], frameon=False)
ax_legends.axis('off')

# Create legend elements for genes
genes_unicos = genes.drop_duplicates(subset=['Gen', 'Color'])
gene_legend_elements = [
    Patch(facecolor=row['Color'], edgecolor='none', label=row['Gen'])
    for idx, row in genes_unicos.iterrows()
]

# Add gene legend
legend_genes = ax_legends.legend(
    handles=gene_legend_elements,
    loc='upper left',
    fontsize=10,
    title='Genes',
    title_fontsize=12,
    frameon=False,
    ncol=1
)

# Create genus color legend at the top of the figure
genus_legend_elements = [
    Patch(facecolor='#FFA500', edgecolor='black', alpha=0.3, label='Alphacoronavirus (A-)'),
    Patch(facecolor='#800080', edgecolor='black', alpha=0.3, label='Betacoronavirus (B-)'),
    Patch(facecolor='#FF0000', edgecolor='black', alpha=0.3, label='Deltacoronavirus (D-)'),
    Patch(facecolor='#008000', edgecolor='black', alpha=0.3, label='Gammacoronavirus (G-)')
]

# Add genus legend at the top center of the figure (aligned with deletion panels)
genus_legend = fig.legend(
    handles=genus_legend_elements,
    loc='upper center',
    fontsize=10,
    title='Genera',
    title_fontsize=12,
    frameon=True,
    ncol=4,
    bbox_to_anchor=(0.60, 0.93),  # Moved to the right to center with deletion panels
    edgecolor='black',
    fancybox=True
)

# Save and show figure
plt.savefig('../figures/Figure_6.pdf', dpi=300, bbox_inches='tight')
plt.savefig('../figures/Figure_6.png', dpi=300, bbox_inches='tight')
print(f"✅ Full figure saved: Figure_6.pdf and Figure_6.png")
plt.close()

# ============================================================================
# GENERATE REDUCED FIGURE WITH ONLY GENUS-LEVEL NODES
# ============================================================================

print("\n🔬 Generating reduced figure with genus-level nodes only...")

# Identify genus common ancestor nodes by finding which node has all descendants of one genus
def get_genus_from_name(name):
    """Extract genus prefix from taxon name"""
    if isinstance(name, str) and '-' in name:
        return name.split('-')[0].strip()
    return None

def get_all_descendants(clade):
    """Get all terminal descendants of a clade"""
    terminals = []
    if clade.is_terminal():
        return [clade]
    for subclade in clade.clades:
        terminals.extend(get_all_descendants(subclade))
    return terminals

# Find the most recent common ancestor (MRCA) for each genus
genus_mrca = {}
for clade in node_positions:
    if not clade.is_terminal():
        # Get all terminal descendants
        descendants = get_all_descendants(clade)
        genera = set()
        for desc in descendants:
            genus = get_genus_from_name(desc.name)
            if genus:
                genera.add(genus)
        
        # If all descendants belong to the same genus, this could be the MRCA
        if len(genera) == 1:
            genus = list(genera)[0]
            # Keep the most ancestral node (the one with most descendants)
            if genus not in genus_mrca or len(descendants) > len(get_all_descendants(genus_mrca[genus])):
                genus_mrca[genus] = clade

print(f"✅ Identified {len(genus_mrca)} genus-level nodes: {list(genus_mrca.keys())}")

# Create mapping for reduced figure: genus name -> node info
# Map genus codes to full names
genus_full_names = {
    'A': 'Alphacoronavirus',
    'B': 'Betacoronavirus',
    'D': 'Deltacoronavirus',
    'G': 'Gammacoronavirus'
}

reduced_nodes = {}
for genus, clade in genus_mrca.items():
    # Convert node name
    node_name = clade.name if clade.name else f"node_{genus}"
    if isinstance(node_name, str) and node_name.startswith("nod_"):
        node_name = node_name.replace("nod_", "node_", 1)
    
    reduced_nodes[genus] = {
        'original_node': node_name,
        'y_position': node_positions[clade],
        'display_name': genus_full_names.get(genus, f"{genus}-Clade")  # Full genus name
    }

# Create reduced figure with larger dimensions
fig_reduced = plt.figure(figsize=(35, 20))  # Larger figure for better visibility

# Create axes for reduced figure - even more space between tree and graphs for names
ax_tree_r = fig_reduced.add_axes([0.005, 0.25, 0.05, 0.65])  # Tree even narrower
ax_del_r = fig_reduced.add_axes([0.22, 0.15, 0.35, 0.75])  # Deletion panel moved even further right for maximum space
ax_del_full_genes_r = fig_reduced.add_axes([0.61, 0.15, 0.35, 0.75], sharey=ax_del_r)  # Right panel adjusted

# This will be redrawn after Y positions are redistributed
# Placeholder for now
ax_tree_r.axis('off')

# Set up y positions and labels for reduced figure
# Redistribute Y positions evenly instead of using original positions
yticks_r = []
yticklabels_r = []
label_colors_r = []

# Sort genera by original position to maintain order
sorted_genera = sorted(reduced_nodes.keys(), key=lambda g: reduced_nodes[g]['y_position'])

# Create evenly spaced Y positions - more separated to match tree branch spacing
n_genera = len(sorted_genera)
if n_genera > 1:
    spacing = 30  # Increased spacing to match tree branch separation
    total_height = (n_genera - 1) * spacing
    start_y = 5  # Shift everything up
    
    # Assign new evenly-spaced Y positions
    new_y_positions = {}
    for i, genus in enumerate(sorted_genera):
        new_y = start_y + (i * spacing)
        new_y_positions[genus] = new_y
        yticks_r.append(new_y)
        yticklabels_r.append(reduced_nodes[genus]['display_name'])
        label_colors_r.append(genus_colors.get(genus, None))
    
    # Update reduced_nodes with new Y positions for plotting
    for genus, new_y in new_y_positions.items():
        reduced_nodes[genus]['plot_y_position'] = new_y
else:
    for genus in sorted_genera:
        yticks_r.append(5)  # Shift single genus up to match new spacing
        yticklabels_r.append(reduced_nodes[genus]['display_name'])
        label_colors_r.append(genus_colors.get(genus, None))
        reduced_nodes[genus]['plot_y_position'] = 5

# Debug: print labels
print(f"🔍 Y ticks for reduced figure: {yticks_r}")
print(f"🔍 Y labels for reduced figure: {yticklabels_r}")
print(f"🔍 Label colors: {label_colors_r}")

# Plot deletions for genus-level nodes
for idx, row in df.iterrows():
    node_name = row['Node']
    
    # Convert nod_ to node_
    if isinstance(node_name, str) and node_name.startswith("nod_"):
        node_name = node_name.replace("nod_", "node_", 1)
    
    # Check if this deletion belongs to any genus MRCA node
    for genus, node_info in reduced_nodes.items():
        if node_name == node_info['original_node']:
            y = node_info['plot_y_position']  # Use redistributed Y position
            del_start = row['Deletion_Start']
            del_end = row['Deletion_End']
            state = row['State']
            
            # Get transparency
            if state == 'ancestral':
                alpha = 0.2
            elif state == 'new':
                alpha = 1.0
            else:
                alpha = 1.0
            
            # Find intersections with genes for pointwise deletions
            for gene in genes_intervals:
                gene_start = gene['Inicio']
                gene_end = gene['Fin']
                overlap_start = max(del_start, gene_start)
                overlap_end = min(del_end, gene_end)
                if overlap_start < overlap_end:
                    color = gene['Color']
                    ax_del_r.hlines(y, overlap_start, overlap_end, colors=color, linewidth=12, alpha=alpha)  # Very thick lines

# Plot complete gene deletions for genus-level nodes
for idx, row in df_full_gene_dels.iterrows():
    node_name = row['Node']
    
    # Convert nod_ to node_
    if isinstance(node_name, str) and node_name.startswith("nod_"):
        node_name = node_name.replace("nod_", "node_", 1)
    
    # Check if this deletion belongs to any genus MRCA node
    for genus, node_info in reduced_nodes.items():
        if node_name == node_info['original_node']:
            y = node_info['plot_y_position']  # Use redistributed Y position
            del_start = row['Del_Start']
            del_end = row['Del_End']
            color = row['Color']
            state = row['State']
            
            if state == 'ancestral':
                alpha = 0.2
            elif state == 'new':
                alpha = 1.0
            else:
                alpha = 1.0
            
            ax_del_full_genes_r.hlines(y, del_start, del_end, colors=color, linewidth=12, alpha=alpha)  # Very thick lines

# Calculate gene bar Y position NOW (before drawing tree) so tree limits can match
gene_bar_y_r = max(yticks_r) + 15  # Increased separation and shifted up

# Now draw the tree with the redistributed Y positions showing real relationships
# Correct phylogeny: (G(D(A,B)))
if len(yticks_r) > 0:
    # Get Y positions for each genus
    y_by_genus = {genus: info['plot_y_position'] for genus, info in reduced_nodes.items()}
    
    # Tree topology: (G(D(A,B)))
    # This means: A and B are sisters, D is sister to (A,B), G is basal/outgroup
    
    y_min, y_max = min(yticks_r), max(yticks_r)
    
    # Define X positions for different levels
    x_root = 0
    x_level1 = 1.5  # G branches off here
    x_level2 = 3.0  # D branches off here  
    x_level3 = 4.5  # A and B branch off here
    x_terminal = 6.0  # Terminal tips
    
    if 'A' in y_by_genus and 'B' in y_by_genus and 'D' in y_by_genus and 'G' in y_by_genus:
        y_A = y_by_genus['A']
        y_B = y_by_genus['B']
        y_D = y_by_genus['D']
        y_G = y_by_genus['G']
        
        # Calculate internal node positions
        y_AB = (y_A + y_B) / 2  # A-B common ancestor
        y_DAB = (y_D + y_AB) / 2  # D-(A,B) common ancestor
        y_root = (y_G + y_DAB) / 2  # Root
        
        # Draw root
        ax_tree_r.plot(x_root, y_root, 'o', color='gray', markersize=8, zorder=1)
        
        # Root to main split
        ax_tree_r.plot([x_root, x_level1], [y_root, y_root], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot([x_level1, x_level1], [y_G, y_DAB], color='black', linewidth=3, zorder=1)
        
        # G branch (outgroup)
        ax_tree_r.plot([x_level1, x_terminal], [y_G, y_G], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot(x_terminal, y_G, 'o', color='black', markersize=12, zorder=2)
        
        # To D-(A,B) node
        ax_tree_r.plot([x_level1, x_level2], [y_DAB, y_DAB], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot([x_level2, x_level2], [y_D, y_AB], color='black', linewidth=3, zorder=1)
        
        # D branch
        ax_tree_r.plot([x_level2, x_terminal], [y_D, y_D], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot(x_terminal, y_D, 'o', color='black', markersize=12, zorder=2)
        
        # To A-B node
        ax_tree_r.plot([x_level2, x_level3], [y_AB, y_AB], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot([x_level3, x_level3], [y_A, y_B], color='black', linewidth=3, zorder=1)
        
        # A and B branches
        ax_tree_r.plot([x_level3, x_terminal], [y_A, y_A], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot(x_terminal, y_A, 'o', color='black', markersize=12, zorder=2)
        
        ax_tree_r.plot([x_level3, x_terminal], [y_B, y_B], color='black', linewidth=3, zorder=1)
        ax_tree_r.plot(x_terminal, y_B, 'o', color='black', markersize=12, zorder=2)
    
    # Set limits to match EXACTLY with the deletion panels - force same Y limits
    ax_tree_r.set_xlim(-2.0, x_terminal + 0.5)  # More space on the left for names
    # Force tree Y limits to match deletion panels exactly
    ax_tree_r.set_ylim(ax_del_r.get_ylim())  # Copy exact limits from deletion panel

# Set axis labels - VERY LARGE font with more separation
ax_del_r.set_xlabel('Genomic Position', fontsize=26, fontweight='bold', labelpad=20)  # Increased labelpad
ax_del_r.set_title('Deletions along Sequence', fontsize=28, fontweight='bold', pad=25)  # Increased pad

ax_del_full_genes_r.set_xlabel('Genomic Position', fontsize=26, fontweight='bold', labelpad=20)  # Increased labelpad
ax_del_full_genes_r.set_title('Complete Gene Deletions', fontsize=28, fontweight='bold', pad=25)  # Increased pad

# Much larger tick labels for x-axis
ax_del_r.tick_params(axis='x', labelsize=20)
ax_del_full_genes_r.tick_params(axis='x', labelsize=20)

# Invert Y axis
ax_del_r.invert_yaxis()
ax_del_full_genes_r.invert_yaxis()

# Add gene bar - thicker and more separated (gene_bar_y_r was already calculated earlier)
for idx, row in genes.iterrows():
    x_start = row['Start']
    x_end = row['End']
    color = row['Color']
    ax_del_r.hlines(gene_bar_y_r, x_start, x_end, colors=color, linewidth=16)  # Very thick gene bar
    ax_del_full_genes_r.hlines(gene_bar_y_r, x_start, x_end, colors=color, linewidth=16)  # Very thick gene bar

# First invert Y axis
ax_del_r.invert_yaxis()
ax_del_full_genes_r.invert_yaxis()

# Then adjust Y limits - adjusted for shifted positions
ax_del_r.set_ylim(gene_bar_y_r + 8, -2)  # Inverted limits (max, min) since axis is inverted
ax_del_full_genes_r.set_ylim(gene_bar_y_r + 8, -2)  # Inverted limits

# Expand X limits to the left to show genus names with negative positions
# Get the maximum genomic position from genes
max_genomic_pos = genes['End'].max()
ax_del_r.set_xlim(-5000, max_genomic_pos + 1000)  # Expand far to the left for genus names
ax_del_full_genes_r.set_xlim(-5000, max_genomic_pos + 1000)  # Expand far to the left

# Set ticks but NO labels (we'll add them manually with text)
ax_del_r.set_yticks(yticks_r)
ax_del_r.set_yticklabels([])  # Empty - will add with text
ax_del_full_genes_r.set_yticks(yticks_r)
ax_del_full_genes_r.set_yticklabels([])  # Empty labels

# Manually add genus names as text labels with colored backgrounds - positioned ultra far to the left
# Use transform to mix axis coordinates (x) with data coordinates (y)
from matplotlib.transforms import blended_transform_factory
trans = blended_transform_factory(ax_del_r.transAxes, ax_del_r.transData)

for i, (y_pos, label, bg_color) in enumerate(zip(yticks_r, yticklabels_r, label_colors_r)):
    # Position text in axis coordinates (0-1) for X, data coordinates for Y
    x_text_pos = -0.40  # Negative axis coordinates position (to the left of the plot area)
    
    # Add text with background color - VERY LARGE
    text_obj = ax_del_r.text(x_text_pos, y_pos, label, 
                             fontsize=28, fontweight='bold',
                             ha='left', va='center',
                             transform=trans,  # Use blended transform
                             bbox=dict(boxstyle='round,pad=1.0', 
                                     facecolor=bg_color if bg_color else 'white', 
                                     alpha=0.3,  # More transparent for lighter color
                                     edgecolor='black',
                                     linewidth=2) if bg_color else None,
                             zorder=100,
                             clip_on=False)  # Don't clip text outside axis
    
print(f"✅ Added {len(yticklabels_r)} genus labels manually")

# Add gene legend (but NOT genus legend - that info is in the node names) - closer to graphs
ax_legends_r = fig_reduced.add_axes([0.98, 0.15, 0.015, 0.75], frameon=False)  # Closer to graphs
ax_legends_r.axis('off')

legend_genes_r = ax_legends_r.legend(
    handles=gene_legend_elements,
    loc='upper left',
    fontsize=24,  # EVEN LARGER font for gene names
    title='Genes',
    title_fontsize=26,  # EVEN LARGER title
    frameon=True,
    ncol=1,
    edgecolor='black',
    fancybox=True
)

# Save reduced figure
plt.savefig('../figures/Figure_6_reduced.pdf', dpi=300, bbox_inches='tight')
plt.savefig('../figures/Figure_6_reduced.png', dpi=300, bbox_inches='tight')
print(f"✅ Reduced figure saved: Figure_6_reduced.pdf and Figure_6_reduced.png")
print(f"   Showing {len(reduced_nodes)} genus-level nodes instead of {len(yticks)} total nodes")

plt.close(fig_reduced)