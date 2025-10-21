#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simple script to generate phylogenetic tree with bootstrap values and colors
"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from Bio import Phylo
import os

def load_tree_with_bootstrap():
    """Load tree with bootstrap values from pre-computed file"""
    
    bootstrap_tree_file = '../data/treerooted_with_bootstrap.nwk'
    
    if not os.path.exists(bootstrap_tree_file):
        print("❌ Bootstrap tree file not found!")
        print(f"   Expected: {bootstrap_tree_file}")
        print("\n📋 Please run tree_apply_bootstrap.py first:")
        print("   cd data")
        print("   python3 ../scripts/tree_apply_bootstrap.py \\")
        print("     --ref treerooted.nwk \\")
        print("     --boot RAxML_tree.raxml.bootstraps \\")
        print("     --aliases names_correpondence.txt \\")
        print("     --out treerooted_with_bootstrap.nwk \\")
        print("     --log bootstrap_mapping.log")
        raise FileNotFoundError(f"Bootstrap tree file not found: {bootstrap_tree_file}")
    
    print("🔧 Loading tree with bootstrap values...")
    tree = Phylo.read(bootstrap_tree_file, 'newick')
    
    # Convert confidence from string (internal node labels) to numeric confidence values
    for node in tree.get_nonterminals():
        if node.name and node.name.isdigit():
            node.confidence = int(node.name)
            node.name = None  # Remove the name after converting to confidence
    
    internal_nodes = list(tree.get_nonterminals())
    nodes_with_support = sum(1 for node in internal_nodes if node.confidence is not None)
    print(f"✅ Tree loaded with {nodes_with_support}/{len(internal_nodes)} nodes with bootstrap support")
    
    return tree

def define_genus_colors():
    """Define colors for each genus"""
    return {
        'A': '#FF8C00',  # Orange - Alphacoronavirus
        'B': '#8B008B',  # Purple - Betacoronavirus
        'D': '#FF0000',  # Red - Deltacoronavirus
        'G': '#008000',  # Green - Gammacoronavirus
        'U': '#000000',  # Black - Unclassified
    }

def draw_tree(tree):
    """Draws the tree with colors and bootstrap - clean design in English"""
    
    # Define colors
    genus_colors = define_genus_colors()
    
    # Create figure with extra width for legends on the right
    fig, ax = plt.subplots(figsize=(50, 60))
    
    # Get genera
    genera_used = set()
    for clade in tree.get_terminals():
        if clade.name and '-' in clade.name:
            genus = clade.name.split('-')[0]
            genera_used.add(genus)
    
    print(f"🎨 Genera found: {sorted(genera_used)}")
    
    # Draw tree - only terminal taxon names, not internal nodes
    # Show ALL bootstrap values, even 0
    def get_branch_label(clade):
        if clade.confidence is not None:
            return f"{int(clade.confidence)}"
        return ""
    
    Phylo.draw(tree, 
              axes=ax,
              label_func=lambda x: x.name if x.is_terminal() else "",
              do_show=False,
              branch_labels=get_branch_label,
              show_confidence=False)
    
    # Color taxa and bootstrap
    texts_processed = 0
    bootstrap_shown = 0
    
    for text in ax.texts:
        if hasattr(text, 'get_text'):
            label_text = text.get_text()
            
            # If it's a taxon name (must contain a hyphen)
            if '-' in label_text and not label_text.replace('.', '').isdigit():
                # Clean spaces at beginning and end
                label_text_clean = label_text.strip()
                genus = label_text_clean.split('-')[0]
                color = genus_colors.get(genus, '#CCCCCC')
                text.set_color(color)
                text.set_fontweight('bold')
                text.set_fontsize(32)
                texts_processed += 1
            
            # If it's a bootstrap value (any digit, including single digits like 0-9)
            elif label_text.replace('.', '').replace('-', '').isdigit():
                try:
                    value = int(label_text)
                    # Color all bootstrap values with a gradient scale (0-100)
                    if value >= 90:
                        color_bootstrap = 'darkgreen'
                        bg_color = 'lightgreen'
                    elif value >= 75:
                        color_bootstrap = 'darkblue'
                        bg_color = 'lightblue'
                    elif value >= 50:
                        color_bootstrap = 'darkorange'
                        bg_color = 'lightyellow'
                    elif value >= 25:
                        color_bootstrap = 'darkred'
                        bg_color = 'lightcoral'
                    else:  # 0-24
                        color_bootstrap = 'purple'
                        bg_color = 'lavender'
                    
                    text.set_color(color_bootstrap)
                    text.set_fontweight('bold')
                    text.set_fontsize(26)
                    text.set_bbox(dict(boxstyle="round,pad=0.5", facecolor=bg_color, alpha=0.8, edgecolor=color_bootstrap))
                    bootstrap_shown += 1
                except ValueError:
                    pass  # Not a valid number, skip
    
    print(f"📝 Taxa colored: {texts_processed}")
    print(f"🔢 Bootstrap values shown: {bootstrap_shown}")
    
    # No title
    
    # Create separate legends for genera and bootstrap
    legend_genera = []
    legend_bootstrap = []
    
    # Genera
    genus_names = {
        'A': 'Alphacoronavirus',
        'B': 'Betacoronavirus', 
        'D': 'Deltacoronavirus',
        'G': 'Gammacoronavirus',
        'U': 'Unclassified'
    }
    
    for genus in sorted(genera_used):
        color = genus_colors.get(genus, '#CCCCCC')
        name = genus_names.get(genus, f'{genus}-coronavirus')
        legend_genera.append(mpatches.Patch(color=color, label=f'{genus}: {name}'))
    
    # Bootstrap with correct colors (colored background as in tree)
    legend_bootstrap.append(mpatches.Patch(facecolor='lightgreen', edgecolor='darkgreen', linewidth=2, label='Bootstrap ≥90'))
    legend_bootstrap.append(mpatches.Patch(facecolor='lightblue', edgecolor='darkblue', linewidth=2, label='Bootstrap 75-89'))
    legend_bootstrap.append(mpatches.Patch(facecolor='lightyellow', edgecolor='darkorange', linewidth=2, label='Bootstrap 50-74'))
    legend_bootstrap.append(mpatches.Patch(facecolor='lightcoral', edgecolor='darkred', linewidth=2, label='Bootstrap 25-49'))
    legend_bootstrap.append(mpatches.Patch(facecolor='lavender', edgecolor='purple', linewidth=2, label='Bootstrap <25'))
    
    # Get figure for handling legends
    fig = ax.figure
    
    # Genera legend aligned to the right (top)
    legend_genera_obj = fig.legend(handles=legend_genera, loc='upper right', 
                                   bbox_to_anchor=(0.98, 0.95), 
                                   title="Genera", 
                                   fontsize=30, title_fontsize=32,
                                   handlelength=2, handleheight=1.5, 
                                   frameon=True, edgecolor='black', fancybox=False)
    
    # Bootstrap legend aligned to the right (below the previous one)
    legend_bootstrap_obj = fig.legend(handles=legend_bootstrap, loc='upper right', 
                                     bbox_to_anchor=(0.98, 0.65), 
                                     title="Bootstrap Support", 
                                     fontsize=30, title_fontsize=32,
                                     handlelength=2, handleheight=1.5, 
                                     frameon=True, edgecolor='black', fancybox=False)
    
    # Add both legends to the figure
    ax.add_artist(legend_genera_obj)
    ax.add_artist(legend_bootstrap_obj)
    
    # Hide axes and ticks
    ax.axis('off')
    ax.set_xticks([])
    ax.set_yticks([])
    
    # Adjust layout
    plt.tight_layout()
    
    # Save
    pdf_name = '../figures/Figure_1.pdf'
    png_name = '../figures/Figure_1.png'
    
    plt.savefig(pdf_name, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(png_name, dpi=300, bbox_inches='tight', facecolor='white')
    
    print(f"🎉 Clean tree generated: {pdf_name} and {png_name}")
    
    plt.close()
    return pdf_name

def main():
    print("🌳 === PHYLOGENETIC TREE GENERATOR WITH BOOTSTRAP ===\n")
    
    # Load tree with bootstrap values
    tree = load_tree_with_bootstrap()
    
    # Draw tree
    print("\n🎨 Drawing tree...")
    file = draw_tree(tree)
    
    if file:
        print(f"\n✅ === PROCESS COMPLETED ===")
        print(f"🎨 Tree with colors by genus")
        print(f"📊 Bootstrap values visible")
        print(f"🌳 Complete phylogenetic structure")
        print(f"📋 Explanatory legend")
        print(f"🖼️  Saved in PDF and PNG")
    else:
        print("❌ Could not generate tree")

if __name__ == "__main__":
    main()
