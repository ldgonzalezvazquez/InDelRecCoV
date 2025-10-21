#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Python script to create beautiful recombination plots from LDhat results
Replaces the R scripts with a more elegant design including gene annotations
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from matplotlib.patches import Rectangle
import seaborn as sns

# Set style for beautiful plots
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

def define_gene_positions():
    """Define gene positions and colors using the exact same palette as the genome figure"""
    
    # Define genes with their positions - handle genes that appear in multiple regions
    genes = {
        '5\'UTR': [0, 1016],
        'ORF1ab': [1017, 30926],
        'NS2': [31842, 32687],
        'HE': [33115, 34450],
        'ORF2': [34466, 35821],
        'S': [35971, 43557],
        '4.9/4.8kDaNSP': [45085, 45389],
        'NSP3/ORF3': [45469, 46925],
        'NS4/15kDaNSP': [47014, 47758],
        'ORF4/5': [47932, 49740],
        'E': [49820, 50176],
        'M': [50246, 51778],
        'NS6': [52427, 52753],
        'ORFx_1': [52938, 53222],  # First ORFx
        'ORF5': [53302, 54772],
        'ORF6': [54776, 55480],
        'ORF7_1': [55488, 56167],  # First ORF7
        'ORF8_1': [56320, 56817],  # First ORF8
        'ORF9_1': [56892, 57350],  # First ORF9
        'ORF10_1': [57353, 57985], # First ORF10
        'N2': [58097, 58789],
        'N': [58790, 60911],
        'IP': [61550, 62218],
        'p10': [62905, 63180],
        'ORF7/3x-lp': [63198, 65642],
        'ORF8_2': [65975, 66637],  # Second ORF8
        'ORF9_2': [67783, 68080],  # Second ORF9
        'ORF10_2': [68236, 68531], # Second ORF10
        'ORF11': [68532, 69074],
        'ORF': [69075, 69287],
        'ORFx_2': [69288, 69850],  # Second ORFx
        'ORFy': [69851, 70384],
        '3\'UTR': [70385, 71920]
    }
    
    # Use the EXACT same color palette as the genome figure
    custom_colors = [
        '#FF0000', '#008000', '#FFA500', '#0000FF', '#FFD700', '#800080', '#00FA9A', '#FF69B4', '#808000',
        '#6495ED', '#D2691E', '#5F9EA0', '#FF00FF', '#A52A2A', '#7FFF00', '#8A2BE2', '#F08080', '#4B0082',
        '#FF7F50', '#008080', '#CD5C5C', '#ADFF2F', '#A6CEE3', '#33A02C', '#DEB887', '#000080', '#800000', '#1F78B4', '#FF0000'
    ]
    
    # Create a list of unique base gene names (without _1, _2 suffixes)
    base_gene_names = []
    for gene_name in genes.keys():
        base_name = gene_name.split('_')[0] if '_' in gene_name else gene_name
        if base_name not in base_gene_names:
            base_gene_names.append(base_name)
    
    # Assign colors to base gene names first (like in the genome figure)
    base_gene_colors = {}
    for i, base_name in enumerate(base_gene_names):
        color_index = i % len(custom_colors)
        base_gene_colors[base_name] = custom_colors[color_index]
    
    # Now assign colors to all gene instances
    gene_colors = {}
    for gene_name in genes.keys():
        base_name = gene_name.split('_')[0] if '_' in gene_name else gene_name
        gene_colors[gene_name] = base_gene_colors[base_name]
    
    # Make UTRs the same color (first color from palette)
    utr_color = custom_colors[0]
    gene_colors['5\'UTR'] = utr_color
    gene_colors['3\'UTR'] = utr_color
    
    return genes, gene_colors

def load_recombination_data(filename='../data/res_full.txt'):
    """Load recombination data from LDhat results"""
    try:
        # Read the data
        data = pd.read_csv(filename, sep='\t')
        print(f"✅ Loaded {len(data)} positions from {filename}")
        return data
    except Exception as e:
        print(f"❌ Error loading {filename}: {e}")
        return None

def create_gene_bar(ax, genes, gene_colors, y_position=0):
    """Create a gene annotation bar below the main plot with unique colors"""
    
    bar_height = 0.15  # Reduced from 0.3 to make bars thinner
    y_bottom = y_position - bar_height/2
    y_top = y_position + bar_height/2
    
    print(f"🎨 Drawing {len(genes)} genes in the bar...")
    
    # Draw genes
    for gene_name, (start, end) in genes.items():
        color = gene_colors[gene_name]
        print(f"   • {gene_name}: {start:,} - {end:,} | Color: {color}")
        
        # Create rectangle for gene
        rect = Rectangle((start, y_bottom), end-start, bar_height, 
                        facecolor=color, edgecolor='black', linewidth=0.2, alpha=0.8)
        ax.add_patch(rect)
        
        # NO MORE GENE LABELS ON TOP - only in legend
    
    print(f"✅ Gene bar completed")

def create_recombination_plot(data, genes, gene_colors, output_prefix='recombination_plot'):
    """Create the main recombination plot with gene annotations"""
    
        # Create figure with subplots - increased height for better spacing
    fig, (ax_main, ax_genes) = plt.subplots(2, 1, figsize=(20, 14), 
                                             gridspec_kw={'height_ratios': [8, 1]},
                                             sharex=True)
    
    # Main recombination plot
    positions = data['Loci'].values
    mean_rho = data['Mean_rho'].values
    l95 = data['L95'].values
    u95 = data['U95'].values
    
    # Plot confidence intervals as filled area
    ax_main.fill_between(positions, l95, u95, alpha=0.3, color='lightblue', 
                         label='95% Confidence Interval')
    
    # Plot mean recombination rate
    ax_main.plot(positions, mean_rho, color='darkblue', linewidth=1, 
                 label='Mean of Population Recombination Rate (ρ)')
    
    # Customize main plot - NO TITLE
    ax_main.set_ylabel('Population Recombination Rate (ρ)', fontsize=24, fontweight='bold')
    ax_main.grid(True, alpha=0.3)
    ax_main.legend(fontsize=22, loc='upper right')
    ax_main.tick_params(axis='both', which='major', labelsize=22)
    
    # Set x-axis limits
    ax_main.set_xlim(0, max(positions))
    
    # Gene annotation bar
    create_gene_bar(ax_genes, genes, gene_colors)
    
    # Customize gene bar
    ax_genes.set_xlabel('Genome Position (bp)', fontsize=24, fontweight='bold', labelpad=15)
    ax_genes.set_ylabel('Genes', fontsize=24, fontweight='bold')
    ax_genes.set_ylim(-0.2, 0.2)  # Adjust vertical space for gene bar
    ax_genes.set_xlim(0, max(positions))
    ax_genes.set_yticks([0])
    ax_genes.set_yticklabels([''])
    ax_genes.tick_params(axis='x', which='major', labelsize=22)
    
    # Add grid for gene positions
    ax_genes.grid(True, alpha=0.2, axis='x')
    
    # Create a compact multi-column legend at the top - handle duplicate genes
    legend_elements = []
    seen_genes = set()
    
    for gene_name, color in gene_colors.items():
        # Extract base gene name (remove _1, _2 suffixes)
        base_name = gene_name.split('_')[0] if '_' in gene_name else gene_name
        
        # Only add to legend if we haven't seen this base gene name before
        if base_name not in seen_genes:
            legend_elements.append(mpatches.Patch(color=color, label=base_name))
            seen_genes.add(base_name)
    
    # Create a separate legend box above the main plot area - extended to full width
    fig.legend(handles=legend_elements, loc='upper center', 
              title='Gene Annotations', fontsize=20, title_fontsize=24, # Adjusted font sizes
              ncol=len(legend_elements)//3, bbox_to_anchor=(0.52, 0.92), frameon=True, 
              columnspacing=1.0, handletextpad=0.5) # Extended legend across plot width
    
        # Adjust layout to prevent legend overlap and make space for external legend
    plt.tight_layout()
    
    # Additional adjustment for legend spacing - make room for legend above
    plt.subplots_adjust(top=0.75)  # Leave space above for legend
    
        # Save plots
    pdf_name = "../figures/Figure_3.pdf"
    png_name = "../figures/Figure_3.png"
    
    plt.savefig(pdf_name, dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(png_name, dpi=300, bbox_inches='tight', facecolor='white')
    
    print(f"🎉 Plots saved: {pdf_name} and {png_name}")
    
    plt.close()  # Close instead of show to prevent interactive window
    
    return fig

def create_summary_statistics(data):
    """Create summary statistics for the recombination analysis"""
    
    print("\n📊 RECOMBINATION ANALYSIS SUMMARY")
    print("=" * 50)
    
    # Basic statistics
    total_positions = len(data)
    polymorphic_positions = len(data[data['Mean_rho'] > 0])
    non_polymorphic_positions = total_positions - polymorphic_positions
    
    print(f"Total genome positions: {total_positions:,}")
    print(f"Polymorphic sites: {polymorphic_positions:,}")
    print(f"Non-polymorphic sites: {non_polymorphic_positions:,}")
    
    # Recombination statistics
    if polymorphic_positions > 0:
        mean_rho_poly = data[data['Mean_rho'] > 0]['Mean_rho'].mean()
        max_rho = data['Mean_rho'].max()
        min_rho = data[data['Mean_rho'] > 0]['Mean_rho'].min()
        
        print(f"\nRecombination Rate Statistics:")
        print(f"Mean ρ (polymorphic sites): {mean_rho_poly:.4f}")
        print(f"Maximum ρ: {max_rho:.4f}")
        print(f"Minimum ρ (polymorphic): {min_rho:.4f}")
    
    # Genome coverage
    genome_length = data['Loci'].max()
    print(f"\nGenome length: {genome_length:,} bp")
    
    return {
        'total_positions': total_positions,
        'polymorphic_positions': polymorphic_positions,
        'genome_length': genome_length
    }

def main():
    """Main function to run the recombination analysis"""
    
    print("🧬 CORONAVIRUS RECOMBINATION ANALYSIS")
    print("=" * 50)
    
    # Define gene positions and colors
    genes, gene_colors = define_gene_positions()
    print(f"✅ Defined {len(genes)} genes with unique colors")
    
    # Load recombination data
    data = load_recombination_data()
    if data is None:
        return
    
    # Create summary statistics
    stats = create_summary_statistics(data)
    
    # Create the main plot
    print("\n🎨 Creating recombination plot...")
    fig = create_recombination_plot(data, genes, gene_colors, 
                                  output_prefix='coronavirus_recombination_analysis')
    
    print("\n✅ Analysis completed successfully!")
    print("📁 Check the generated PDF and PNG files")

if __name__ == "__main__":
    main()
