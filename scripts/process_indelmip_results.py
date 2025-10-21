#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Process InDelMIP output to identify and classify deletion events
Analyzes binary sequences to detect ancestral vs. lineage-specific deletions
"""

import os
from Bio import SeqIO
from ete3 import Tree
from progress.bar import Bar

# Function to select file
def select_file(extensions, description):
    """Shows a list of files with specified extensions and allows selection."""
    files = [f for f in os.listdir() if any(f.endswith(ext) for ext in extensions)]
    if not files:
        print(f"No {description} files found.")
        return None

    print(f"\nSelect a {description} file to process:")
    for i, file in enumerate(files):
        print(f"{i + 1}: {file}")

    while True:
        try:
            selection = int(input(f"Enter the number of the {description} file to use: "))
            if 1 <= selection <= len(files):
                return files[selection - 1]
            else:
                print("Invalid number. Please try again.")
        except ValueError:
            print("Invalid input. Please enter a number.")

# Use predefined files from data/ directory
fasta_file = "../data/outputmip.fasta"
tree_file = "../data/treerooted.nwk"

# Verify files exist
if not os.path.exists(fasta_file) or not os.path.exists(tree_file):
    print(f"Error: Required files not found in ../data/")
    print(f"  - {fasta_file}")
    print(f"  - {tree_file}")
    exit()

# Function to detect deletions in a sequence
def detect_deletions(seq):
    """Detects continuous regions of 0s (deletions) in a sequence."""
    deletions = []
    start = None
    for i, char in enumerate(seq):
        if char == '0' and start is None:
            start = i
        elif char == '1' and start is not None:
            deletions.append((start + 1, i))
            start = None
    if start is not None:
        deletions.append((start + 1, len(seq)))
    return deletions

# Function to determine if a deletion is similar to one in the ancestral node
def is_similar_deletion(start, end, ancestor_start, ancestor_end, error_tolerance_percent):
    """Determines if a deletion is similar to one in the ancestral node, considering error percentage."""
    size = end - start + 1
    ancestor_size = ancestor_end - ancestor_start + 1
    event_mean = (start + end) / 2
    ancestor_mean = (ancestor_start + ancestor_end) / 2
    tolerance = size * (error_tolerance_percent / 100)

    if abs(event_mean - ancestor_mean) <= tolerance and abs(size - ancestor_size) <= tolerance:
        return True
    return False

# Function to preprocess sequence and merge small regions of '1's within '0's
def preprocess_sequence(seq, max_one_run_error_percent):
    """Replaces small regions of '1's surrounded by '0's with '0's."""
    seq_list = list(seq)
    n = len(seq_list)

    i = 0
    while i < n:
        if seq_list[i] == '0':
            i += 1
        else:
            start_one = i
            while i < n and seq_list[i] == '1':
                i += 1
            end_one = i - 1
            # Check if there are '0's to the left and right
            left = start_one - 1 >= 0 and seq_list[start_one - 1] == '0'
            right = end_one + 1 < n and seq_list[end_one + 1] == '0'
            if left and right:
                # Calculate size of adjacent '0's
                # Size of '0's to the left
                left_zero_start = start_one - 1
                while left_zero_start >= 0 and seq_list[left_zero_start] == '0':
                    left_zero_start -= 1
                left_zero_size = start_one - left_zero_start - 1
                # Size of '0's to the right
                right_zero_end = end_one + 1
                while right_zero_end < n and seq_list[right_zero_end] == '0':
                    right_zero_end += 1
                right_zero_size = right_zero_end - end_one - 1
                # Minimum size of adjacent '0's
                min_zero_size = left_zero_size + right_zero_size
                # Calculate maximum allowed size of '1's region
                max_one_run_size = int(min_zero_size * (max_one_run_error_percent / 100))
                if (end_one - start_one + 1) <= max_one_run_size:
                    # Replace '1's with '0's
                    for j in range(start_one, end_one + 1):
                        seq_list[j] = '0'
    return ''.join(seq_list)

# Load sequences
print("Loading FASTA sequences...")
sequences = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}

# Load tree from specified file
tree = Tree(tree_file, format=1)

# *** Step 1: Generate graph of independent events vs. preprocessing error percentage ***

# Function to count independent events with variable preprocessing
def count_independent_deletions_preprocess(tree, sequences, max_one_run_error_percent):
    """Counts independent deletions in each tree node varying preprocessing."""
    independent_deletions = 0
    # Preprocess sequences with given percentage
    preprocessed_sequences = {}
    for seq_name, seq in sequences.items():
        preprocessed_sequences[seq_name] = preprocess_sequence(seq, max_one_run_error_percent)
    nodes = list(tree.traverse("postorder"))
    for node in nodes:
        # Root node is now included
        node_seq = preprocessed_sequences.get(node.name)
        ancestor_seq = preprocessed_sequences.get(node.up.name) if node.up else ''
        if not node_seq:
            continue
        node_deletions = detect_deletions(node_seq)
        independent_deletions += len(node_deletions)
    return independent_deletions

# Generate graph of independent events vs. preprocessing error percentage
error_values_preprocess = list(range(0, 101, 10))
independent_deletions_counts_preprocess = []

print("\nCalculating independent events for different preprocessing error percentages...")
bar = Bar('Calculating', max=len(error_values_preprocess))
for max_one_run_error_percent in error_values_preprocess:
    count = count_independent_deletions_preprocess(tree, sequences, max_one_run_error_percent)
    independent_deletions_counts_preprocess.append(count)
    bar.next()
bar.finish()

# Generate text graph for preprocessing results
max_count_preprocess = max(independent_deletions_counts_preprocess)
print("\nGraph: Independent Deletion Events vs. Preprocessing Error Percentage\n")
print("Preprocessing Error Percentage | Number of Independent Deletion Events")
print("-" * 80)
for error, count in zip(error_values_preprocess, independent_deletions_counts_preprocess):
    bar_graph = "*" * int((count / max_count_preprocess) * 50) if max_count_preprocess > 0 else ""
    print(f"{error:>34}% | {bar_graph} ({count})")
print("\n")

# *** Step 2: Request preprocessing error percentage from user ***
print("The preprocessing error percentage determines how large regions of '1's within deletions can be to be ignored.")
print("A higher percentage will allow ignoring larger regions of '1's within deletions.")
max_one_run_error_percent = float(input("Enter the selected preprocessing error percentage (%): "))

# *** Step 3: Apply preprocessing with selected percentage ***
# Preprocess sequences and store them
print("\nApplying preprocessing to sequences...")
preprocessed_sequences = {}
for seq_name, seq in sequences.items():
    preprocessed_sequences[seq_name] = preprocess_sequence(seq, max_one_run_error_percent)

# *** Step 4: Generate graph of independent events vs. error tolerance percentage using preprocessed sequences ***

# Function to count independent events varying error tolerance percentage
def count_independent_deletions(tree, sequences, error_tolerance_percent):
    """Counts independent deletions in each tree node using preprocessed sequences."""
    independent_deletions = 0
    nodes = list(tree.traverse("postorder"))
    for node in nodes:
        # Root node is now included
        node_seq = sequences.get(node.name)
        ancestor_seq = sequences.get(node.up.name) if node.up else ''
        if not node_seq:
            continue
        node_deletions = detect_deletions(node_seq)
        ancestor_deletions = detect_deletions(ancestor_seq)
        for start, end in node_deletions:
            similar_event_found = False
            if ancestor_seq:
                for ancestor_start, ancestor_end in ancestor_deletions:
                    if is_similar_deletion(start, end, ancestor_start, ancestor_end, error_tolerance_percent):
                        similar_event_found = True
                        break
            if not similar_event_found:
                independent_deletions += 1
    return independent_deletions

# Generate graph of deletion events vs. error tolerance percentage
error_values = list(range(0, 101, 10))
independent_deletions_counts = []

print("\nCalculating independent events for different error tolerance percentages (with preprocessing applied)...")
bar = Bar('Calculating', max=len(error_values))
for error_tolerance_percent in error_values:
    count = count_independent_deletions(tree, preprocessed_sequences, error_tolerance_percent)
    independent_deletions_counts.append(count)
    bar.next()
bar.finish()

# Generate text graph for results
max_count = max(independent_deletions_counts)
print("\nGraph: Independent Deletion Events vs. Error Tolerance Percentage (With Preprocessing Applied)\n")
print("Error Tolerance Percentage | Number of Independent Deletion Events")
print("-" * 70)
for error, count in zip(error_values, independent_deletions_counts):
    bar_graph = "*" * int((count / max_count) * 50) if max_count > 0 else ""
    print(f"{error:>26}% | {bar_graph} ({count})")
print("\n")

# *** Step 5: Request error tolerance percentage from user ***
print("The error tolerance percentage determines how similar two deletions must be to be considered the same.")
print("A higher percentage will allow considering deletions with larger differences as similar.")
print("Both the start and end positions of deletions are considered, as well as deviation to either side.")
error_tolerance_percent = float(input("Enter the selected error tolerance percentage (%): "))

# *** Step 6: Generate graph of event sizes using preprocessed sequences and selected error_tolerance_percent ***
# Function to count and obtain sizes of independent events
def count_event_sizes_preprocessed(tree, sequences, error_tolerance_percent):
    """Counts sizes of independent events detected with preprocessing."""
    event_sizes = []
    nodes = list(tree.traverse("postorder"))
    for node in nodes:
        # Root node is now included
        node_seq = sequences.get(node.name)
        ancestor_seq = sequences.get(node.up.name) if node.up else ''
        if not node_seq:
            continue
        node_deletions = detect_deletions(node_seq)
        ancestor_deletions = detect_deletions(ancestor_seq)
        for start, end in node_deletions:
            similar_event_found = False
            if ancestor_seq:
                for ancestor_start, ancestor_end in ancestor_deletions:
                    if is_similar_deletion(start, end, ancestor_start, ancestor_end, error_tolerance_percent):
                        similar_event_found = True
                        break
            if not similar_event_found:
                event_sizes.append(end - start + 1)
    return event_sizes

# Get event sizes with selected preprocessing and tolerance percentage
event_sizes_preprocessed = count_event_sizes_preprocessed(tree, preprocessed_sequences, error_tolerance_percent)

# Define size ranges for graph
size_range = 10  # Size of each range interval
max_event_size = max(event_sizes_preprocessed) if event_sizes_preprocessed else 0
ranges = range(0, max_event_size + size_range, size_range)
size_counts_preprocessed = {f"{r}-{r + size_range - 1}": 0 for r in ranges}

# Count events in each range
for size in event_sizes_preprocessed:
    range_val = (size // size_range) * size_range
    range_key = f"{range_val}-{range_val + size_range - 1}"
    size_counts_preprocessed[range_key] += 1

# Generate size range graph with preprocessing
print("\nGraph: Deletion Event Sizes by Ranges (With Selected Preprocessing and Tolerance)\n")
print("Size Range | Frequency")
print("-" * 40)
max_freq = max(size_counts_preprocessed.values()) if size_counts_preprocessed else 0
for range_key, freq in size_counts_preprocessed.items():
    bar_graph = "*" * int((freq / max_freq) * 50) if max_freq > 0 else ""
    print(f"{range_key:>15} | {bar_graph} ({freq})")

# Step 7: Request minimum deletion size
print("\nThe minimum deletion size allows filtering small events that may not be relevant.")
min_event_size = int(input("Enter the minimum deletion size to consider for output files: "))

# *** Step 7: Execute event writing with selected parameters and preprocessed sequences ***
def check_ancestral_deletions(tree, sequences, error_tolerance_percent, min_event_size):
    """Writes detected events to output files."""
    with open("all_events.txt", mode="w") as all_file, open("new_events.txt", mode="w") as new_file:
        all_file.write("Node\tDeletion_Start\tDeletion_End\tSize\tState\n")
        new_file.write("Node\tDeletion_Start\tDeletion_End\tSize\tState\n")

        print("Writing events to files...")
        nodes = list(tree.traverse("postorder"))
        bar = Bar('Writing nodes', max=len(nodes))
        for node in nodes:
            # Root node is now included
            node_seq = sequences.get(node.name)
            ancestor_seq = sequences.get(node.up.name) if node.up else ''
            if not node_seq:
                bar.next()
                continue
            node_deletions = detect_deletions(node_seq)
            ancestor_deletions = detect_deletions(ancestor_seq)
            for start, end in node_deletions:
                size = end - start + 1
                if size < min_event_size:
                    continue
                if ancestor_seq:
                    state = "New" if not any(is_similar_deletion(start, end, ancestor_start, ancestor_end, error_tolerance_percent) for ancestor_start, ancestor_end in ancestor_deletions) else "Ancestral"
                else:
                    state = "New"
                all_file.write(f"{node.name}\t{start}\t{end}\t{size}\t{state}\n")
                if "New" in state:
                    new_file.write(f"{node.name}\t{start}\t{end}\t{size}\t{state}\n")
            bar.next()
        bar.finish()
    print("\nResults have been written to 'all_events.txt' and 'new_events.txt'.")

# Execute event writing
check_ancestral_deletions(tree, preprocessed_sequences, error_tolerance_percent, min_event_size)
