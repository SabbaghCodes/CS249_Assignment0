import argparse
import gc
import time
import tracemalloc
import sys

def read_fasta_pattern(file_path):
    """
    Reads the pattern from a FASTA file, removes headers, and returns the concatenated sequence in uppercase.
    
    Arguments:
        file_path (str): The path to the FASTA file containing the pattern.
        
    Returns:
        str: The concatenated pattern sequence in uppercase.
    """
    pattern = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):  # Skip header lines
                pattern.append(line.upper())
    return ''.join(pattern)

def compute_lps(pattern):
    """
    Computes the LPS (Longest Prefix Suffix) array used in the KMP algorithm for efficient pattern matching.
    
    Arguments:
        pattern (str): The pattern to be searched.
        
    Returns:
        list: The LPS array representing the longest prefix suffix for the pattern.
    """
    m = len(pattern)
    lps = [0] * m
    j = 0  # length of previous lps
    for i in range(1, m):
        while j > 0 and pattern[i] != pattern[j]:
            j = lps[j - 1]
        if pattern[i] == pattern[j]:
            j += 1
            lps[i] = j
    return lps

def log_message(f, message):
    """
    Logs a message to both the console and an output file.
    
    Arguments:
        f (file object): The output file to log the message.
        message (str): The message to be logged.
    """
    print(message)
    f.write(message + "\n")

def search_exact_matches(f, chromosome, genome_sequence, pattern):
    """
    Searches for exact matches of the pattern in the genome sequence using the KMP algorithm.
    
    Arguments:
        f (file object): The output file to log the matches.
        chromosome (str): The chromosome name where the search is performed.
        genome_sequence (str): The genome sequence to search within.
        pattern (str): The pattern to search for.
        
    Returns:
        int: The total number of exact matches found.
    """
    total_matches = 0
    lps = compute_lps(pattern)  # Compute LPS array for pattern
    pattern_len = len(pattern)
    j = 0  # KMP state
    
    # KMP-based search
    for i in range(len(genome_sequence)):
        while j > 0 and genome_sequence[i] != pattern[j]:
            j = lps[j - 1]
        if genome_sequence[i] == pattern[j]:
            j += 1
        else:
            j = 0
        if j == pattern_len:  # Match found
            start = i - pattern_len + 1
            end = start + pattern_len - 1
            log_message(f, f"{chromosome}\t{start}\t{end}\tExact match\n")
            total_matches += 1
            j = lps[j - 1]
    
    return total_matches

def count_exact_matches(fasta_file, pattern, output_file):
    """
    Counts the total number of exact matches of the pattern in the genome sequence.
    
    Arguments:
        fasta_file (str): The path to the genome FASTA file.
        pattern (str): The pattern to search for.
        output_file (str): The file to store the results.
    """
    start_time = time.time()  # Track the start time
    tracemalloc.start()  # Start memory tracking
    
    with open(output_file, 'w') as f:
        total_matches = 0
        
        with open(fasta_file, 'r') as genome_file:
            chromosome = None
            sequence = []
            
            # Process genome file and find matches
            for line in genome_file:
                line = line.strip()
                if line.startswith('>'):  # New chromosome header
                    if chromosome is not None:
                        total_matches += search_exact_matches(f, chromosome, ''.join(sequence), pattern)
                    chromosome = line[1:].split()[0]  # Extract chromosome name
                    sequence = []
                else:
                    sequence.append(line.upper())  # Append sequence lines
            
            if chromosome is not None:  # Final match processing
                total_matches += search_exact_matches(f, chromosome, ''.join(sequence), pattern)
    
        peak_memory, _ = tracemalloc.get_traced_memory()  # Get peak memory usage
        elapsed_time = time.time() - start_time  # Calculate elapsed time
        tracemalloc.stop()  # Stop memory tracking
        log_message(f, f"Exact matches found: {total_matches}")
        log_message(f, f"Time taken: {elapsed_time:.4f} seconds")
        log_message(f, f"Peak memory usage: {peak_memory / 1024:.2f} KB")

def search_approximate_matches(f, chromosome, genome_sequence, pattern):
    """
    Searches for approximate matches (with at most 1 mismatch) of the pattern in the genome sequence.
    
    Arguments:
        f (file object): The output file to log the matches.
        chromosome (str): The chromosome name where the search is performed.
        genome_sequence (str): The genome sequence to search within.
        pattern (str): The pattern to search for.
        
    Returns:
        int: The total number of approximate matches found.
    """
    total_matches = 0
    pattern_len = len(pattern)
    lps = compute_lps(pattern)  # Compute LPS array for pattern
    
    # Function to check if a substring is an approximate match (1 mismatch allowed)
    def is_approximate_match(text, pat, start):
        mismatches = 0
        mismatch_type = None
        i, j = start, 0
        while j < len(pat) and i < len(text):
            if text[i] == pat[j]:
                i += 1
                j += 1
            else:
                mismatches += 1
                if mismatches > 1:
                    return False, None  # More than 1 mismatch, not an approximate match
                # Check for possible substitutions, insertions, or deletions
                if i + 1 < len(text) and j + 1 < len(pat) and text[i + 1] == pat[j + 1]:
                    i += 1
                    j += 1  # Substitution
                    mismatch_type = "Substitution"
                elif i + 1 < len(text) and text[i + 1] == pat[j]:
                    i += 1  # Deletion 
                    mismatch_type = "Deletion" 
                elif j + 1 < len(pat) and text[i] == pat[j + 1]:
                    j += 1  # Insertion
                    mismatch_type = "Insertion"
                else:
                    return False, None
        if j < len(pat):
            mismatches += (len(pat) - j)
        return mismatches <= 1, mismatch_type  # Return true if 1 or fewer mismatches
    
    j = 0  # KMP state
    for i in range(len(genome_sequence)):
        while j > 0 and genome_sequence[i] != pattern[j]:
            j = lps[j - 1]
        if genome_sequence[i] == pattern[j]:
            j += 1
        else:
            j = 0
        start_idx = i - pattern_len + 1
        if j == pattern_len:  # Exact match
            match_type = "Exact match"
        else:
            is_one_mismatch, mismatch_type = is_approximate_match(genome_sequence, pattern, start_idx)
            if start_idx >= 0 and is_one_mismatch:
                match_type = f"Approximate match: {mismatch_type}"  
            else:
                continue  # Skip if neither exact nor approximate
        
        log_message(f, f"{chromosome}\t{start_idx}\t{start_idx + pattern_len - 1}\t{match_type}\n")
        total_matches += 1
        j = lps[j - 1] if j > 0 else 0
    
    return total_matches

def count_approximate_matches(fasta_file, pattern, output_file):
    """
    Counts the total number of approximate matches (1 mismatch allowed) of the pattern in the genome sequence.
    
    Arguments:
        fasta_file (str): The path to the genome FASTA file.
        pattern (str): The pattern to search for.
        output_file (str): The file to store the results.
    """
    start_time = time.time()  # Track the start time
    tracemalloc.start()  # Start memory tracking
    
    with open(output_file, 'w') as f:
        total_matches = 0
        
        with open(fasta_file, 'r') as genome_file:
            chromosome = None
            sequence = []
            
            # Process genome file and find matches
            for line in genome_file:
                line = line.strip()
                if line.startswith('>'):  # New chromosome header
                    if chromosome is not None:
                        total_matches += search_approximate_matches(f, chromosome, ''.join(sequence), pattern)
                    chromosome = line[1:].split()[0]  # Extract chromosome name
                    sequence = []
                else:
                    sequence.append(line.upper())  # Append sequence lines
            
            if chromosome is not None:  # Final match processing
                total_matches += search_approximate_matches(f, chromosome, ''.join(sequence), pattern)
    
        peak_memory, _ = tracemalloc.get_traced_memory()  # Get peak memory usage
        elapsed_time = time.time() - start_time  # Calculate elapsed time
        tracemalloc.stop()  # Stop memory tracking
        log_message(f, f"Total matches with at most 1 mismatch: {total_matches}")
        log_message(f, f"Time taken: {elapsed_time:.4f} seconds")
        log_message(f, f"Peak memory usage: {peak_memory / 1024:.2f} KB")
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()  # Command-line argument parser
    parser.add_argument("--genome", required=True, help="Path to genome FASTA file")  # Genome file path
    parser.add_argument("--pattern", required=True, help="Path to pattern FASTA file")  # Pattern file path
    parser.add_argument("--output", required=True, help="Output file for results")  # Output file path
    parser.add_argument("--approx", action='store_true', help="Use approximate matching")  # Use approximate matching if set
    args = parser.parse_args()  # Parse arguments
    
    # Open output file and start logging
    with open(args.output, 'w') as f:
        log_message(f, f"Loading genome sequence from {args.genome}")
        log_message(f, f"Reading pattern sequence from {args.pattern}")
        pattern_seq = read_fasta_pattern(args.pattern)  # Read the pattern sequence
        log_message(f, f"Pattern of length {len(pattern_seq)} loaded\n")
        
        gc.collect()  # Perform garbage collection
        # Count exact or approximate matches based on user's choice
        if args.approx:
            count_approximate_matches(args.genome, pattern_seq, args.output)
        else:
            count_exact_matches(args.genome, pattern_seq, args.output)
