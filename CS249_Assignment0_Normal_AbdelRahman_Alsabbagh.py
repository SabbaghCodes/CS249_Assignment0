import argparse
import gc
import time
import tracemalloc

def read_fasta_sequence(file_path):
    """
    Reads a genome FASTA file and returns a dictionary where the keys are chromosome names 
    and the values are the corresponding sequences.
    
    Arguments:
        file_path (str): The path to the FASTA file containing genome sequences.
    
    Returns:
        dict: A dictionary with chromosome names as keys and genome sequences as values.
    """
    genome = {}
    current_chromosome = None
    current_sequence = []
    
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):  
                if current_chromosome is not None:
                    genome[current_chromosome] = ''.join(current_sequence)
                current_chromosome = line[1:].split()[0]  # extracting chromosome name
                current_sequence = []
            else:
                current_sequence.append(line.upper())
    
    if current_chromosome is not None:
        genome[current_chromosome] = ''.join(current_sequence) 
    
    return genome

def read_fasta_pattern(file_path):
    """
    Reads a pattern FASTA file and returns the pattern as a single string.
    
    Arguments:
        file_path (str): The path to the FASTA file containing the pattern sequence.
    
    Returns:
        str: The pattern sequence as a single concatenated string.
    """
    pattern = []
    with open(file_path, 'r') as f:
        for line in f:
            line = line.strip()
            if not line.startswith('>'):
                pattern.append(line.upper())
    return ''.join(pattern)

def log_message(f, message):
    """
    Logs a message by printing it to the console and writing it to a file.
    
    Arguments:
        f (file object): The file object to write the message to.
        message (str): The message to log.
    """
    print(message)
    f.write(message + "\n")

def compute_lps(pattern):
    """
    Computes the Longest Prefix Suffix (LPS) array for the given pattern.
    This array is used in the Knuth-Morris-Pratt (KMP) string matching algorithm.
    
    Arguments:
        pattern (str): The pattern string to compute the LPS array for.
    
    Returns:
        list: A list representing the LPS array for the given pattern.
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

def count_exact_matches(genomes, pattern, output_file):
    """
    Counts the exact matches of the pattern in the genome sequences using the KMP algorithm.
    Writes the results (chromosome, start, end, match type) to the output file.
    
    Arguments:
        genomes (dict): A dictionary of genome sequences with chromosome names as keys.
        pattern (str): The pattern string to search for in the genome sequences.
        output_file (str): The path to the output file to write the results.
    """
    start_time = time.time()
    tracemalloc.start()
    
    with open(output_file, 'w') as f:
        total_matches = 0
        for chromosome, genome_sequence in genomes.items():
            lps = compute_lps(pattern)
            pattern_len = len(pattern)
            j = 0  # KMP state
            
            for i in range(len(genome_sequence)):
                while j > 0 and genome_sequence[i] != pattern[j]:
                    j = lps[j - 1]
                if genome_sequence[i] == pattern[j]:
                    j += 1
                else:
                    j = 0
                if j == pattern_len:
                    start = i - pattern_len + 1
                    end = start + pattern_len - 1
                    f.write(f"{chromosome}\t{start}\t{end}\tExact match\n")
                    total_matches += 1
                    j = lps[j - 1]
    
    peak_memory, _ = tracemalloc.get_traced_memory()
    elapsed_time = time.time() - start_time
    tracemalloc.stop()
    log_message(f, f"Exact matches found: {total_matches}") 
    log_message(f, f"Time taken: {elapsed_time:.4f} seconds")
    log_message(f, f"Peak memory usage: {peak_memory / 1024:.2f} KB")

def count_approximate_matches(genomes, pattern, output_file):
    """
    Counts the approximate matches (with at most 1 mismatch) of the pattern in the genome sequences.
    Writes the results (chromosome, start, end, match type) to the output file.
    
    Arguments:
        genomes (dict): A dictionary of genome sequences with chromosome names as keys.
        pattern (str): The pattern string to search for in the genome sequences.
        output_file (str): The path to the output file to write the results.
    """
    start_time = time.time()
    tracemalloc.start()
    
    with open(output_file, 'w') as f:
        total_matches = 0
        pattern_len = len(pattern)
        lps = compute_lps(pattern)
        
        def is_approximate_match(text, pat, start):
            """
            Checks if a match of the pattern in the genome sequence is approximate, 
            allowing at most 1 mismatch.
            
            Arguments:
                text (str): The genome sequence to check the pattern against.
                pat (str): The pattern string to check.
                start (int): The starting index of the match.
            
            Returns:
                bool: True if the match is approximate, False otherwise.
            """
            mismatches = 0
            match_type = ""
            i, j = start, 0
            while j < len(pat) and i < len(text):
                if text[i] == pat[j]:
                    i += 1
                    j += 1
                else:
                    mismatches += 1
                    if mismatches > 1:
                        return False, ""
                    if i + 1 < len(text) and j + 1 < len(pat) and text[i + 1] == pat[j + 1]:
                        i += 1
                        j += 1  # substitution
                        match_type = "Substitution"
                    elif i + 1 < len(text) and text[i + 1] == pat[j]:
                        i += 1  # insertion
                        match_type = "Insertion"
                    elif j + 1 < len(pat) and text[i] == pat[j + 1]:
                        j += 1  # deletion 
                        match_type = "Deletion"
                    else:
                        return False, ""
            if j < len(pat):
                mismatches += (len(pat) - j)
            return (mismatches <= 1), match_type
        
        for chromosome, genome_sequence in genomes.items():
            j = 0  # KMP state
            for i in range(len(genome_sequence)):
                while j > 0 and genome_sequence[i] != pattern[j]:
                    j = lps[j - 1]
                if genome_sequence[i] == pattern[j]:
                    j += 1
                else:
                    j = 0
                start_idx = i - pattern_len + 1
                if j == pattern_len:
                    match_type = "Exact match"
                else:
                    approximate, match_type = is_approximate_match(genome_sequence, pattern, start_idx)
                    if not approximate:
                        continue
                
                f.write(f"{chromosome}\t{start_idx}\t{start_idx + pattern_len - 1}\t{match_type}\n")
                total_matches += 1
                j = lps[j - 1] if j > 0 else 0
    
    peak_memory, _ = tracemalloc.get_traced_memory()
    elapsed_time = time.time() - start_time
    tracemalloc.stop()
    log_message(f, f"Total matches with at most 1 mismatch: {total_matches}")
    log_message(f, f"Time taken: {elapsed_time:.4f} seconds")
    log_message(f, f"Peak memory usage: {peak_memory / 1024:.2f} KB")

if __name__ == "__main__":
    """
    This script matches a given pattern against a genome sequence using either exact matching 
    or approximate matching (with at most 1 mismatch). The results are written to the specified output file.
    
    Arguments:
        --genome: Path to the genome FASTA file
        --pattern: Path to the pattern FASTA file
        --output: Path to the output file where results will be saved
        --approx: (Optional) Flag to perform approximate matching instead of exact matching
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--genome", required=True, help="Path to genome FASTA file")
    parser.add_argument("--pattern", required=True, help="Path to pattern FASTA file")
    parser.add_argument("--output", required=True, help="Output file for results")
    parser.add_argument("--approx", action='store_true', help="Use approximate matching")
    
    args = parser.parse_args()
    with open(args.output, 'w') as f:
        log_message(f, f"Reading genome sequence from {args.genome}")
        genome_seq = read_fasta_sequence(args.genome)
        log_message(f, f"{len(genome_seq)} chromosomes loaded\n")
        
        log_message(f, f"Reading pattern sequence from {args.pattern}")
        pattern_seq = read_fasta_pattern(args.pattern)
        log_message(f, f"Pattern of length {len(pattern_seq)} loaded\n")
        
        gc.collect()
        if args.approx:
            count_approximate_matches(genome_seq, pattern_seq, args.output)
        else:
            count_exact_matches(genome_seq, pattern_seq, args.output)
