import unittest

def compute_lps(pattern):
    """
    Computes the Longest Prefix Suffix (LPS) array for the given pattern.
    
    The LPS array is used in the Knuth-Morris-Pratt (KMP) string matching algorithm
    to avoid redundant comparisons. It stores the length of the longest proper prefix
    which is also a suffix for each prefix of the pattern.
    
    Args:
        pattern (str): The pattern to compute the LPS array for.
        
    Returns:
        list: The LPS array where each element represents the length of the longest 
              proper prefix which is also a suffix for the corresponding prefix of the pattern.
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

def search_exact_matches(text, pattern):
    """
    Searches for exact matches of the given pattern in the text using the KMP algorithm.
    
    The function uses the LPS array to efficiently search for all occurrences of the 
    pattern in the text.
    
    Args:
        text (str): The text in which to search for the pattern.
        pattern (str): The pattern to search for in the text.
        
    Returns:
        list: A list of starting indices of exact matches of the pattern in the text.
    """
    total_matches = 0
    lps = compute_lps(pattern)
    pattern_len = len(pattern)
    j = 0  # KMP state
    matches = []
    
    for i in range(len(text)):
        while j > 0 and text[i] != pattern[j]:
            j = lps[j - 1]
        if text[i] == pattern[j]:
            j += 1
        else:
            j = 0
        if j == pattern_len:
            start = i - pattern_len + 1
            matches.append(start)
            total_matches += 1
            j = lps[j - 1]
    
    return matches

def search_approximate_matches(text, pattern):
    """
    Searches for approximate matches (with at most 1 mismatch) of the given pattern in the text.
    
    The function attempts to find matches by considering insertions, deletions, and substitutions
    as part of the mismatch handling, and it returns the starting indices of the matches along with 
    the type of mismatch ('I' for insertion, 'D' for deletion, 'S' for substitution, 'E' for exact).
    
    Args:
        text (str): The text in which to search for the pattern.
        pattern (str): The pattern to search for in the text.
        
    Returns:
        list: A list of tuples where each tuple contains the starting index of the match and the type of match:
              'E' for exact match, 'I' for insertion, 'D' for deletion, 'S' for substitution.
    """
    total_matches = 0
    pattern_len = len(pattern)
    matches = []
    
    def is_approximate_match(text, pat, start):
        """
        Checks if a mismatch at the given starting position is within the allowed limit 
        (at most 1 mismatch) and returns the mismatch type if applicable.
        
        Args:
            text (str): The text in which to search for the pattern.
            pat (str): The pattern to search for.
            start (int): The starting index in the text to check for an approximate match.
        
        Returns:
            tuple: A tuple containing a boolean indicating whether the match is valid and the mismatch type 
                   ('S', 'D', 'I', or None if no mismatch is found).
        """
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
                    return False, None
                if i + 1 < len(text) and j + 1 < len(pat) and text[i + 1] == pat[j + 1]:
                    i += 1
                    j += 1  # substitution
                    mismatch_type = "S"
                elif i + 1 < len(text) and text[i + 1] == pat[j]:
                    i += 1  # deletion 
                    mismatch_type = "D"
                elif j + 1 < len(pat) and text[i] == pat[j + 1]:
                    j += 1  # insertion
                    mismatch_type = "I"
                else:
                    return False, None
        if j < len(pat):
            mismatches += (len(pat) - j)
        return mismatches <= 1, mismatch_type
    
    for i in range(len(text) - pattern_len + 1):
        if text[i:i + pattern_len] == pattern:
            matches.append((i, "E"))
            total_matches += 1
        else:
            match_found, mismatch_type = is_approximate_match(text, pattern, i)
            if match_found:
                matches.append((i, mismatch_type))
                total_matches += 1
    
    return matches

class TestPatternMatching(unittest.TestCase):
    """
    Unit test cases for the pattern matching functions (exact and approximate).
    """
    
    def test_exact_match(self):
        """
        Test case for exact matches of the pattern in the text.
        """
        text = "ACGTACGTACGT"
        pattern = "ACGT"
        self.assertEqual(search_exact_matches(text, pattern), [0, 4, 8])
    
    def test_exact_match_single_occurrence(self):
        """
        Test case for a single occurrence of the exact match.
        """
        text = "TTTACGTAAA"
        pattern = "ACGT"
        self.assertEqual(search_exact_matches(text, pattern), [3])
    
    def test_exact_match_no_occurrence(self):
        """
        Test case for no exact matches of the pattern in the text.
        """
        text = "TTTTTT"
        pattern = "ACGT"
        self.assertEqual(search_exact_matches(text, pattern), [])
    
    def test_approximate_match_1(self):
        """
        Test case for approximate matches with different mismatch types.
        """
        text = "ACGTACCTACGT"
        pattern = "ACGT"
        self.assertEqual(search_approximate_matches(text, pattern), [(0, 'E'), (1, 'I'), (4, 'S'), (7, 'D'), (8, 'E')])
    
    def test_approximate_match_2(self):
        """
        Test case for approximate matches with a different set of mismatches.
        """
        text = "ACGTAGCTACGT"
        pattern = "ACGT"
        self.assertEqual(search_approximate_matches(text, pattern), [(0, 'E'), (1, 'I'), (7, 'D'), (8, 'E')])
    
    def test_approximate_match_3(self):
        """
        Test case for approximate matches with multiple mismatches in the text.
        """
        text = "ACGTTACGTACGT"
        pattern = "ACGT"
        self.assertEqual(search_approximate_matches(text, pattern), [(0, 'E'), (1, 'I'), (4, 'D'), (5, 'E'), (6, 'I'), (8, 'D'), (9, 'E')])
    
    def test_approximate_match_4(self):
        """
        Test case for approximate matches with mismatches at different positions.
        """
        text = "ACGACGTACGT"
        pattern = "ACGT"
        self.assertEqual(search_approximate_matches(text, pattern), [(2, 'D'), (3, 'E'), (4, 'I'), (6, 'D'), (7, 'E')])
    
    def test_approximate_match_no_occurrence(self):
        """
        Test case for no approximate matches of the pattern in the text.
        """
        text = "TTTTTT"
        pattern = "ACGT"
        self.assertEqual(search_approximate_matches(text, pattern), [])

if __name__ == "__main__":
    unittest.main()
