# SuffixTreeLempelZiv
Farzad Farnoud, California Institute of Technology, September 2014

Uses code from Mark Nelson

## Synopsis
This program gives the suffix tree and lempel-ziv factorization of a sequence.
In adition it can give the length of the longest prefix for each position that has appeard before.

## Code Example
For help, use STree_LZF -h or STree_LZF --help

        Gives the suffix tree and Lempel-Ziv factorization of the input string

USAGE 1: SCREEN MODE

        STree_LZF -s input-sequence

USAGE 2: FILE MODE

        STree_LZF -f input-file [-b start-line] [-e end-line] [-t suffix-tree-output-file]
        [-l LZ-factorization-output-file [-p]]

NOTES for SCREEN MODE:

        The input is given as a sequence, and all outputs are printed on screen.
	Example:
	        >STree_LZF -s banana
		Checking uniqueness of last character...Last character not unique. Adding '$'.
		
		Initializing suffix tree...done!
		
		Constructing suffix tree of sequence of length 7...done!
		Suffix tree has 11 nodes.
		
		Computing LZF...done!
		
		Suffix tree:
		Start  End  Suf  First Last  String
		0     1  -1    -1      6  banana$
		4     2  -1     3      6  na$
		6     3  -1     3      6  na$
		8     4   6     1      3  na
		4     5  -1     5      6  $
		0     6   8     1      3  na
		6     7  -1     5      6  $
		0     8   0     0      1  a
		8     9  -1     5      6  $
		0    10  -1     5      6  $
		done!
			
		Block lengths in Lempel-Ziv factorization:
		1 1 1 3 1
		
		The i'th number is the length of longest prefix starting at i that also starts at j<i:
		0 0 0 3 2 1 0
		done!
		
		Press Enter to exit...
		
NOTES for FILE MODE:

        The input is given as a file. The input sequence is the sequence starting on 
		line start-line and ending on end-line. Newline characters are removed from the sequence.
	If -t is given, the suffix tree is printed in file suffix-tree-output-file.
	If -l is given block length of the Lempel-Ziv factorization are given in file 
		LZ-factorization-output-file.
	If -p is used, a list is appended to LZ-factorization-output-file, where the 
		i'th number is the length  of longest prefix starting at i that also starts at j<i.
