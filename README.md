# SuffixTreeLempelZiv

Gives the suffix tree and lempel-ziv factorization of a sequence

For help, use -h or --help

NAME
        STree_LZF
        Gives the suffix tree and Lempel-Ziv factorization of the input string
USAGE 1:
        STree_LZF -s input-sequence
USAGE 2:
        STree_LZF -f input-file [-b start-line] [-e end-line] [-t suffix-tree-output-file]
        [-l LZ-factorization-output-file [-p]]
NOTES:
        -b and -e indicate the line numbers within which the input file is processed.
        If -p is used, a list is given where the i'th number is the length  of longest
                prefix starting at i that also starts at j<i.
        If -s is used, all outputs are printed on screen.
