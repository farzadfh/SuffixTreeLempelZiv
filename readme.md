# SuffixTreeLempelZiv

## Synopsis
This program gives the suffix tree and Lempel-Ziv factorization of a sequence. In adition it can give the length of the longest prefix for each position that has appeard before. The main motivation is to get Lempel-Ziv factorization in linear time.

## Code Example
For help, use `STree_LZF -h` or `STree_LZF --help`.

The code can be used in two ways. First, when the input sequence is given on screen. In this case the ouput is also printed on screen:
```
> STree_LZF -s banana

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
``` 
Note that if the final character is not unique $ is added. Based on the output, the blocks of the Lempel-Ziv factorization of `banana$` are `|b|a|n|ana|$|` with lengths 1, 1, 1, 3, 1.

In the second usage, the input is a file and the output is also printed in files. 


_Farzad Farnoud, http://farnoud.info_  
_California Institute of Technology, September 2014_  
_The output routine uses code from Mark Nelson_  
