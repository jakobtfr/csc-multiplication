### Sparse matrix multiplication of two matrices in Compressed Sparse Column (CSC) - Format
How to use: \
Input: Two matrices in .txt files in the following format:
  * rows, cols
  * values
  * row_indices
  * col_ptr

Example:
  * 4,4
  * 5,0.5,6,1,3
  * 0,2,1,1,3
  * 0,2,3,4,5

Output in the same format.

#### CLI commands:
entirely optional, no commands will use a standard value for the execution.\
Some commands:\
 -V <Number>          Specify the version of the function.\
 -a <Filename>        Specify file containing Matrix a.\
 -b <Filename>        Specify file containing Matrix b.\
 -o <Filename>        Specify the output file .\
use -h to get a detailed overview
