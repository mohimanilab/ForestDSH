# ForestDSH

# Compiling

For compiling on the cluster, use the following

/usr/bin/g++ ForestDSH.cpp -o ForestDSH -std=c++0x -O3

/usr/bin/g++ lambda.cpp -o lambda -std=c++0x -O3

/usr/bin/g++ lsh\_simulated.cpp common\_funcs.cpp -o lsh_simulated -std=c++0x -O3

/usr/bin/g++ minhash\_simulated.cpp common\_funcs.cpp -o minhash_simulated -std=c++0x -O3

# Running

ForestDSH takes in 8 arguments with additional options for a data file, Q matrix file, and a flag for R/CR based data generation.

`usage: ./ForestDSH <N> <S> <C1> <C2> <C3> <matrix> <b> [<data>] [-qmatrix=<qmatrix>][-test_r_cr]`

N is the number of pairs to generate

S is the length of each string

C1, C2, and C3 are the constant multipliers that affect whether or not a bucket is accepted or killed

matrix should be a csv file that represents the probability matrix. NOTE: probability matrices with more than 10 rows or columns are currently unsupported.

b is the number of bands to use. If it is -1, then the program will use the number of bands necessary to reach 99% true positive rate.

data is an optional argument that should be a csv file that represents the strings. If left blank, N pairs of S length strings will be randomly generated based on the matrix. Otherwise, the format should look like the following:

X1,Y1

X2,Y2

...

XN, YN

* Furthermore, S should be the length of the longest string. 

qmatrix is an optional argument that should be a csv file representing the probability matrix Q for pairs that are independent of each other. The utility in the ProbEst folder can be used to generate these matricies, as well as the P matrix, based on L1 norms R and CR. See folder and documentation for details.

test_r_cr is a flag that causes data to be generated according to the R / CR L1 norm. Currently, we always create pairs with L1 norm less than R. Note, currently the R value is hardcoded in the cpp file.


Example usage: ./ForestDSH 2000 10000 0.5 0.03125 0.03125 matrix_csv/matrix_1.csv -1

Example usage with data file: ./ForestDSH 2000 10000 0.5 0.03125 0.03125 matrix_csv/matrix_1.csv -1 test/sample_data.txt

Output

ForestDSH outputs two lines in csv format as well as the total time it took to run (greater than the time the algorithm took since it includes the time needed to generate the data).

The headers are fairly self-explanatory.

Example output:

C1,C2,C3,Total complexity,Complexity of hashing,Complexity of tree construction,Complexity of FP,#bands needed,Overall empirical TP,TP at one band,Theoretical TP in one band,Overall empirical FP, empirical FP in one band,Theoretical FP in one band,#buckets,#nodes in the tree,total number of x nodes,total number of y nodes,gamma_x,gamma_y,mapped_x,mapped_y

0.5,0.03125,0.03125,2235.54,1340,208,687.54,19,8335,0.219399,0.218657,37501,0.000493434,0.000495448,58936,241468,11499,65402,10.2911,2.21452,391686,83928

Total run time: 5100


Example usage (With R/CR data generation, user defined Q matrix):

./ForestDSH 2000 1000 0.5 0.03125 0.03125 ProbEst/matrix_S1000_R100_L2.csv -1 -qmatrix=ProbEst/matrix_S1000_R400_L2.csv -test_r_cr

minhash\_simulated takes in the following as arguments:

`./minhash_simulated <N> <b> <r> <p1> <p2> <p3> <p4> [<data>]`

N is the number of data points

b and r are the tuned hyperparameters

p1, p2, p3, and p4 are the values from the probability matrix.

data is an optional argument that should be a csv file that represents the strings. If left blank, N pairs of S length strings will be randomly generated based on the matrix. The format of the csv should be the same as the one described in the ForestDSH section.

Example usage

./minhash\_simulated 2000 800 8 0.345  0 0.31 0.345 sample_data.txt

The output is in the following format:

[N] [r] [b] [expected fp rate] [observed fp rate] [observed tp rate] [time for hashing] [time for bucketing and finding matches] [total time] [GOOD TP if observed tp rate >= 0.99, blank otherwise]

lsh\_simulated can be ran in the exact same manner as minhash\_simulated.

lhfast and mhfast are the same as lsh\_simulated and minhash\_simulated respectively, but they use an estimate of 0.015ms for each positive check instead of performing checks on each positive. This speeds up the actual runtime at the cost of overestimating the time required.
