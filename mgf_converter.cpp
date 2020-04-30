//
//  main.cpp
//  DSHTree
//
//  Created by Sean Chang on 4/25/19.
//  Copyright © 2019 Sean Chang. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cmath>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <chrono>
#include <stdio.h>
#include <string.h>

typedef std::chrono::high_resolution_clock Clock;

//Keeping track of the number of accepts, rejects, alpha, and beta
int accepts = 0;
int rejects = 0;
double TFP = 0.0;
double TTP = 0.0;
double gx = 0.0;
double gy = 0.0;
int totalX = 0;
int totalY = 0;

//Debugging purposes
int numNodes = 0;
int xnodes = 0;
int ynodes = 0;
int deepest = 0;

//Toggle this boi
double C1;

double p1;
double p2;
double p3;
double p4;

std::vector<bool> found;
int unique;
std::vector<std::vector<bool>> ufp;
int funique;

//2D vector that holds doubles
typedef std::vector<std::vector<double>> arr;
typedef std::vector<std::vector<long double>> larr;

//p and q are 2D vectors that represent a probability matrix
struct probs_header {
    larr p, q;
    std::vector<long double> px, py, qx, qy;
};
typedef struct probs_header probs;

//Special string structure so we can tell which string is paired with which string
struct istring {
    std::string str;
    int id;
};

//Represents a bucket
//xstr and ystr are the strings at that bucket.
//xmatches and ymatches represent the strings that go into the bucket.
struct bucket {
    std::string xstr, ystr;
    std::vector<istring*> xmatches, ymatches;
};

//Tree node
//Contains a vector of its children, what string it represents, and all the buckets at that node
struct tnode {
    std::vector<tnode*> children;
    std::string str;
    std::vector<bucket*> bvec;
};

//Copied from common_funcs.cpp
//Xcode linker why don't you work
#define MODULO 100000

probs* set_pq_values(larr pmatrix) {
    
    probs *ret = new probs();
    larr qmatrix;
    
    ret -> p = pmatrix;
    
    std::vector<long double> rowsums(pmatrix.size(), 0.0);
    std::vector<long double> colsums(pmatrix.size(), 0.0);
    
    for(int i = 0; i < pmatrix.size(); i++) {
        for(int j = 0; j < pmatrix.size(); j++) {
            rowsums.at(i) += pmatrix.at(i).at(j);
            colsums.at(j) += pmatrix.at(i).at(j);
        }
    }
    
    ret -> px = rowsums;
    ret -> py = colsums;
    
    for(int i = 0; i < pmatrix.size(); i++) {
        std::vector<long double> newrow;
        for(int j = 0; j < pmatrix.size(); j++) {
            newrow.push_back(rowsums.at(i) * colsums.at(j));
        }
        qmatrix.push_back(newrow);
    }
    
    ret -> q = qmatrix;
    
    std::vector<long double> qx(qmatrix.size(), 0.0);
    std::vector<long double> qy(qmatrix.size(), 0.0);
    
    for(int i = 0; i < qmatrix.size(); i++) {
        for(int j = 0; j < qmatrix.size(); j++) {
            qx.at(i) += qmatrix.at(i).at(j);
            qy.at(j) += qmatrix.at(i).at(j);
        }
    }
    
    ret -> qx = qx;
    ret -> qy = qy;
    
    return ret;
}

//The probability that x and y were generated given p
double matchProbability(probs *prob, std::string x, std::string y) {
    double product = 1.0;
    for(int i = 0; i < x.length(); i++) {
        int xind = std::stoi(x.substr(i, 1));
        int yind = std::stoi(y.substr(i, 1));
        product *= (prob -> p).at(xind).at(yind);
    }
    return product;
}

double matchProbabilityQ(probs *prob, std::string x, std::string y) {
    double product = 1.0;
    for(int i = 0; i < x.length(); i++) {
        int xind = std::stoi(x.substr(i, 1));
        int yind = std::stoi(y.substr(i, 1));
        product *= (prob -> q).at(xind).at(yind);
    }
    return product;
}

std::string ionToString(std::vector<std::pair<int, std::string>> strInfo, int S) {
    std::string str = "";
    for(int i = 0; i < S; i++) {
        str += "0";
    }
    for(int j = 0; j < strInfo.size(); j++) {
        str.replace(strInfo.at(j).first, 1, strInfo.at(j).second);
    }
    return str;
}

std::string convertString(std::string s) {
    for(int i = 0; i < s.length(); i++) {
        if(s.substr(i, 1).compare("0") == 0) {
            s.replace(i, 1, "1");
        }
        else if(s.substr(i, 1).compare("2") == 0 || s.substr(i, 1).compare("3") == 0) {
            s.replace(i, 1, "0");
        }
    }
    return s;
}

void convertData(probs *prob, std::string file, std::string out, double threshold, int S) {
    std::string line;
    std::ifstream data (file);
    int id = 0;
    bool isX = true;
    std::vector<istring*> xstrings, ystrings;
    if(data.is_open()) {
        while(getline (data, line) ) {
            if(line == "BEGIN IONS") {
                std::vector<std::pair<int, std::string>> strInfo;
                while(getline (data, line)) {
                    if(line == "END IONS") {
                        std::string entry = ionToString(strInfo, S);
                        if(isX) {
                            istring *x = new istring();
                            x -> str = entry;
                            x -> id = id;
                            xstrings.push_back(x);
                            isX = false;
                        }
                        else {
                            double pval = matchProbability(prob, xstrings.at(id) -> str, entry);
                            double qval = matchProbabilityQ(prob, xstrings.at(id) -> str, entry);
                            double val = pval / qval;
                            if(isnan(val) || val < threshold) {
                                xstrings.pop_back();
                            }
                            else {
                                istring *y = new istring();
                                y -> str = entry;
                                y -> id = id;
                                ystrings.push_back(y);
                                id++;
                            }
                            isX = true;
                        }
                        break;
                    }
                    int index = (int)(stod(line.substr(0, line.find("\t"))) + 0.5);
                    std::string value = (line.substr(line.length() - 1, 1));
                    strInfo.push_back(std::make_pair(index, value));
                }
                
            }
        }
        data.close();
        std::cout << "Done reading data\n";
        std::cout << "Number of x-strings: " << xstrings.size() << "\n";
        std::cout << "Number of y-strings: " << ystrings.size() << "\n";
        
        std::ofstream outfile;
        outfile.open(out);
        for(int i = 0; i < xstrings.size(); i++) {
            std::string xstr = xstrings.at(i) -> str;
            std::string ystr = ystrings.at(i) -> str;
            
            outfile << xstr << "," << ystr << "\n";
        }
        outfile.close();
    }
}

larr read_matrix(std::string matrix_filename){
    std::ifstream file(matrix_filename);
    larr p;
    const int MAXLINE = 256;
    char line[MAXLINE];
    char delim[] = ",";
    char *token;
    while(file) {
        std::vector<long double> newRow;
        file.getline(line, MAXLINE);
        if(line[0] == '\0') {
            break;
        }
        token = strtok(line, delim);
        newRow.push_back(atof(token));
        while((token = strtok(NULL, delim)) != NULL) {
            newRow.push_back(atof(token));
        }
        p.push_back(newRow);
    }
    return p;
}


int main(int argc, const char * argv[]) {
    if(argc < 5) {
        std::cout << "usage: ./mgf_converter <matrix> <infile> <outfile> <threshold> <dimensions>\n";   
    }
    larr pmatrix = read_matrix(argv[1]);
     
    probs *prob = set_pq_values(pmatrix);
    
    convertData(prob, argv[2], argv[3], atof(argv[4]), atoi(argv[5]));
    
    
    return 0;
}
