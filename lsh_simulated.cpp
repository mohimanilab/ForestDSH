#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <vector>
#include <set>
#include <map>
#include <unordered_map>
#include <chrono>
#include <string>
#include <memory>
#include "common_funcs.hpp"
#include <string.h>

using namespace std;

typedef chrono::high_resolution_clock Clock;

// mass spectrometry related variables
string specfile = "PRM_plus.mgf";
string pepfile = "uniprot-proteome_homo_sapiens_reference_proteome_5_23_2016.fasta";
vector <int> ids{4403, 4164, 4162, 5157, 4163, 4910, 4911, 4404, 4165, 4405};    

// peptide names
vector <string> pep_names;

// The 2 sets of vectors that need to be compared
vector <vector <bool> > X;
vector <vector <bool> > Y;

// Length of a vector
int T = 10000;

// Number of vectors in each set
int N;

// The structure of all probabilities
probs P;

// Hashes and LSH parameters
int r, b;
vector <vector <int> > hashx;
vector <vector <int> > hashy;
vector <map <vector<bool>, pair<vector <int>, vector <int> > > > buckets;
vector <set<int> > match_hashes;

// fp and tp rates;
double fp, tp, expfp;


vector <vector<int> > lsh_hasher(vector<vector<bool> > vecs, int num_hashes) {
    
    vector <vector <int> > hashes;
    int num_vecs = vecs.size();
    int length_vec = vecs[0].size();
    int seed = 1234;

    for (int i = 0; i < num_hashes; i++) {
        int randnum = my_rand(&seed) % length_vec;
        
        vector <int> cur_hashes;
        for (int j = 0; j < num_vecs; j++) {
            cur_hashes.push_back(vecs[j][randnum]);
        }
        hashes.push_back(cur_hashes);
    }

    return hashes;
}


void get_expfp() {

    double s = P.p1X * P.p1Y + P.p0X * P.p0Y;
    expfp = 1 - pow(1 - pow(s, r), b);
    cout << "aa " << s << " " << expfp << endl;
}


void hashing_and_bucketing() {

    for (int i = 0; i < X.size(); i++) {
        int sum = 0;
        for (int j = 0; j < X[i].size(); j++) {
            if (X[i][j] == 1) sum++;
        }
        // cout<<sum<<endl;
    }

    buckets.resize(b);
    // cout<<"Hashing X\n";
    hashx = lsh_hasher(X, r * b);
    // cout<<"Hashing Y\n";
    hashy = lsh_hasher(Y, r * b); 
    
    // cout<<"Bucketing X\n";
    for (int i = 0; i < X.size(); i++) {
        vector <int> hashes;
        for (int j = 0; j < r * b; j++) hashes.push_back(hashx[j][i]);

        vector <bool> cur_vec;
        int bucket_num = 0;
        for (int j = 0; j < r * b; j++) {
            cur_vec.push_back(hashes[j]);

            if (j % r == r - 1) {
                
                if (buckets[bucket_num].find(cur_vec) == buckets[bucket_num].end()) {
                    vector <int> a{i};
                    vector <int> b;
                    buckets[bucket_num].insert(make_pair(cur_vec, make_pair(a, b)));
                }
                else buckets[bucket_num][cur_vec].first.push_back(i);

                cur_vec.clear();
                bucket_num++;
            }
        }
    }
    
    // cout<<"Bucketing Y\n";
    for (int i = 0; i < Y.size(); i++) {
        vector <int> hashes;
        for (int j = 0; j < r * b; j++) hashes.push_back(hashy[j][i]);
        
        vector <bool> cur_vec;
        int bucket_num = 0;
        for (int j = 0; j < r * b; j++) {
            cur_vec.push_back(hashes[j]);

            if (j % r == r - 1) {
                
                if (buckets[bucket_num].find(cur_vec) != buckets[bucket_num].end())
                    buckets[bucket_num][cur_vec].second.push_back(i);
                
                cur_vec.clear();
                bucket_num++;
            }
        }
    } 
}


void obtain_all_matches() {
     
    match_hashes.resize(N);
    long long total_matches = 0;
    long long real_matches = 0;
    
    for (int i = 0; i < b; i++) {
        for (auto it = buckets[i].begin(); it != buckets[i].end(); ++it) {            
            vector <int> vecsX = it->second.first;
            vector <int> vecsY = it->second.second;
            for (int j = 0; j < vecsX.size(); j++) {
                for (int j2 = 0; j2 < vecsY.size(); j2++) {
                    match_hashes[vecsX[j]].insert(vecsY[j2]);
                }
            }
        }
    }

    // find real and total matches
    for (int i = 0; i < N; i++) {
        total_matches += match_hashes[i].size();
        for (auto it = match_hashes[i].begin(); it != match_hashes[i].end(); ++it) {
            int dp = 0;
            int j = *it;
            for (int t = 0; t < T; t++) {
                if (X[i][t] && Y[j][t]) dp++;
            }
            if (dp > double(0.7 * P.p11 * T) && i == j) real_matches++;
        }
    }
    
    //cout<<total_matches<<" matches found\n";
    //cout<<real_matches<<" real matches found\n";
    
    fp = double(total_matches-real_matches) / (double (N) * double (N-1));
    tp = double(real_matches) / double(N);
}


void obtain_msgf_matches() {
    
    cout<<"Obtaining matches...\n";
    
    match_hashes.resize(X.size());
    long long total_matches = 0;
    long long real_matches = 0;
   
    for (int i = 0; i < b; i++) {
        for (auto it = buckets[i].begin(); it != buckets[i].end(); ++it) {            
            
            vector <bool> thehash = it->first;
            int sum = 0;
            for (int j = 0; j < thehash.size(); j++) sum += thehash[j];
            
            vector <int> vecsX = it->second.first;
            vector <int> vecsY = it->second.second;
            for (int j = 0; j < vecsX.size(); j++) {
                //cout<<vecsX[j]<<"\t";
                for (int j2 = 0; j2 < vecsY.size(); j2++) {
                    match_hashes[vecsX[j]].insert(vecsY[j2]);
                }
            }
            //cout<<vecsY.size()<<"\t";
            //cout<<sum<<endl;
        }
    }

    int max_spec_len = X[0].size()/2; 

    // find real and total matches
    for (int i = 0; i < X.size(); i++) {
        total_matches += match_hashes[i].size();
        cout<<total_matches<<endl;
        for (auto it = match_hashes[i].begin(); it != match_hashes[i].end(); ++it) {
            int j = *it;
            string pep = pep_names[j];
            double* t_prm = nullptr;
            double* t_srm = nullptr;
            double t_mass;
            //theoretical_vector(0, pep, t_prm, t_srm, t_mass);
            int dp = 0;
            for (int k = 0; k < pep.length() - 1; k++) dp += X[i][t_prm[k]];
            for (int k = 0; k < pep.length() - 1; k++) dp += X[i][t_srm[k] + max_spec_len];
        }
    }
}

vector<bool> str_to_bool(char *s) {
    vector<bool> v;
    for(int i = 0; i < strlen(s); i++) {
	if(s[i] == '0') {
	    v.push_back(false);
	}
	else {
	    v.push_back(true);
	}
    }
    return v;
}

int main(int argc, char** argv) {

    if (argc < 8) {
        cout<<"Too few arguments\n";
        return 0;
    }

    N = atoi(argv[1]);
    b = atoi(argv[2]);
    r = atoi(argv[3]);

    P = set_pq_values(atof(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]));
    get_expfp();

    // we dont need N and p00, p01, p10, p11 if we get mass spectrometry data, so those
    // parameters are just ignored in that case.

    //cout<<"Creating data\n";
    auto cstart = Clock::now(); 
    
    // The simulated data case
    
    pair <vector<vector<bool> >, vector<vector<bool> > > vecs;
    if(argc < 9) {
	vecs = make_data(P, N, T);
	X = vecs.first;
	Y = vecs.second;
    }
    else {
	cout << "Using sample data\n";
	char delim[] = ",";
	string d = argv[8];
	ifstream data(d);
	const int MAXDATA = 2 * 2000 + 2;
	char dline[MAXDATA];
	char *dtoken;
	while(data) {
	    data.getline(dline, MAXDATA);
	    if(dline[0] == '\0') {
		break;
	    }
	    dtoken = strtok(dline, delim);
	    vector<bool> new_x = str_to_bool(dtoken);
	    X.push_back(new_x);
	    dtoken = strtok(NULL, delim);
	    vector<bool> new_y = str_to_bool(dtoken);
	    Y.push_back(new_y);
	}
    }

    

    // The mass spectrometry case 
    
    /*pair<pair <vector<vector<bool> >, vector<vector<bool> > >, vector<string> > vecs;
    vecs = get_mass_spec_data(specfile, pepfile, ids);
    X = vecs.first.first;
    Y = vecs.first.second;
    pep_names = vecs.second;
    */

    auto cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n";
    
    //cout<<"Making hashes and putting them in buckets\n";
    cstart = Clock::now();
    hashing_and_bucketing();
    cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n"; 
    double t1 = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();

    //cout<<"Creating matches\n";
    cstart = Clock::now();
    obtain_all_matches();
    //obtain_msgf_matches();
    cend = Clock::now();
    //cout<<"Done in "<<chrono::duration_cast<chrono::milliseconds>(cend-cstart).count()<<" ms\n\n"; 
    double t2 = chrono::duration_cast<chrono::milliseconds>(cend-cstart).count();

    print_results(N, r, b, expfp, fp, tp, t1, t2);

    return 0;
}
