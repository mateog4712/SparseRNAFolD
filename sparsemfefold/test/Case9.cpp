#include <string>
#include <fstream>
#include <vector>
#include <iostream>
#include <tuple>

#include <iomanip>

#include <sstream>

#include <limits>
#include <iterator>
#include <cstring>
#include <cassert>
#include <numeric>
#include <map>
//#define CATCH_CONFIG_MAIN
#include "catch.hpp"

// Make random input structures for all of the sequences given in the file
using namespace std;

string make_structure(string seq);
bool check_Pseudoknot(vector<tuple<int,int> > used, int i, int j);

static int pairs[8][8] =
  /* _  A  C  G  U  X  K  I */
{ { 0, 0, 0, 0, 0, 0, 0, 0 },
  { 0, 0, 0, 0, 5, 0, 0, 5 },
  { 0, 0, 0, 1, 0, 0, 0, 0 },
  { 0, 0, 2, 0, 3, 0, 0, 0 },
  { 0, 6, 0, 4, 0, 0, 0, 6 },
  { 0, 0, 0, 0, 0, 0, 2, 0 },
  { 0, 0, 0, 0, 0, 1, 0, 0 },
  { 0, 6, 0, 0, 5, 0, 0, 0 } };

int main(int argc,char **argv) {
    vector<string> names;
    vector<string> seqs;
    string family;

    string file = "/home/mgray7/output2/fasta/";
    cout << "Put Family: ";
    cin >> family;
    file = file + family + ".txt";
    ifstream in(file);

    string str;
    int i = 0;
    while(getline(in,str)){
    // cout << str << endl;
        if(i%2 == 0){
            names.push_back(str);
        }
        else if(i%2==1){
            seqs.push_back(str);
        }
        ++i;
    }  
    in.close();

    // string intermediary = make_structure(seqs[0]);
    // cout << intermediary << endl;
    // exit(0);

    string command = "./build/SparseMFEFold -d1";
    // string command2 = "/home/mgray7/includes/vienna/bin/RNAfold -d1 --noPS ";
    string command2 = "/home/mgray7/RNAfold/build/src/RNAFold -d1 --noPS ";
    double score = 0;
    int size = seqs.size();
    // exit(0);
    string filename = "/home/mgray7/output2/structures/" + family +".txt";
    ofstream out3(filename);
    std::vector<int> wrong;
    for(int i =0;i<size;++i){
        string intermediary = make_structure(seqs[i]);
        
        out3 << names[i] << endl;
        out3 << intermediary << endl;
    
    }
    out3.close();
    return 0;
}

string make_structure(string seq){
    std::map<char,uint> base;
    base['A']=1;
    base['C']=2;
    base['G']=3;
    base['U']=4;

    int len = seq.length();
    string structure (len,'.');
    int k = 0;
    vector<tuple<int,int> > used;
    int count = 0;
    while(k<.5*log2(len)){
        
        std::random_device dev;
        std::mt19937 rng(dev());
        std::uniform_int_distribution<std::mt19937::result_type> dist6(0,len); // distribution in range [0, len]
        int i = dist6(rng);
        std::uniform_int_distribution<std::mt19937::result_type> dist7(i,len); // distribution in range [0, len]
        int j = dist7(rng);
        
        bool knot = check_Pseudoknot(used,i,j);

        if(j-i>3 && !knot){
            if(pairs[base[seq[i]]][base[seq[j]]] > 0){
                used.push_back(make_tuple(i,j));
                structure[i] = '(';
                structure[j] = ')';

                ++k;
            }
            
            
        }
        count++;
        if(count>250){
            count = 0;
            k = 0;
            used.clear();
            structure = string(len,'.');
        }
    }

    return structure;

}

bool check_Pseudoknot(vector<tuple<int,int> > used, int i, int j){
  for(int m = 0; m<used.size();++m){
        int k = get<0>(used[m]);
        int l = get<1>(used[m]);
        // cout << i << "\t" << k << "\t" << l << "\t" << j << endl;
        if(i==k || j==l || i==l || j==k) return true;
        if((i < k  && j > k && j < l) || (i < l  && j > l && i > k)) return true;
    }
  return false;
}