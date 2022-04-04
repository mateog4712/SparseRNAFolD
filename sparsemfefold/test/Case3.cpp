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

    string command = "./build/src/SparseMFEFold ";
    string command2 = "/home/mgray7/RNAfold/build/src/RNAFold ";
    double score = 0;
    int size = seqs.size();
    // exit(0);
    string filename = "/home/mgray7/output2/results/" + family +".txt";
    ofstream out3(filename);
    for(int i =0;i<size;++i){
        cout << i << endl;
        string intermediary = make_structure(seqs[i]);
        cout << intermediary << endl;
        
        string commands = command + " -r \"" + intermediary + "\" " + seqs[i] +  " > out.txt";
        system(commands.c_str());

        ifstream in1("out.txt");
        getline(in1,str);
        getline(in1,str);
        in1.close();
        // exit(0);
        string structure = str.substr(0,seqs[i].length());
        

        double energy = stod(str.substr(seqs[i].length()+2,str.length()-seqs[i].length()-1));
        // exit(0);
        ofstream out1("/home/mgray7/RNAfold/in.txt");
        out1 << names[i] << endl;
        out1 << seqs[i] << endl;
        out1 << intermediary << endl;
        out1.close();
        string infile = "/home/mgray7/RNAfold/in.txt";
    
        string commands2 = command2 + "--enforceConstraint --constraint<" + infile + " > /home/mgray7/RNAfold/out.txt";
        system(commands2.c_str());

        ifstream in2("/home/mgray7/RNAfold/out.txt");
        getline(in2,str);
        getline(in2,str);
        getline(in2,str);
        in2.close();
        // exit(0);
        // exit(0);
        // cout << str.length() << " " << seqs[i].length()+2 << " " << str.length()-seqs[i].length()-1;
        string structure2 = str.substr(0,seqs[i].length());
        // exit(0);
        
        double energy2 = stod(str.substr(seqs[i].length()+2,str.length()-seqs[i].length()-1));
        // exit(0);
        out3 << names[i] << endl;
        out3 << seqs[i] << endl;
        out3 << intermediary << endl;
        out3 << structure << endl;
        out3 << structure2 << endl;
        // for(int k=0;k<structure.length();++k){
        //     if(structure[k] == '[') structure[k] = '(';
        //     if(structure[k] == ')') structure[k] = ')';
        // }
        // exit(0);
        
        out3 << energy << "\t" << energy2 << "\t" << structure.compare(structure2) << endl;

        out3 << endl;
        if (energy==energy2) score++;
        // structure.compare(structure2) == 0
    
    }
    out3 << score << " out of " << size << endl;
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
    while(k<5){
        
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