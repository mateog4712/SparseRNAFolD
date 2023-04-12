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
// Used to output files to bpseq
// Used to do comparisons between bpseqs

int main(int argc,char **argv) {
    
    string file = "/home/mgray7/output2/linearfold/fasta/rnasoft.txt";

    string file_M = "/home/mgray7/output2/linearfold/pred/";

    ifstream in(file);
    
    string str;
    int i = 0;
    vector<string> names;
    vector<string> seqs;
    while(getline(in,str)){
        if(i%2==0){
            string name = str.substr(1,str.length()-1);
            names.push_back(name);
        }
        if(i%2==1){
            seqs.push_back(str);
        }
        ++i;
    }
    in.close();

    string fileO = "";
    string program = "sparse";
    string program2 = "linear";
    string command = "./build/SparseMFEFold -d2 ";
    string command2p1 = "echo ";
    string command2p2 = " | ../../LinearFold/linearfold -V -d 2";
    for(int j=0;j<names.size();++j){
        // fileO = file_M + program + "/" + names[j] + ".txt";
        // ofstream out(fileO);
        
        
        // string commands = command + seqs[j] +  " > out.txt";
        // // cout << commands << endl;
        // // exit(0);
        // system(commands.c_str());

        // ifstream in1("out.txt");
        // getline(in1,str);
        // getline(in1,str);
        // in1.close();

        // string structure = str.substr(0,seqs[j].length());
        // double energy = stod(str.substr(seqs[j].length()+2,str.length()-1));

        // out << ">" <<  names[j] << endl;
        // out << seqs[j] << endl;
        // out << structure << endl;
        

        // out.close();

        string command2 = command2p1 + seqs[j] + command2p2 +  " > out.txt";
        system(command2.c_str());

        ifstream in2("out.txt");
        getline(in2,str);
        getline(in2,str);
        in2.close();

        string structure = str.substr(0,seqs[j].length());
        double energy = stod(str.substr(seqs[j].length()+2,str.length()-1));

        fileO = file_M + program2 + "/" + names[j] + ".txt";
        ofstream out1(fileO);

        out1 << ">" <<  names[j] << endl;
        out1 << seqs[j] << endl;
        out1 << structure << endl;

        out1.close();
    
    }

    return 0;
}

