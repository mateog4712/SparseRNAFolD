#include <string>
#include <fstream>
#include <vector>
#include <iostream>
//#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

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

string command = "./build/src/SparseMFEFold -m ";
string command2 = "/home/mgray7/RNAfold/build/src/RNAFold ";
double score = 0;
int size = seqs.size();
// exit(0);
ofstream out3("results.txt");
for(int i =0;i<size;++i){
    string commands = command + seqs[i] +  " > out.txt";
    system(commands.c_str());

    ifstream in1("out.txt");
    getline(in1,str);
    getline(in1,str);
    in1.close();

    string structure = str.substr(0,seqs[i].length());
    double energy = stod(str.substr(seqs[i].length()+2,str.length()-1));
    
    ofstream out1("/home/mgray7/RNAfold/in.txt");
    out1 << names[i] << endl;
    out1 << seqs[i] << endl;
    out1.close();
    string infile = "/home/mgray7/RNAfold/in.txt";
    
    string commands2 = command2 + "< " + infile + " > /home/mgray7/RNAfold/out.txt";
    system(commands2.c_str());

    ifstream in2("/home/mgray7/RNAfold/out.txt");
    getline(in2,str);
    getline(in2,str);
    getline(in2,str);
    in2.close();

    string structure2 = str.substr(0,seqs[i].length());
    double energy2 = stod(str.substr(seqs[i].length()+2,str.length()-1));
    out3 << names[i] << endl;
    out3 << seqs[i] << endl;
    out3 << structure << endl;
    out3 << structure2 << endl;
    for(int k=0;k<structure.length();++k){
        if(structure[k] == '[') structure[k] = '(';
        if(structure[k] == ')') structure[k] = ')';
    }
    out3 << energy << "\t" << energy2 << "\t" << structure.compare(structure2) << endl;

    out3 << endl;
    if (energy==energy2) score++;
    // structure.compare(structure2) == 0
    
}
out3 << score << " out of " << size << endl;
out3.close();
return 0;
}