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
// Used to change fasta into individual files
int main(int argc,char **argv) {
    
    string file = "/home/mgray7/output2/fasta/rnasoft.txt";

    string file_M = "/home/mgray7/output2/linearfold/rnasoft/";
    string fileO = "";

    ifstream in(file);
    
    string str;
    int i = 0;
    vector<string> names;
    vector<string> seqs;
    vector<string> stru;
    while(getline(in,str)){
        if(i%3==0){
            string name = str.substr(1,str.length()-1);
            names.push_back(name);
        }
        if(i%3==1){
            seqs.push_back(str);
        }
        if(i%3 == 2){
            stru.push_back(str);
        }
        ++i;
    }
    in.close();

    for(int j=0;j<names.size();++j){
        fileO = file_M + names[j] + ".txt";
        ofstream out(fileO);
        out << ">" <<  names[j] << endl;
        out << seqs[j] << endl;
        out << stru[j] << endl;
        out.close();
    }


    return 0;
}