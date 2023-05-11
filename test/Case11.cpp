#include <string>
#include <fstream>
#include <sstream>
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

 // Make file containing lengths of sequences
using namespace std;

int main(int argc,char **argv) {
    vector<string> names;
    vector<string> seqs;

    string file = "/home/mgray7/output2/fasta/";

    string family;

    cout << "Put family: ";
    cin >> family;

    string filename = file + family + ".txt";

    ifstream in(filename);
    string str;
    int i = 0;
    while(getline(in,str)){
        if(i%2 == 0){
            names.push_back(str);
        }
        else if(i%2==1){
            seqs.push_back(str);
        }
        ++i;
    }
    in.close();

    cout << seqs.size() << endl;
    string output = "/home/mgray7/output2/lengths/" + family + ".txt";
    ofstream out(output);

    for(i = 0;i<seqs.size();++i){
        out << seqs[i].length() << endl;
    }
    out.close();
    return 0;
}