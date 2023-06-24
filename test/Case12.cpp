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

 // run results for file
using namespace std;

int main(int argc,char **argv) {

    string fileI = "/home/mgray7/output2/linearfold/fasta/rnasoft_non_redundant_removed.txt";
    ifstream in(fileI);
    string folderO = "/home/mgray7/output2/linearfold/pred/sparse2/";
    string str;
    int i =0;
    string name;
    while(getline(in,str)){
        if(i%2 == 0){
            name = str.substr(1,str.length());
        } else{
            string fileO = folderO + name + ".txt";
            ofstream out(fileO);
            string command = "./build/SparseRNAFolD -d1 " + str + " > out.txt";
            system(command.c_str());
            ifstream in2("out.txt");
            string str2;
            getline(in2,str2);
            string seq = str2;
            getline(in2,str2);
            string struc = str2;
            struc = struc.substr(0,seq.length());
            in2.close();
            out << ">" <<name << endl;
            out << seq << endl;
            out << struc << endl;
            out.close();

        }
        ++i;
    }
    in.close();
    return 0;
}