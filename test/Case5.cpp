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
// Used to change files from rnastrand into readable format
int main(int argc,char **argv) {
    
    string file = "/home/mgray7/output2/fasta/rnasoft.txt";

    string fileO = "/home/mgray7/output2/fasta/rnasoftO.txt";

    ifstream in(file);
    ofstream out(fileO);
    string str;
    string seq;
    string structure;
    while(getline(in,str)){
        if(str[2] == 'F'){
            // cout << str << endl;
            // exit(0);
            string name = str.substr(7,str.length()-7-3);
            out << ">" << name << endl;
            getline(in,str);
            getline(in,str);
            getline(in,str);
        }
        else{

            if(str == "" && seq != "" && structure!=""){
                transform(seq.begin(), seq.end(), seq.begin(), ::toupper);
                out << seq << endl;
                out << structure << endl;
                seq = "";
                structure = "";
                // out << endl;
            }
            else{
                if(str[0] == '(' || str[0] == ')' || str[0] == '.'){
                    structure = structure + str;
                } else if(str== "") {
                }
                else {
                    seq = seq + str;
                }
                
                
            }
        }
    }
    in.close();
    out.close();


    return 0;
}