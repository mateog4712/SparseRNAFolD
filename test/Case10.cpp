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

//Take all time and memory numbers from the families and compile them into a single file
using namespace std;

int main(int argc,char **argv) {
    vector<string> names;
    vector<string> seqs;

    string version = "";
    cout << "Put folder: ";
    cin >> version;

    string start = "/home/mgray7/output2/comparison/";
    string folder = start + version + "/";

    string file1 = folder + "ase.txt";
    string file2 = folder + "crw.txt";
    string file3 = folder + "pdb.txt";
    string file4 = folder + "rfa.txt";
    string file5 = folder + "srp.txt";
    string file6 = folder + "tmr.txt";

    
    string fileOT = folder + "time.txt";
    string fileOM = folder + "memory.txt";

    ofstream out1(fileOT);
    ofstream out2(fileOM);
    
    ifstream in1(file1);
    ifstream in2(file2);
    ifstream in3(file3);
    ifstream in4(file4);
    ifstream in5(file5);
    ifstream in6(file6);


    string str;
    while(getline(in1,str)){
        istringstream iss(str);
        string time; 
        iss >> time;
        string memory; 
        iss >> memory;
        out1 << time << endl;
        out2 << memory << endl;
    }  
    in1.close();
     while(getline(in2,str)){
        istringstream iss(str);
        string time; 
        iss >> time;
        string memory; 
        iss >> memory;
        out1 << time << endl;
        out2 << memory << endl;
    }  
    in2.close();
     while(getline(in3,str)){
        istringstream iss(str);
        string time; 
        iss >> time;
        string memory; 
        iss >> memory;
        out1 << time << endl;
        out2 << memory << endl;
    }  
    in3.close();
     while(getline(in4,str)){
        istringstream iss(str);
        string time; 
        iss >> time;
        string memory; 
        iss >> memory;
        out1 << time << endl;
        out2 << memory << endl;
    }  
    in4.close();
     while(getline(in5,str)){
        istringstream iss(str);
        string time; 
        iss >> time;
        string memory; 
        iss >> memory;
        out1 << time << endl;
        out2 << memory << endl;
    }  
    in5.close();
     while(getline(in6,str)){
        istringstream iss(str);
        string time; 
        iss >> time;
        string memory; 
        iss >> memory;
        out1 << time << endl;
        out2 << memory << endl;
    }  
    in6.close();
    out1.close();
    out2.close();

    
    return 0;
}