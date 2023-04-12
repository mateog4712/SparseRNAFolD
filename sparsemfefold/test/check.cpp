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

int main(int argc,char **argv) {
    string family;

    string file = "/home/mgray7/output2/results/";
    cout << "Put Family: ";
    cin >> family;
    file = file + family + ".txt";
    ifstream in(file);

    string str;
    int i = 0;
    while(getline(in,str)){
    // cout << str << endl;
        if(i%7 == 5){
            std::vector<std::string> vec;
            
            // std::cout << str << std::endl;
            istringstream iss(str);
            std::string T;
            int energy1;
            int energy2;
            int wrong;
            int j = 0;
            while(getline(iss,T,'\t')){
                // std::cout << T << std::endl;
                vec.push_back(T);
            }
            energy1 = atoi(vec[0].c_str());
            energy2 = atoi(vec[1].c_str());
            std::cout << i/7 << std::endl;
            // std::cout << energy1 << " " << energy2 << std::endl;
            if(energy1 != energy2){
                std::cout << "At " << i/7 << std::endl;
            }
            // exit(0);
        }
        ++i;
    }  
    in.close();

   
}



