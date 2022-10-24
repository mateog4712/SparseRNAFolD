#include <string>
#include <fstream>
#include <vector>
#include <iostream>
//#define CATCH_CONFIG_MAIN
#include "catch.hpp"

using namespace std;

int main(int argc,char **argv) {

cout << "Put file: ";
string filename;
cin >> filename;


string file = "/home/mgray7/output2/fasta/" + filename + ".txt";

string str;

ifstream in(file);

int i=0;
int num = 0;
int max_length=0;
while(getline(in,str)){
++i;
if(i%2==0){
++num;
int length = str.length();
if(length>max_length) max_length = length;

}
}
in.close();


cout << "number of sequences = " << num << endl << "max length is " << max_length << std::endl;




return 0;
}