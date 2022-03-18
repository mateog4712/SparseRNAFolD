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
vector<string> istructs;
vector<string> ostructs;
vector<double> energies;


ifstream in("/home/mgray7/sparsemfe/sparsemfefold/example.txt");

string str;
int i = 0;
while(getline(in,str)){
    // cout << str << endl;
    if(i%5 == 0){
        names.push_back(str);
    }
    else if(i%5==1){
        seqs.push_back(str);
    }
    else if(i%5==2){
        istructs.push_back(str);
    }
    else if(i%5==3){
        ostructs.push_back(str);
    }
    else if(i%5==4){
        energies.push_back(stod(str));
    }
    ++i;
}
in.close();

// cout << energies[0] << endl;

string command = "./build/src/SparseMFEFold ";

double score = 0;
int size = seqs.size();
for(int i =0;i<size;++i){
    string commands = command + "-r \"" + istructs[i] + "\" "+ seqs[i] +  " > out.txt";
    cout << i << endl;
    system(commands.c_str());
    ifstream in1("out.txt");
    getline(in1,str);
    getline(in1,str);
    in1.close();
    string structure = str.substr(0,seqs[i].length());
    double energy = stod(str.substr(seqs[i].length()+2,str.length()-1));
    if(energy == energies[i]) score+=1;
    else{
        cout << names[i] << endl;
        if(structure != ostructs[i]) cout << structure << endl << ostructs[i] << endl;
        if(energy != energies[i]) cout << energy << " != " << energies[i] << endl;
    }
}



cout << score << " out of " << size << endl;


return 0;
}


// TEST_CASE( "Structure is predicted", "[seq]"){

// string seq ="CCCAAACCAAAGGAAGGG";
// string structure = "(________________)";
// string restricted = "([[___[[___]]__]])";

// string seq2 = "AACUUGUCUUAGCUUUGCAGUCGAGUU";
// string structure2 = "(((.....................)))";
// string restricted2 = "((((((.((.........)).))))))";


// CHECK(iterativeFold(seq2,structure2) == restricted2);
// CHECK(iterativeFold(seq,structure) == restricted);




//}