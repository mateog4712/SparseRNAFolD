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

void getResults(string program){

string fileI = "../../output2/linearfold/results/" + program + ".txt";

ifstream in(fileI);
double f_measure = 0;
double sensitivity = 0;
double ppv = 0;
int count = 0;
string str;


getline(in,str);
while(getline(in,str)){
if(str == "" || str[0] == ' ') continue;
istringstream ss(str);
string nothing = "";
double i = 0;
ss >> nothing;
ss >> i;
if(i != -1){
f_measure+=i;
ss >> i;
sensitivity+=i;
ss >> i;
ppv+=i;
count+=1;
}
}
f_measure/=count;
sensitivity/=count;
ppv/=count;
string file = "../../output2/linearfold/results/" + program  + "O.txt";


ofstream out(file,ofstream::app);
out << f_measure << "\t" << sensitivity << "\t" << ppv << endl;
out.close();

}

void runFmeasure(string program){
    /** Updates the perl file for BPseq and runs it on the files from the family inputted **/
    string line;
    string folderI;
    string folderO;
    string dirN; 

    string calcf;
    string calcfU;
    
    calcf = "../../calculate_f_measure.txt";
    calcfU = "../../calculate_f_measureU.pl";
    
   
    ifstream cf(calcf.c_str());
    ofstream cfU(calcfU.c_str(), ofstream::trunc);

    // Starts a counter for reading the file
    int counter = 1;
    while (getline(cf, line)) {
        // If we've reached the 7th line of the txt file, place the line which look at the correct folder
        if (counter == 98) {
            cfU << "$crr_path = \"/home/mgray7/output2/linearfold/corrbp" << "/\";" << endl;
        } // The same is done with the 8th line. It is replaced with the correct folder location
        else if (counter == 99) {
            //change this line to fix problem
            cfU << "$chk_path = \"/home/mgray7/output2/linearfold/predbp/" << program << "/\";" << endl;
        } else if (counter == 152) {
            cfU << "$experiment = \"/home/mgray7/output2/linearfold/results/" << program << ".txt\";" << endl;
        } else {
            // Send the exact line from the txt file into the perl file
            cfU << line << endl;
        }
        // Increment the counter
        counter++;

    }
    // Closes the files
    cf.close();
    cfU.close();

    // // Runs the updated perl file
    string command;
     command = "perl ../../calculate_f_measureU.pl";
    
    system(command.c_str());

}


void runPerl(string program) {
/** Updates the perl file for BPseq and runs it on the files from the family inputted **/
    string line;
    string folderI;
    string folderO;
    string dirN;

    // Runs the following parts twice: once for predictions and once for the correct info
    int i = 0;
    if(program == "rnasoft") i = 1;
    // Takes the converted perl file and opens it for reading and a new perl file for writing
    // When done as a perl file, the getline function messes up and skips lines
    string bpseq = "../../bpseq.txt";
    string bpseqU = "../../bpseqU.pl";
    ifstream bp(bpseq.c_str());
    ofstream bpU(bpseqU.c_str(), ofstream::trunc);

    // defines the folders that are going to be outputted to
    if (i == 0) {
        folderI = "pred";
        folderO = "predbp";
    } else {
        folderI = "rnasoft";
        folderO = "corrbp";
    }

    // Starts a counter for reading the file
    int counter = 1;
    while (getline(bp, line)) {
        // If we've reached the 7th line of the txt file, place the line which look at the correct folder
        if (counter == 7) {
            if(i==0){
                bpU << "$inputpath = \"/home/mgray7/output2/linearfold/" << folderI << "/" << program << "/\";" << "#input path" << endl;

            }
            else{
                bpU << "$inputpath = \"/home/mgray7/output2/linearfold/" << folderI << "/\";" << "#input path" << endl;

            }
        } // The same is done with the 8th line. It is replaced with the correct folder location
        else if (counter == 8) {
            if(i==1){
                bpU << "$outputpath = \"/home/mgray7/output2/linearfold/" << folderO << "/\";" << "#output path" << endl;
            }
            else {
                bpU << "$outputpath = \"/home/mgray7/output2/linearfold/" << folderO << "/" << program << "/\";" << "#output path" << endl;
            }
        } else {
            // Send the exact line from the txt file into the perl file
            bpU << line << endl;
        }
        // Increment the counter
        counter++;
    }
    // Close the files
    bp.close();
    bpU.close();
    // Run the updated perl file
    string command = "perl ../../bpseqU.pl";
    system(command.c_str());
    
}


int main(int argc,char **argv) {
    
    string program = "sparse";
    // runPerl(program);

    // runFmeasure(program);
    getResults(program);



    return 0;
}

