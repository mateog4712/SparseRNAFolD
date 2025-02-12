# SparseRNAFolD

#### Description:
Software implementation of SparseRNAFolD.      
SparseRNAFolD is a time- and space-efficient sparsified minimum free energy folding of pseudoknot-free RNAs

#### Supported OS: 
Linux, macOS


### Installation:  
Requirements: A compiler that supports C++11 standard (tested with g++ version 4.9.0 or higher), Pthreads, and CMake version 3.1 or greater.    

[CMake](https://cmake.org/install/) version 3.1 or greater must be installed in a way that HFold can find it.    
To test if your Mac or Linux system already has CMake, you can type into a terminal:      
```
cmake --version
```
If it does not print a cmake version greater than or equal to 3.1, you will have to install CMake depending on your operating system.

[Linux instructions source](https://geeksww.com/tutorials/operating_systems/linux/installation/downloading_compiling_and_installing_cmake_on_linux.php)

#### Steps for installation   
1. [Download the repository](https://github.com/mateog4712/SparseRNAFolD.git) and extract the files onto your system.
2. From a command line in the root directory (where this README.md is) run
```
cmake -H. -Bbuild
cmake --build build
```   
If you need to specify a specific compiler, such as g++, you can instead run something like   
```
cmake -H. -Bbuild -DCMAKE_CXX_COMPILER=g++
cmake --build build
```   
This can be useful if you are getting errors about your compiler not having C++17 features.

Help
========================================

```
Usage: SparseRNAFolD[options] [input sequence]
```

Read input file from cmdline; predict minimum free energy and optimum structure using the time- and space-efficient MFE RNA folding algorithm.

```
  -h, --help             Print help and exit
  -V, --version          Print version and exit
  -v, --verbose          Turn on verbose output
  -m, --mark-candidates  Represent candidate base pairs by square brackets
  -r, --input-structure  Give a restricted structure as an input structure
  -d, --dangles=INT      How to treat \"dangling end\" energies for bases adjacent to helices in free ends and multi-loops (default=`2')
  -P, --paramFile        Read energy parameters from paramfile, instead of using the default parameter set.
      --noGC             Turn off garbage collection and related overhead
```

Remarks:
    The default parameter file is DP09. This can be changed via -P and specifying the parameter file you would like

#### Example:
    assume you are in the SparseRNAFolD directory
    ./build/src/SparseRNAFolD GGGGAAAACCCC
    ./build/src/SparseRNAFolD GGGGAAAACCCC -r "((........))" -d1

```
./build/src/SparseRNAFolD -m -v UAACUUAGGGGUUAAAGUUGCAGAUUGUGGCUCUGAAAACACGGGUUCGAA

UAACUUAGGGGUUAAAGUUGCAGAUUGUGGCUCUGAAAACACGGGUUCGAA
.[[(((..(..[((.[[[([[[...]]])]]].))]...)..)))]].... (-6.00)

TA cnt:165
TA max:167
TA av:167
TA rm:6

Can num:109
Can cap:118
TAs num:165
TAs cap:169
```
    
## Results
Results can be found at https://github.com/mateog4712/SparseRNAFolD-RawData

## Questions
For questions, you can email mateo2@ualberta.ca
