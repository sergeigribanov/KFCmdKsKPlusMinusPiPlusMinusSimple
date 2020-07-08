# KFCmdKsKpi
A simple (pi+pi-)(K+-pi-+) hypothesis test (2 vertices).

## Dependencies
1. Eigen3: http://eigen.tuxfamily.org/index.php?title=Main_Page
2. ROOT: https://root.cern.ch
3. CCGO: https://github.com/sergeigribanov/ccgo
4. KFBase: https://github.com/sergeigribanov/KFBase
5. KFCmd: https://github.com/sergeigribanov/KFCmd

## Building and running

To build this ROOT scripts you need to do following steps:
1. Build and install all dependencies.
2. Copy rootlogon.C form KFCmd installation directory (see env subdirectory) to workin directory.
3. Copy TrPh.h and TrPh.C files to working directory
4. Run CERN ROOT: root -l
5. .L TrPh.C++
6. Open an input file with the tr_ph tree: auto fl = TFile::Open(..., "read")
7. TrPh x(tr_ph)
8. x.Loop("output_file.root", 1.3) // 1.3 is a value of a magnetic field
