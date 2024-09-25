# PythiaJet

Project for the coursework on selected topics in heavy-ion physics 

The project is based on pp event generation in 5 TeV in PYTHIA and FASTJET clustering

Link to the packages:

1. PYTHIA 8.312 [https://www.pythia.org/download/pythia83/pythia8312.tgz]
2. FASTJET 3.4.3 [https://fastjet.fr/repo/fastjet-3.4.3.tar.gz]
3. ROOT 3.32.00 [https://github.com/root-project/root/tree/v6-32-00-patches]

How to run on Mac: 
clang++ -std=c++17 -o <output> <input>.cc -I/pythia/pythia8312/include -L/pythia/pythia8312/lib -lpythia8 `root-config --cflags --glibs` `/fastjet/fastjet-install/bin/fastjet-config --cxxflags --libs --plugins`

