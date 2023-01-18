#include <math.h>
#include <map>
#include <iostream>
#include <string>
#include <vector>
#include <stdio.h>
#include <cstdlib>
#include <time.h>
#include <fstream>
#include <sstream>

using namespace std;
#pragma once

const double PI = acos(-1);

struct XYZ { 
    double x;
    double y;
    double z; 
};
struct Atom {
    XYZ pnt;
    string residue;
    string atom_name;
    int id;
    int aa_id;
    double radius;
    double sasa;
};

struct Probe {
    XYZ pnt;
    double radius;
};


class SASA{
    public:

    vector<struct Atom> atoms;  // array of atoms
    vector<struct Atom> pair;   // pair list
    Probe p;
    double total_sasa;  // total solvent accessible surface area
    char* infile;   // input data file name
    int Nprobe; // Number of random probe positions

    // functions
    double dist(XYZ a, XYZ b);  
    bool accessible(void);
    void trim(string& s);
    void load_pdb_file(void);
    void load_my_file(void);
    void getPair(Atom target);
    void checkProbePosition(void);
    void run(void);

    SASA(char* file){
        total_sasa  = 0;// Initialize solvent accessible surface area variable
        infile=file;
        load_my_file();
        p.radius=1.4; // probe default radius
        Nprobe=2000;
    };
    ~SASA(void){};
};




