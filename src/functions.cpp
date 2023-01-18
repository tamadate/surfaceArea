#include "sasa.hpp"

// Deterimine probe is aaccessible
bool SASA::accessible(void)
{
    for (auto &a : pair){
        double d1=dist(p.pnt,a.pnt);
        double d2=p.radius+a.radius;
        if(d1<d2) return false;
    }
    return true;
}


// Find all atoms within radius d of target
void SASA::getPair(Atom target)
{
    pair.clear();
    double d = target.radius + 2*p.radius + 2;
    for (auto &a : atoms){
        if(dist(target.pnt, a.pnt) < d){
            if(target.id == a.id) continue;
            pair.push_back(a);
        }
    }
}


// distance calculation
double SASA::dist(XYZ a, XYZ b){
    double dx=a.x-b.x;
    double dy=a.y-b.y;
    double dz=a.z-b.z;
    double r2=dx*dx+dy*dy+dz*dz;
    return sqrt(r2);
}


// distance calculation
void SASA::checkProbePosition(void){
    FILE* f=fopen("probePosition.dat","a");
    fprintf(f,"%e\t%e\t%e\n", p.pnt.x, p.pnt.y, p.pnt.z);
    fclose(f);
}


// Extract one Atom struct for each atom in the PDB file
void SASA::trim(string& s) {
    string::size_type pos = s.find_last_not_of(' ');
    if(pos != string::npos) {
        if (s.length() != pos + 1)
            s.erase(pos + 1);
            pos = s.find_first_not_of(' ');
        if(pos != 0)
            s.erase(0, pos);
    }
    else s="";
}

void SASA::load_pdb_file(void)
{
    map<string, float> radii;
    radii["C"] = 1.7;
    radii["H"] = 1.2;
    radii["O"] = 1.52;
    radii["N"] = 1.55;
    radii["S"] = 1.8;

    string line;
    ifstream stream(infile);
	string str;
	while(getline(stream,line)) {
        if(line.substr(0,4) != "ATOM") continue;
        struct Atom atom;
        atom.id = atoi(line.substr(6,5).c_str());
        atom.atom_name = line.substr(12,3);
        atom.residue   = line.substr(17,3);
        atom.aa_id     = atoi(line.substr(22,4).c_str());
        atom.pnt.x     = atof(line.substr(30, 8).c_str());
        atom.pnt.y     = atof(line.substr(38, 8).c_str());
        atom.pnt.z     = atof(line.substr(46, 8).c_str());
        trim(atom.atom_name);
        trim(atom.residue);
        atoms.push_back(atom);
    }
}



void SASA::load_my_file(void)
{
    ifstream stream(infile);
	string str;
	while(getline(stream,str)) {
        if (str.length()==0) continue;
        string tmp;
		istringstream stream2(str);
		// Atomic
        int loop=0;
        Atom atom;
        while(getline(stream2,tmp,'\t')) {
            if (loop==0) atom.id=stoi(tmp);
            if (loop==1) atom.pnt.x=stod(tmp);
            if (loop==2) atom.pnt.y=stod(tmp);
            if (loop==3) atom.pnt.z=stod(tmp);
            if (loop==4) atom.atom_name=tmp;
            if (loop==5) atom.radius=stod(tmp);
            if (loop==6) atom.residue=tmp;
            if (loop==7) atom.aa_id=stod(tmp);
            loop++;
        }
        atoms.push_back(atom);
    }
}
