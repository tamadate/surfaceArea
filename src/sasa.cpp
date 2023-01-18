#include "sasa.hpp"

void SASA::run(void){
    // For each atom, calculate SASA
    for (auto &a : atoms){
        getPair(a); // Get nearby atoms 

        double rij=a.radius+p.radius;
        int Naccess=0;

        // Count number of accessible positions on solvent surface
        for(int j = 0; j < Nprobe; j++){
            // Get a random point on a sphere surface centered at target with radius r
            double phi=4*PI*((double)rand()/RAND_MAX)-2*PI;
            double z=2*rij*((double)rand()/RAND_MAX)-rij;
            double rz=sqrt(rij*rij-z*z);
            double x=rz*cos(phi);
            double y=rz*sin(phi);
            p.pnt.x=x+a.pnt.x;
            p.pnt.y=y+a.pnt.y;
            p.pnt.z=z+a.pnt.z;
            if(accessible()) {
                Naccess++;
                checkProbePosition();
            }
        }

        // Calculate SASA (accessible / total)
        double ratio = float(Naccess)/Nprobe;
        double Stotal = 4*PI*rij*rij;
        double sasa = ratio*Stotal;
        total_sasa += sasa;
        a.sasa = sasa;
    }  

    // Print atomic statistics
    FILE* f=fopen("eachAtom.dat","w");
    for(auto &a : atoms) fprintf(f,"%d\t%s\t%f\n", a.id, a.residue.c_str(), a.sasa);
    fclose(f);
    printf("Total solvent accessible surface area: %.f\n", total_sasa); // Print totals
}