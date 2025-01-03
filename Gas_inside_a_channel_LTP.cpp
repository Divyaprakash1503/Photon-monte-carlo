//gas inside a channel //final// Linear temperature profile
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;

const double PI = 3.14159265358979323846;

// Random number generator in the range [0, 1)
double rand_double() {
    return rand() / (RAND_MAX + 1.0);
}

// Generate a random number in a given range [min, max]
double rand_range(double min, double max) {
    return min + (max - min) * rand_double();
}

// Normalize direction vector (u, v, w)
void normalize(double& u, double& v, double& w) {
    double magnitude = sqrt(u * u + v * v + w * w);
    u /= magnitude;
    v /= magnitude;
    w /= magnitude;
}

// Generate a random direction for the photon
void random_direction(double& u, double& v, double& w) {
    double theta = acos(1 - (2 * rand_double()));
    double phi = 2 * PI * rand_double();
    u = sin(theta) * cos(phi);
    v = sin(theta) * sin(phi);
    w = cos(phi);
    normalize(u, v, w);
}

void random_Sdirection(double& u, double& v,double& w)
{
   double theta = acos(1 - (2 * rand_double()));
   double phi =asin(sqrt(rand_double()));
   u = sin(theta) * cos(phi);
   u=fabs(u);
   v = sin(theta) * sin(phi);
   w = cos(phi);
   normalize(u, v, w);
   }

   void random_SRdirection(double& u, double& v,double& w)
{
   double theta = acos(1 - (2 * rand_double()));
   double phi =asin(sqrt(rand_double()));
   u = sin(theta) * cos(phi);
   u=fabs(u);
   u=-u;
   v = sin(theta) * sin(phi);
   w = cos(phi);
   normalize(u, v, w);
   }
// Photon structure
struct Photon {
    double x, y, z, u, v, w;
    double energy;
};

const double kappa = 1; // Absorption coefficient
const double step_size=0.005;
const double sigma = 5.67 * pow(10, -8); // Stefan-Boltzmann constant
const double X = 1, Y = 1, Z = 1;        // Geometry of problem
double vol = (Z * Y * step_size);
double area=(Y*Z);

// Photon tracking function
void track_photon(Photon& photon, int& cell, int& path_multiplier, double& E0, double& tow, double* E_Abs, double& cds, double& B, double& B0)
{
    double Ds = 0;
    double oldx, oldy, oldz;
    double x, y, z;

    for (double i = 0; i <= 100; i += step_size) {
        // Boundary condition handling
        if (photon.z <= 0) photon.z = 1;
        else if (photon.z >= 1) photon.z = 0;

        if (photon.y >= 1) photon.y = 0;
        else if (photon.y <= 0) photon.y = 1;

        // If photon is absorbed or goes out of bounds
        if (tow < 0.000001 || photon.x > 1 || photon.x < 0) {
            break; // Photon gets absorbed
        }

        // Update photon position
        oldx = photon.x;
        oldy = photon.y;
        oldz = photon.z;
        photon.x += step_size * photon.u;
        photon.y += step_size * photon.v;
        photon.z += step_size * photon.w;

      //  cout << "Photon " << i << ", X = " << photon.x << ", Y = " << photon.y << ", Z = " << photon.z << endl;

        // Backward photon calculation
        if (photon.u < 0) {
            if (photon.x <= B0) {
                Ds = step_size * (path_multiplier + 1) + cds;
                double IF = (oldx - B0) / (oldx - photon.x);
                x = B0;
                y = oldy - IF * (oldy - photon.y);
                z = oldz - IF * (oldz - photon.z);
                cds = sqrt(pow(photon.x - x, 2) + pow(photon.y - y, 2) + pow(photon.z - z, 2));
                Ds -= cds;
               // cout << "DS for cell = " << Ds << endl;
                path_multiplier = 0;
                B0 -= step_size;

                double alpha = 1 - exp(-kappa * (Ds));
                tow *= (1 - alpha);
                double E1 = E0 * (1 - alpha);
                double delE = E0 - E1;
                //cout << "tow = " << tow << endl;
               // cout << "E0 = " << E0 << "   E1 = " << E1 << endl;
               // cout << "delE = " << delE << endl;
                E0 = E1;

                E_Abs[cell] = delE + E_Abs[cell];
                cell--;
            } else {
                path_multiplier++;
            }
        } else {
            // Forward photon calculation
            if (photon.x >= B) {
                Ds = step_size * (path_multiplier + 1) + cds;
                double IF = (B - oldx) / (photon.x - oldx);
                x = B;
                y = oldy + IF * (photon.y - oldy);
                z = oldz + IF * (photon.z - oldz);
                cds = sqrt(pow(photon.x - x, 2) + pow(photon.y - y, 2) + pow(photon.z - z, 2));
                Ds -= cds;
               // cout << "DS for cell = " << Ds << endl;
                path_multiplier = 0;
                B += step_size;

                double alpha = 1 - exp(-kappa * (Ds));
                tow *= (1 - alpha);
                double E1 = E0 * (1 - alpha);
                double delE = E0 - E1;
                //cout << "tow = " << tow << endl;
               // cout << "E0 = " << E0 << "   E1 = " << E1 << endl;
               // cout << "delE = " << delE << endl;
                E0 = E1;

                E_Abs[cell] = delE + E_Abs[cell];
                cell++;
            } else {
                path_multiplier++;
            }
        }
    }
}

int main() {
    srand(static_cast<unsigned int>(time(0)));


    double E_Abs[201] = {0};
    double E_Emit[201] = {0};
    double Q[201] = {0};
    int n = 10; // Number of photons per cell

    // Pre-calculate emission energy for each cell
    int count = 1;
    for (double J = step_size / 2; J <= 1; J += step_size) {
        double T = 500 + 1000*J;
        E_Emit[count] = (4 * kappa * sigma * pow(T, 4) * vol) ;
        count++;
    }

    // Output the emitted energy
    for (int J1 = 1; J1 <= 200; J1++)
    {
        cout << "E_Emit of " << J1 << " = " << E_Emit[J1] << endl;
    }

    // Loop through each cell
    for (int S = 1; S <= 200; S++)
    {
        cout << "Cell = " << S << endl;

        double Boundry_0 = (S * step_size) - step_size;
        double Boundry = S * step_size;
        int cell = S;
        int N = S - 1;
	double P= ((N*0.005)+(cell*0.005))/2;
        // Loop through each volumetric photon in the cell
        for (int j = 1; j <= n; j++)
        {
            cout << "Photon No = " << j << endl;

            // Initialize photon with random position and direction
            Photon photon = {rand_range( (N * 0.005), 0.005 + (N * 0.005)), rand_double(), rand_double(), 0.0, 0.0, 0.0, 1.0};
            random_direction(photon.u, photon.v, photon.w);
             double B0=Boundry_0;
             double B=Boundry;
		cell=S;

            double T = 500 + 1000*P;
            double E0 = (4 * kappa * vol * sigma * pow(T, 4)) / n;
            double tow = 1;
            int path_multiplier = 0;
            double cds = 0;

            // Track photon
            track_photon(photon, cell, path_multiplier, E0, tow, E_Abs, cds, B, B0);

            //cout << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," << endl;
        }

    }

    for (int j1 = 1; j1 <= n; j1++)
        {
            cout << "surface left wall Photon No = " << j1 << endl;
            // Initialize surface photon from left wall with random position and direction
            Photon SLphoton = {0, rand_double(), rand_double(), 0.0, 0.0, 0.0, 1.0};
            random_Sdirection(SLphoton.u, SLphoton.v, SLphoton.w);
            double B0=0;
            double B=step_size;
            int cell=1;
            double T = 500;
            double E0 = (area * sigma * pow(T, 4)) / n;
            double tow = 1;
            int path_multiplier = 0;
            double cds = 0;
             // Track photon
            track_photon(SLphoton, cell, path_multiplier, E0, tow, E_Abs, cds, B, B0);
           // cout << ",,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,," << endl;

           // Initialize surface photon  from right wall with random position and direction
            cout << "surface right wall Photon No = " << j1 << endl;
            Photon SRphoton = {1, rand_double(), rand_double(), 0.0, 0.0, 0.0, 1.0};
            random_SRdirection(SRphoton.u, SRphoton.v, SRphoton.w);
              B0=1-step_size;
               B=1;
		 cell=200;
             double TR = 1500;
              E0 = (area * sigma * pow(TR, 4)) / n;
              tow = 1;
              path_multiplier = 0;
              cds = 0;
             // Track photon
            track_photon(SRphoton, cell, path_multiplier, E0, tow, E_Abs, cds, B, B0);
            cout << ":::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::" << endl;
            }

    // Print energy absorption for each cell
    for (int I = 1; I <= 200; I++) {
        cout << "E_abs[" << I << "] = " << E_Abs[I] << endl;
    }

    // Calculate and print the source term for each cell
    for (int m = 1; m <= 200; m++) {
        Q[m] = (E_Emit[m] - E_Abs[m]) / vol;
        cout << "Q[" << m << "] = " << Q[m] / 1000 << endl;
    }


    ofstream outfile("LP100.txt");
    for (int m1=1;m1<=200;m1++) //cell
    {

outfile<<m1*0.005<<" "<<Q[m1]/1000<<endl;

    }
 outfile.close();
 return 0;
}

