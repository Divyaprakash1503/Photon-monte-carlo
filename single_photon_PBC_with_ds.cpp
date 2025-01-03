
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <fstream>
using namespace std;
const double PI = 3.14159265358979323846;

// Random number generator in the range [0, 1)
double rand_double()
 {return rand() / (RAND_MAX + 1.0);}
// Generate a random number in a given range [min, max]
double rand_range(double min, double max)
{return min + (max - min) * rand_double();}

// Photon structure
struct Photon
{
    double x, y, z,u,v,w; // where u is velocity in x direction , v in y and , w in z.
    double energy;
                    };

void normalize (double& u, double& v, double& w)
{
    double magnitude = sqrt( u*u + v*v + w*w );
    u = u/magnitude;
    u=fabs(u);
    v /= magnitude;
    w /= magnitude;
}
 void random_direction(double& u, double& v, double& w)
 {
    double phi = rand_double()*PI -(PI/2); // doubt 1
    double  theta= 2 * PI * rand_double();

     u = cos(theta) * cos(phi);
     v =sin(theta) * cos(phi);
     w = sin(phi);

     normalize(u,v,w);
 }

int main()
{
srand(static_cast<unsigned int>(time(0))); // random number changes for every iteration
//const int num_photons = 100;
   double tow = 1;
   const double tol = 0.001;
   const double X=1;
   const double Y=1; // Geometry of problem
   const double Z=1;
   const double kappa = 1;//absorption coefficient
   int path_multiplier=0;
   double Ds=0;
   double S=0.005;
   // Initialize photon with random position and random direction in the range [0, 1]
        Photon photon = {0, rand_double(), rand_double(), 0.0, 0.0, 0.0, 1.0};
        random_direction(photon.u, photon.v, photon.w);
        cout << "starting ="<<photon.x<<","<<photon.y<<","<<photon.z<<endl;
        ofstream outfile("single_photon_PEB.txt");

            double step_size=0.005;

            for (double i =0 ;i <= 10; i+=step_size)
            {


                 //Boundry condition
                 if ( photon.z<=0 )

                {photon.z = 1;}  // bottom boundry

                 else if(photon.z>=Z)
            {
               photon.z = 0; }

              if ( photon.y<=0 )

                {photon.y = 1;}  // side boundry

                 else if(photon.y>=Y)
            {
               photon.y = 0; }


            if ( tow<tol || photon.x>=X || photon.x<0 )
            {
                break; // photon gets absorbed
            }

            else{
              // Update photon position
          photon.x += step_size * photon.u;
          photon.y += step_size * photon.v;
          photon.z += step_size * photon.w;
            }

//// calculation of path length
            if(photon.x >= S)
            {
                Ds = step_size * (path_multiplier+1);
                cout<<"DS for cell ="<<Ds<<endl;
                path_multiplier=0;
                S=S+step_size;

            }
            else{path_multiplier = path_multiplier +1;}

            double alpha = 1- exp(-kappa *(Ds));
                 tow = tow*(1-alpha);

           // random_direction(photon.u, photon.v, photon.w);

    cout <<  "Photon " << i<<","<<photon.x<<"  "<< ", " << photon.y << ", " << photon.z<<", "<<alpha<<","<<tow<<endl;
    outfile<<"Photon " << i<<","<<photon.x << ", " << photon.y << ", " << photon.z<<", "<<alpha<<","<<tow<<endl;
            }// end of for loop

         outfile.close();

}


