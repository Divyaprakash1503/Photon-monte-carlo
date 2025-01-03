//working on this//success increase and decrease of N
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <ctime>
using namespace std;

const double PI = 3.14159265358979323846;
const int SE = 50;       // Number of surface elements
const int N = 200000;         // Number of photons emitted per segment
const double step_size = 0.05;
const double sigma = 5.67 * pow(10, -8);
const double R = 0.5;
const double A = PI * (2 * R);

// Random number generators
double rand_double() {
    return rand() / (RAND_MAX + 1.0);
}
double rand_range(double min, double max) {
    return min + (max - min) * rand_double();
}

// Photon structure to store the position and direction
struct Photon {
    double x, y, z;
    double u, v, w;
};

void Initialization(Photon &photon, double R, double segment_angle_start, double segment_angle_end, double E, double* Q) {
    photon.z = 0.0;
    double theta = rand_range(segment_angle_start, segment_angle_end);
    photon.x = R * cos(theta);
    photon.y = R * sin(theta) + 1.0;

    // Lambertian emission for direction
    double phi = 2 * PI * rand_double();
    double sin_theta = sqrt(rand_double());
    double cos_theta = sqrt(1 - cos_theta * cos_theta);

    photon.u = sin_theta * cos(phi);
    photon.v = sin_theta * sin(phi);
    photon.w = cos_theta;

    //cout << "Starting Position: (" << photon.x << ", " << photon.y << ", " << photon.z << ") "<<endl;
         //<< "Direction: (" << photon.u << ", " << photon.v << ", " << photon.w << ")\n";

    for (double k = 0; k <= 100; k += step_size) 
    {
        photon.x += step_size * photon.u;
        photon.y += step_size * photon.v;
        //cout << "Photon Position at step " << k << ": (" << photon.x << ", " << photon.y << ", " << photon.z << ")\n";

        double IF = sqrt(pow(photon.x, 2) + pow(photon.y - 1, 2));
        if (IF <= R) {
           // cout << "Photon absorbed inside radius.\n";
            break;
        }

        if (photon.y <= 0 && (photon.x >= -1.5 && photon.x <= 1.5)) 
        {
            //cout << "Energy E = " << E << endl;
            int cell = 1;
            for (double L = -1.5; L <= 3; L += 0.025)
             {
                if (photon.x >= L && photon.x < L + 0.025) {
                    Q[cell] += E;
                    break;
                }
                cell++;
            }

            break;
        } else if (photon.y >= 3 || (photon.x >= 1.5 || photon.x <= -1.5)) {
           // cout << "Photon out of bounds.\n";
            break;
        }
    }
}

int main() {
    srand(static_cast<unsigned int>(time(0)));

    const double segment_angle = 2 * PI / SE;
    double Q[121] = {0};  // Energy array

    for (int i = 1; i <= SE; ++i) {
        double segment_angle_start = i * segment_angle;
        double segment_angle_end = (i + 1) * segment_angle;

        for (int j = 0; j < N; ++j) {
            Photon photon;
            double T = 1000;
            double E = (A * 1) / N;
            Initialization(photon, R, segment_angle_start, segment_angle_end, E, Q);
        }
        cout<<"segment"<<"   "<<i<<"   "<<"done"<<endl;
    }

    for (int I = 1; I <= 120; I++) {
        cout << "Q "<< -1.5+(I*0.025) <<"    " << Q[I] << endl;
    }

    return 0;
}

