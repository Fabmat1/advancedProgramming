#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <fstream>
#include <iostream>

using namespace std;

struct GravitySimulation {
    // A struct used for N-body simulations

    double G;
    double Softening;
    int RandomSeed;
    vector<double> AccelerationTime;
    vector<double> x;
    vector<double> y;
    vector<double> z;
    vector<double> vx;
    vector<double> vy;
    vector<double> vz;
    vector<double> m;
    vector<double> ax;
    vector<double> ay;
    vector<double> az;
    vector<double> V;

    default_random_engine generator;
    uniform_real_distribution<double> distribution;

    GravitySimulation(double G=1.0, double Softening=0.05, int RandomSeed=643)
            : G(G), Softening(Softening), RandomSeed(RandomSeed), generator(RandomSeed), distribution(0.0, 1.0)
    {
    }

    // Begin function definitions:

    void CalcAcceleration() {
        auto t_start = chrono::high_resolution_clock::now();
        int Nparticles = x.size();
        vector<vector<double>> dx(Nparticles, vector<double>(Nparticles)),
                               dy(Nparticles, vector<double>(Nparticles)),
                               dz(Nparticles, vector<double>(Nparticles));

        #pragma omp parallel for
        for (int i = 0; i < Nparticles; ++i) {
            for (int j = 0; j < Nparticles; ++j) {
                dx[i][j] = x[i] - x[j];
                dy[i][j] = y[i] - y[j];
                dz[i][j] = z[i] - z[j];
            }
        }

        vector<vector<double>> r2_ij(Nparticles, vector<double>(Nparticles)),
                               r_ij(Nparticles, vector<double>(Nparticles)),
                               m2_ij(Nparticles, vector<double>(Nparticles)),
                               F_ij(Nparticles, vector<double>(Nparticles));

        #pragma omp parallel for
        for (size_t i = 0; i < Nparticles; ++i)
        {
            for (size_t j = 0; j < Nparticles; ++j)
            {
                r2_ij[i][j] = dx[i][j]*dx[i][j] + dy[i][j]*dy[i][j] + dz[i][j]*dz[i][j] + Softening*Softening;
                r_ij[i][j] = sqrt(r2_ij[i][j]);
                m2_ij[i][j] = m[i] * m[j];
                F_ij[i][j] = G * m2_ij[i][j] / r2_ij[i][j] / r_ij[i][j];
            }
        }

        ax.resize(Nparticles);
        ay.resize(Nparticles);
        az.resize(Nparticles);
        V.resize(Nparticles);

        #pragma omp parallel for
        for (int i = 0; i < F_ij.size(); i++)
        {
            double sumx = 0;
            double sumy = 0;
            double sumz = 0;
            double sumV = 0;

            for (int j = 0; j < F_ij[i].size(); j++)
            {
                sumx += F_ij[i][j] * dx[i][j];
                sumy += F_ij[i][j] * dy[i][j];
                sumz += F_ij[i][j] * dz[i][j];
                sumV -= G * m2_ij[i][j] / r_ij[i][j];
            }

            // This assumes m[i] != 0
            ax[i] = -sumx / m[i];
            ay[i] = -sumy / m[i];
            az[i] = -sumz / m[i];
            V[i] = sumV / m[i];
        }

        auto t_end = chrono::high_resolution_clock::now();
        AccelerationTime.push_back(chrono::duration<double>(t_end - t_start).count());
    }

    void InitializeHernquistHalo(double TotalMass=1.0, double ScaleRadius=1.0, int Nparticles=384) {
        x.resize(Nparticles);
        y.resize(Nparticles);
        z.resize(Nparticles);
        vx.resize(Nparticles);
        vy.resize(Nparticles);
        vz.resize(Nparticles);
        m.resize(Nparticles);
        ax.resize(Nparticles);
        ay.resize(Nparticles);
        az.resize(Nparticles);
        V.resize(Nparticles);

        vector<double> M_particles(Nparticles), r_particles(Nparticles), Theta_particles(Nparticles), Phi_particles(Nparticles);
        vector<double> EscapeVel_particles(Nparticles), v_particles(Nparticles);

        #pragma omp parallel for
        for (int i = 0; i < Nparticles; ++i) {
            M_particles[i] = distribution(generator) * TotalMass;
            r_particles[i] = ScaleRadius / (sqrt(TotalMass / M_particles[i]) - 1.0);
            Theta_particles[i] = acos(2.0 * distribution(generator) - 1.0);
            Phi_particles[i] = distribution(generator) * M_PI * 2;

            // set x, y, z
            x[i] = r_particles[i] * sin(Theta_particles[i]) * cos(Phi_particles[i]);
            y[i] = r_particles[i] * sin(Theta_particles[i]) * sin(Phi_particles[i]);
            z[i] = r_particles[i] * cos(Theta_particles[i]);

            // set velocities
            EscapeVel_particles[i] = sqrt(2 * G * M_particles[i] / r_particles[i]);
            v_particles[i] = distribution(generator) * 0.1 * EscapeVel_particles[i];
            Theta_particles[i] = acos(2.0 * distribution(generator) - 1.0);
            Phi_particles[i] = distribution(generator) * M_PI * 2;

            // set vx, vy, vz
            vx[i] = v_particles[i] * sin(Theta_particles[i]) * cos(Phi_particles[i]);
            vy[i] = v_particles[i] * sin(Theta_particles[i]) * sin(Phi_particles[i]);
            vz[i] = v_particles[i] * cos(Theta_particles[i]);
            m[i] = TotalMass/Nparticles;
        }

        CalcAcceleration();
    }
    void RunSimulation(double dt = 0.01, double tmax = 30.0, bool save_out=false) {
        // Leapfrog offset
        CalcAcceleration();
        #pragma omp parallel for
        for (int i = 0; i < vx.size(); i++) {
            vx[i] += ax[i] * dt / 2.0;
            vy[i] += ay[i] * dt / 2.0;
            vz[i] += az[i] * dt / 2.0;
        }

        // main loop
        double t = 0;
        int IntegerTimestep = 0;
        while (t < tmax) {
            #pragma omp parallel for
            for (int i = 0; i < x.size(); i++) {
                x[i] += vx[i] * dt;
                y[i] += vy[i] * dt;
                z[i] += vz[i] * dt;
            }

            CalcAcceleration();
            #pragma omp parallel for
            for (int i = 0; i < vx.size(); i++) {
                vx[i] += ax[i] * dt;
                vy[i] += ay[i] * dt;
                vz[i] += az[i] * dt;
            }

            if (IntegerTimestep % 10 == 0 && save_out){
                ofstream outfile;
                outfile.open("cpp_Nbody_coords.csv", ios_base::app); // append instead of overwrite
                for (int i = 0; i < x.size()-1; ++i) {
                    outfile << x[i] << ";" << y[i] << ";" << z[i] << ";";
                }
                outfile << x[x.size()-1] << ";" << y[x.size()-1] << ";" << z[x.size()-1];
                outfile << "\n";
            }
            t += dt;
            IntegerTimestep += 1;
        }
    }
};


int main(){
    if (FILE* file = fopen("cpp_Nbody_times.csv", "r")) {
        fclose(file); // Close the file handle
        remove("cpp_Nbody_times.csv"); // Delete the file
    }

    if (FILE* file = fopen("cpp_Nbody_coords.csv", "r")) {
        fclose(file); // Close the file handle
        remove("cpp_Nbody_coords.csv"); // Delete the file
    }

    cout << "Generating initial plot.." << endl;

    GravitySimulation grav = GravitySimulation();
    grav.InitializeHernquistHalo(1, 1, 512);
    grav.RunSimulation(0.01, 30, true);


    for (int i = 5; i < 12; ++i) {

        GravitySimulation grav = GravitySimulation();
        grav.InitializeHernquistHalo(1, 1, pow(2, i));
        grav.RunSimulation(0.01, 2);

        double mean_acc = 0;
        for (double i : grav.AccelerationTime) {
            mean_acc += i;
        }
        mean_acc /= grav.AccelerationTime.size();

        cout << "N = " << pow(2, i) << ": " << mean_acc << endl;

        ofstream outfile;
        outfile.open("cpp_Nbody_times.csv", ios_base::app); // append instead of overwrite
        outfile << mean_acc << "\n";
    }
    return 0;
}