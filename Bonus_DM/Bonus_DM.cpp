// I chose to ignore the instructions for task 3ff, lets see how it works out
#include <iostream>
#include <random>
#include <chrono>
#include <immintrin.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <tuple>


using namespace std;


tuple<vector<double>, vector<double>, vector<double>, vector<double>> readDataFromFile(const string& filename, const string& separator) {

    vector<double> col1, col2, col3, col4;
    ifstream file(filename);
    string line;

    while (getline(file, line)) {
        if (line.empty() || line[0] == '#'){
            continue;
        }

        istringstream iss(line);
        string token;
        vector<string> tokens;

        while (getline(iss, token, separator[0])) {
            tokens.push_back(token);
        }

        if (tokens.size() == 4) {
            col1.push_back(stod(tokens[0]));
            col2.push_back(stod(tokens[1]));
            col3.push_back(stod(tokens[2]));
            col4.push_back(stod(tokens[3]));
        }
    }

    return make_tuple(col1, col2, col3, col4);
}


vector<double> generate_linspace(int length, double start, double end){
    vector<double> result(length);

    // Calculate the increment value
    double increment = (end - start) / (length - 1);

    // Vectorize the loop using AVX2
    #pragma omp parallel for
    for (int i = 0; i < length; i += 4) {
        __m256d indices = _mm256_set_pd(i + 3, i + 2, i + 1, i);
        __m256d increments = _mm256_set1_pd(increment);
        __m256d values = _mm256_fmadd_pd(indices, increments, _mm256_set1_pd(start));

        _mm256_storeu_pd(&result[i], values);
    }

    return result;
}


vector<double> mass_shells(vector<double> x,vector<double> y,vector<double> z,vector<double> m, vector<double> r){
    vector<double> result(r.size(), 0.0);
    double dist;

    result[0] = 0.0;

    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i) {
        dist = sqrt(x[i]*x[i]+y[i]*y[i]+z[i]*z[i]);
        for (int j = 0; j < r.size(); ++j) {
            if (r[j] > dist){
                #pragma omp atomic
                result[j] += m[i];
                break;
            }
        }
    }

    return result;
}


vector<double> volume_shells(vector<double> r){
    vector<double> result(r.size(), 0.0);

    result[0] = 0.0;

    #pragma omp parallel for
    for (int i = 1; i < r.size(); ++i) {
        result[i] = 4.0/3.0*M_PI*(r[i]*r[i]*r[i]-r[i-1]*r[i-1]*r[i-1]);
    }

    return result;
}


void save_file(vector<double> vector1, vector<double> vector2, vector<double> vector3, vector<double> vector4, const std::string &filepath) {
    ofstream outputFile(filepath);

    for (int i = 0; i < vector1.size(); i++) {
        outputFile << vector1[i] << "\t" << vector2[i] << "\t" << vector3[i] << "\t" << vector4[i] << "\n";
    }

    outputFile.close();
    std::cout << "Saved output." << std::endl;
}



int main() {
    cout << "Reading File..." << endl;
    vector<double> x, y, z, m;
    tie(x, y, z, m) = readDataFromFile("../Bonus_DM/GalaxyFromIllustrisTNG50Dark_DM_Subhalo852966.txt", " ");

    cout << "Calculating..." << endl;
    vector<double> r;
    r = generate_linspace(40, 0., 300.);

    vector<double> masses = mass_shells(x, y, z, m, r);
    vector<double> volumes = volume_shells(r);
    vector<double> densities(r.size());

    #pragma omp parallel for
    for (int i = 0; i < r.size(); ++i) {
        densities[i] = masses[i]/volumes[i];
    }

    cout << "Saving Output..." << endl;

    save_file(r, masses, volumes, densities, "DM_output.txt");

    return 0;
}