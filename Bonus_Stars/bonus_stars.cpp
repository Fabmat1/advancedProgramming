#include <iostream>
#include <random>
#include <chrono>
#include <immintrin.h>
#include <omp.h>
#include <fstream>
#include <sstream>
#include <tuple>
#include <iomanip>
#include <Eigen/Eigenvalues>

using namespace std;


tuple<vector<double>, vector<double>, vector<double>, vector<double>, vector<double>, vector<double>> readDataFromFile(const string& filename, const string& separator) {

    vector<double> col1, col2, col3, col4, col5, col6;
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

        if (tokens.size() == 6) {
            col1.push_back(stod(tokens[0]));
            col2.push_back(stod(tokens[1]));
            col3.push_back(stod(tokens[2]));
            col4.push_back(stod(tokens[3]));
            col5.push_back(stod(tokens[4]));
            col6.push_back(stod(tokens[5]));
        }
    }

    return make_tuple(col1, col2, col3, col4, col5, col6);
}


double diag_moi(vector<double> vec1, vector<double> vec2, vector<double> m){
    double result;

    #pragma omp parallel for
    for (int i = 0; i < vec1.size(); ++i) {
        #pragma omp atomic
        result += m[i]*(vec1[i]*vec1[i]+vec2[i]*vec2[i]);
    }

   return result;
}


double nondiag_moi(vector<double> vec1, vector<double> vec2, vector<double> m){
    double result;

    #pragma omp parallel for
    for (int i = 0; i < vec1.size(); ++i) {
    #pragma omp atomic
        result -= m[i]*(vec1[i]*vec2[i]);
    }

    return result;
}


vector<vector<double>> calculateEigen(const std::vector<std::vector<double>>& matrix) {
    // Convert the matrix to an Eigen MatrixXd object
    Eigen::Matrix3d eigenMatrix;
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            eigenMatrix(i, j) = matrix[i][j];
        }
    }

    // Compute the eigenvalues and eigenvectors
    Eigen::EigenSolver<Eigen::Matrix3d> solver(eigenMatrix);
    Eigen::Matrix3cd eigenvectors = solver.eigenvectors();

    // Convert the eigenvectors to vector<vector<double>> format
    vector<vector<double>> result;
    for (int i = 0; i < 3; ++i) {
        vector<double> eigenvector;
        for (int j = 0; j < 3; ++j) {
            eigenvector.push_back(eigenvectors(j, i).real());
        }
        result.push_back(eigenvector);
    }

    return result;
}


vector<double> crossProduct(const std::vector<double>& a, const vector<double>& b) {
    vector<double> result(3);
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = a[2] * b[0] - a[0] * b[2];
    result[2] = a[0] * b[1] - a[1] * b[0];
    return result;
}


vector<vector<double>> rotate_coords(vector<double> x, vector<double> y, vector<double> z, vector<double> m){
    vector<vector<double>> moi_tensor(3, vector<double>(3));

    moi_tensor[0][0] = diag_moi(y, z, m);
    moi_tensor[1][1] = diag_moi(x, z, m);
    moi_tensor[2][2] = diag_moi(x, y, m);

    moi_tensor[0][1] = nondiag_moi(x, y, m);
    moi_tensor[1][0] = moi_tensor[0][1];
    moi_tensor[0][2] = nondiag_moi(x, z, m);
    moi_tensor[2][0] = moi_tensor[0][2];
    moi_tensor[1][2] = nondiag_moi(y, z, m);
    moi_tensor[2][1] = moi_tensor[1][2];

    vector<vector<double>> eigenvector;
    eigenvector = calculateEigen(moi_tensor);

    vector<double> xdir(3);
    vector<double> ydir(3);
    vector<double> zdir(3);

    xdir = eigenvector[0];
    ydir = eigenvector[1];
    zdir = crossProduct(xdir, ydir);

    // For some reason these have flipped signs when compared to the python implementation
    for (int i = 0; i < 3; ++i) {
        xdir[i] = -xdir[i];
        ydir[i] = -ydir[i];
    }

    vector<double> x_prime(x.size());
    vector<double> y_prime(x.size());
    vector<double> z_prime(x.size());

    #pragma omp parallel for
    for (int i = 0; i < x.size(); ++i) {
        x_prime[i] = xdir[0]*x[i]+xdir[1]*y[i]+xdir[2]*z[i];
        y_prime[i] = ydir[0]*x[i]+ydir[1]*y[i]+ydir[2]*z[i];
        z_prime[i] = zdir[0]*x[i]+zdir[1]*y[i]+zdir[2]*z[i];
    }

    vector<vector<double>> result(3, vector<double>(x.size()));
    result[0] = x_prime;
    result[1] = y_prime;
    result[2] = z_prime;

    return result;
}


void save_file(vector<double> vector1, vector<double> vector2, vector<double> vector3, vector<double> vector4, vector<double> vector5, vector<double> vector6, const std::string &filepath) {
    ofstream outputFile(filepath);

    for (int i = 0; i < vector1.size(); i++) {
        outputFile << fixed << setprecision(8) << vector1[i] << "\t" << vector2[i] << "\t" << vector3[i] << "\t" << vector4[i] << "\t" << vector5[i] << "\t" << vector6[i] << "\n";
    }

    outputFile.close();
    std::cout << "Saved output." << std::endl;
}



int main() {
    vector<double> x, y, z, m, m_init, z_formed;
    cout << "Reading data..." << endl;
    tie(x, y, z, m, m_init, z_formed) = readDataFromFile("../Bonus_Stars/GalaxyFromIllustrisTNG50_Stars_Subhalo521803.txt", " ");
    cout << "Read data successfully!" << endl;

    cout << "Calculating rotated coordinates..." << endl;
    vector<vector<double>> rotated_coords(3, vector<double>(x.size()));
    rotated_coords = rotate_coords(x, y, z, m);
    cout << "Successfully transformed coordinate system!" << endl;

    x = rotated_coords[0];
    y = rotated_coords[1];
    z = rotated_coords[2];

    cout << "Saving transformed coordinate system..." << endl;

    save_file(x, y, z, m, m_init, z_formed, "rotated_stars.csv");

    cout << "All done!" << endl;

    return 0;
}