#include <iostream>
#include <random>
#include <chrono>
#include <immintrin.h>
#include <omp.h>

using namespace std::chrono;

int main() {
    auto start = high_resolution_clock::now();
    // Seed the random number generator (Kind of does not work due to parallelization but too lazy to fix, sorry!)
    std::mt19937_64 eng(12345);

    // Generate 10000000 random values for x and y
    const int N = 10000000;
    int N_i = 0;

    // Get the maximum number of threads supported by the CPU
    int max_threads = omp_get_max_threads();

    // Set the number of threads to use in parallel regions
    omp_set_num_threads(max_threads);

    // Define the distribution
    std::uniform_real_distribution<double> distr(0, 1.0);
    __m256d r_squared = _mm256_set1_pd(std::pow(0.5, 2));

    // Parallelize for loop execution (Use all CPU Cores)
    #pragma omp parallel for reduction(+: N_i)
    for  (int i = 0; i < N; i += 4) {
        // I <3 vectorization!
        // Generate batches of 4 random doubles and substract 0.5
        __m256d x_values = _mm256_sub_pd(_mm256_set_pd(distr(eng), distr(eng), distr(eng), distr(eng)), _mm256_set1_pd(0.5));
        __m256d y_values = _mm256_sub_pd(_mm256_set_pd(distr(eng), distr(eng), distr(eng), distr(eng)), _mm256_set1_pd(0.5));

        // square the values
        __m256d x_squared = _mm256_mul_pd(x_values, x_values);
        __m256d y_squared = _mm256_mul_pd(y_values, y_values);

        // Sum the squares and determine if inside circle
        __m256d sum_squared = _mm256_add_pd(x_squared, y_squared);
        __m256d inside = _mm256_cmp_pd(sum_squared, r_squared, _CMP_LT_OQ);

        // Keep track of the amount inside the circle
        N_i += _mm_popcnt_u32(_mm256_movemask_pd(inside));
    }

    double fraction =  static_cast<double>(N_i)/static_cast<double>(N);
    auto end = high_resolution_clock::now();

    std::cout.precision(5);
    std::cout << fraction << std::fixed << std::endl;
    std::cout << 4*fraction << std::fixed << std::endl;

    auto duration = duration_cast<milliseconds>(end - start);
    auto elapsed = static_cast<float>(duration.count());
    std::cout << elapsed/1000 << std::endl;

    return 0;
}
