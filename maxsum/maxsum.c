#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>


double maxsum(size_t N, size_t ld, double *matrix, int NT){

    size_t subsector = N / NT;
    double max = 0.0;

    omp_set_num_threads(NT);
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        double current_max = 0.0;
        size_t row = 0;
        double sum = 0.0;

        for(size_t i=0; i<subsector; i++){
            
            row = id * subsector + i;
            sum = 0.0;
            
            for(size_t j=0; j<N; j++){
                sum += sqrt(matrix[row * ld + j]);
            }

            if(sum > current_max){
                current_max = sum;
            }
        }

        #pragma omp critical 
        {
            if(current_max > max){
                max = current_max;
            }
        }
    }

    return max;
}