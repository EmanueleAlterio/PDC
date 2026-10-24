//Emanuele Alterio 

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>


double maxsum(int N, int ld, double *matrix, int NT){

    double max = 0.0;

    #pragma omp parallel num_threads(NT)
    {
        int id = omp_get_thread_num();
        double current_max = 0.0;

        #pragma omp for nowait
        for(int row=0; row<N; row++){
            
            double sum = 0.0;
            
            for(int j=0; j<N; j++){
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