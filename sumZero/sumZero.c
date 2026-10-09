#include <stdio.h>
#include <omp.h>

int main(){
    int NT, sum=0, i, id , A;
    omp_lock_t sem;
    omp_init_lock(&sem);

    NT=4;
    omp_set_num_threads(NT);

    #pragma omp parallel private(i, id, A)
    {
        id=omp_get_thread_num();
        A = 1 - 2*(id%2);

        /*
        È possibile utilizzare questa variante al posto dei lock
        
            #pragma omp critical
            {
                sum = sum + A;
            }
        */
        

        omp_set_lock(&sem);
        for (i=0; i<10000; i++){
            sum = sum + A;
        }
        omp_unset_lock(&sem);

        printf("il thread %d ha fatto somma: %d\n", id, sum);
    }

    printf("La somma finale e\': %d", sum);
    return 0;
}
