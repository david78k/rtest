#include <stdio.h>
#include <omp.h>

int main(int argc, char *argv[]) {
    const int N = 10;
    int i, j;
    double a[N][N], sum = 1;
 
    #pragma omp parallel 
	{
	printf("id: %d\n", omp_get_thread_num());
	printf("cores: %d\n", omp_get_num_threads());
	}

    #pragma omp parallel for
    for (i = 0; i < N; i++) {
	printf("id: %d\n", omp_get_thread_num());
	for(j = 0; j < N; j++) {
        	a[i][j] = (i + 1) * (j + 1);
		//printf("a[i] %d\n", a[i]);
		sum += (a[i][j] / 1000000.0);
	}
    }
 
    printf("sum = %f\n", sum);

    return 0;
}
