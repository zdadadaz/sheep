#include <omp.h>
#include <stdio.h>

int main(){
#pragma omp parallel
{
	int ID = omp_get_thread_num();
	printf("hello %d", ID);
	printf("out %d\n", omp_get_num_threads());
}
}
