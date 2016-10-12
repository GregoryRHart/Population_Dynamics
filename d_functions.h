#ifndef D_FUNCTIONS
#define D_FUNCTIONS

#define SUCCESS 0xF00D
#define FAILURE 0xDEAD

int initialize_gpu();
int allocate_device_partsum(double** d_part_sum, int size);
int free_device_partsum(double** d_part_sum);

void acc_sub(double* part_sum, int start, int end, double sub);
int cuda_sub(double* d_part_sum, double* part_sum, int start, int end, double sub);

#endif
