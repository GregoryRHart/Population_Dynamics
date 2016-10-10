#include "d_functions.h"

void 
acc_sub(double* part_sum, int start, int end, double sub)
{
	for(long i_s=start; i_s<end; i_s++)
	{
		part_sum[i_s] -= sub;
	}	
}

