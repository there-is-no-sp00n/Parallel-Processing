#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

double integrate(double x)
{
	return 4/(1 + x * x);
}


int main(int argc, char *argv[])
{
	int num_of_procs, n; //number of processors and number of intervals
	int my_id;
	int i;

	double given_pi, my_pi, y, x1, x2, h, sum, total;

	double begin_time, end_time;
	
	given_pi = 3.141592653589793238462643;	

	num_of_procs = atoi(argv[1]);
	n = atoi(argv[2]);

	double *num_in_proc = (double *)calloc(n, sizeof(double));

	MPI_Status status;
	
	//initialize MPI
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_procs);

	//start the clock
	begin_time = MPI_Wtime();
	
	//the end points of the rectangles
	x1 = ((double) my_id)/((double)num_of_procs);
	x2 = ((double) (my_id +1))/((double)num_of_procs);

	//printf("x1, x2 is %f %f\n", x1, x2);

	//upper limit - lower limt over number of intervals * number of processors
	//h = 1.0/((double)(n * num_of_procs));
	sum = 0.0;

	for(i = 1; i < n; i++)
	{
		y = x1 + (x2 - x1) * (((double) i)/((double) n));
		sum = (double)integrate(y);
	}

	//getting the endpoints integrated
	sum += integrate(x1);
	sum += integrate(x2);
	sum /= 2;
	total = sum;

	if(my_id == 0)
	{
		for(i = 1; i <= num_of_procs; i++)
		{
			if (i == num_of_procs)
			{
				my_pi = num_in_proc[i-1];
				break;
			}
			
			MPI_Recv(&(num_in_proc[i]), 1, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status);
			
		}
		
	}
	else
	{
		MPI_Send(&total, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
	}
	
	if(my_id == 0)
	{
		end_time = MPI_Wtime();
		//for(i =0; i < num_of_procs; i++)
		//{
			//printf("my pi is %f \n", my_pi);
			//my_pi += num_in_proc[i];
			//printf("iiiin number %d %f\n", i, num_in_proc[i]);
		//}
		//my_pi *= h/3;
		
		printf("\n*****************\n");
		printf("\n\nPI is %f and Error is %f \n", my_pi, fabs(my_pi-given_pi));
		printf("Time taken is %f \n\n", (end_time-begin_time));
        	printf("\n*****************\n");
	}

	MPI_Finalize();
}

