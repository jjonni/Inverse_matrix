#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

unsigned long long rdtsc();
double average(double *arr, const int size);
void printM(double **m, int size);
void clearM(double **m, int size);
double** minor_matrix(double **m, int size, int row, int col);
double minor(double **m, int size, int row, int col);
double detM(double **m, int n);
double** transpose_matrix(double **m, int size);
double** attach_matrix(double **m, int size);
double** inverse_matrix(double **m, int size);


int main() {
	omp_set_dynamic(0);
	omp_set_num_threads(10);
	const unsigned long long cpu_HZ = 2400000000ULL;
	unsigned long long start, end;
	int time_repeat = 1;
    double time[time_repeat];
	int m_size = 11;

	double **matrix = (double**)calloc(m_size, sizeof(*matrix)), **inverse_m;
	for (int i = 0; i < m_size; ++i) matrix[i] = (double*)calloc(m_size, sizeof(*matrix[i]));

	double num = -100;
	for (int i = 0; i < m_size; ++i) {
		for (int j = 0; j < m_size; ++j) {
			matrix[i][j] = num;
			if ((i % 2) != 0 && (j % 5) == 0) num /=4;
			if((i % 2) == 0 && (j%3) == 0) num *=2;
			else num += 17;
		}
	} matrix[3][2]=0; matrix[4][0]=0;

	//printf("Matrix:\n\n");
	//printM(matrix, m_size);
	//printf("\n\n");

	for (int i = 0; i < time_repeat; ++i) {
		start = rdtsc();
		// Рассматриваемая операция
		inverse_m = inverse_matrix(matrix, m_size);
		end = rdtsc();
		time[i] = (double)(end - start)/cpu_HZ;
		printf("time = %f\n", time[i]);
	}


	//printf("Inverse matrix:\n\n");
	//printM(inverse_m, m_size);

	clearM(matrix, m_size);
	clearM(inverse_m, m_size);

    printf("time = %f sec\n", average(time, time_repeat));

	return 0;
}

unsigned long long rdtsc() {
	unsigned int lo, hi;
    asm volatile ( "rdtsc\n" : "=a" (lo), "=d" (hi) );

    return ((unsigned long long)hi << 32) | lo;
}

double average(double *arr, const int size) {
	if (size <= 0) {
		fprintf(stderr, "average(): size <= 0");
	}
   	double Sum = 0;


   	for (int i = 0; i < size; ++i) {
   		Sum += arr[i];
   	}
   	return Sum/size;
}

void printM(double **m, int size) {
	if (size == 0) printf("printM(): size of M is 0");

	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			printf("%f ", m[i][j]);
		}
		printf("\n");
	}
}

void clearM(double **m, int size) {
	for (int i = 0; i < size; ++i) {
		free(m[i]);
	}
	free(m);
}

double** minor_matrix(double **m, int size, int row, int col) {
	if (row < 0 || col < 0) {
		fprintf(stderr, "minor(): row < 0 or col < 0\n");
		exit(-1);
	}
	else if (row >= size || col >= size) {
		fprintf(stderr, "minor(): row >= size or col >= size\n");
		exit(-1);
	}


	double **new_m = (double**)calloc(size - 1, sizeof(*new_m));
	for (int i = 0; i < size - 1; ++i) {
		new_m[i] = (double*)calloc(size - 1, sizeof(*new_m[i]));
		if (new_m[i] == NULL) {
			fprintf(stderr, "Allocation failed\n");
			exit(-1);
		}
	}

	int offsetrow = 0, offsetcol = 0;

	int i;
	#pragma omp parallel for collapse(2) shared(offsetrow, offsetcol)
	for (i = 0; i < size; ++i) {
		if (i == row) {
			offsetrow = 1;
			continue;
		}

		offsetcol = 0;

		for (int j = 0; j < size; ++j) {
			if (j == col) {
				offsetcol = 1;
				continue;
			}

			new_m[i - offsetrow][j - offsetcol] = m[i][j];
		}
	}
	//printM(new_m, size_m); printf("\n");
	return new_m;
}

double minor(double **m, int size, int row, int col) {
	double **minor_m = minor_matrix(m, size, row, col);
	double det = detM(minor_m, size - 1);
	clearM(minor_m, size - 1);

	return det;
}

double detM(double **m, int size) {
	if (size == 0) {
		perror("detM(): n == 0");
		exit(-1);
	}
	if (size == 1) {
		return m[0][0];
	}
	else if (size == 2) {
		return m[0][0]*m[1][1] - m[0][1]*m[1][0];
	}

	double det = 0;
	int sign = 1;

	for (int i = 0; i < size; ++i) {
		det += sign*m[0][i] * minor(m, size, 0, i);
		//printf("det = %f\n", det);
		sign = - sign;
	}

	return det;
}

double** transpose_matrix(double **m, int size) {
	double **transposed_m = (double**)calloc(size, sizeof(*transposed_m));
	for (int i = 0; i < size; ++i) transposed_m[i] = (double*)calloc(size, sizeof(transposed_m[i]));

	//#pragma omp parallel for
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			transposed_m[j][i] = m[i][j];
		}
	}
	return transposed_m;
}

double** attach_matrix(double **m, int size) {
	int sign = 1;
	double **attach_m = (double**)calloc(size, sizeof(*attach_m));
	for (int i = 0; i < size; ++i)
		attach_m[i] = (double*)calloc(size, sizeof(attach_m[i]));

	//#pragma parallel for
	for (int i = 0; i < size; ++i) {
		for (int j = 0; j < size; ++j) {
			if ((i+j) % 2 == 0) sign = 1;
			else sign = -1;

			attach_m[i][j] = sign*minor(m, size, i, j);
		}
	}
	//printM(attach_m, size); printf("\n");

	return attach_m;
}

double** inverse_matrix(double **m, int size) {
	double determinant = detM(m, size);
	if (determinant == 0) {
		fprintf(stderr, "Inverse matrix does not exist: detM = 0.\n");
		return NULL;
	} else {
		double **inverse_m = attach_matrix(m, size), **transpose_m;
		//printM(inverse_m, size); printf("\n");

		//#pragma omp parallel for collapse(2) shared(inverse_m)
		for (int i = 0; i < size; ++i) {
			for (int j = 0; j < size; ++j) {
				//printf("Thread %d is processing element inverse_m[%d][%d]\n",
				//		omp_get_thread_num(), i, j);
				inverse_m[i][j] /= determinant;
			}
		}

		transpose_m = transpose_matrix(inverse_m, size);
		clearM(inverse_m, size);
		return transpose_m;
	}
}