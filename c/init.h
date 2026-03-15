#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "sp.h"
#include "rng.h"

int init(struct Sp* sp[], int Nmax, double m, const double x[], const int Nx, const double p[], const int Np, double (*fw)(double, double), int number) {

	printf("Initializing fw...\n");

	double dx = x[1] - x[0];
	double dp = p[1] - p[0];

	int Fw;
	int actualNumber = 0;

	int count = 0;
	
	for (int i = 0; i < Nx; i++) {

		for (int j = 0; j < Np; j++) {

			Fw = (int)round((double)number * fw(x[i], p[j]) * dx * dp);

			actualNumber += abs(Fw);

			if (actualNumber > Nmax) {

				printf("The number of particles has reached the limit.");
				exit(0);
			}

			for (int k = 0; k < abs(Fw); k++) {

				sp[count] = (struct Sp*)malloc(sizeof(struct Sp));
				sp[count]->number = count;
				
				sp[count]->m = m;
				sp[count]->xcell = i;
				sp[count]->x = x[i] + (rnd() - 0.5) * dx;
				sp[count]->pcell = j;
				sp[count]->p = p[j] + (rnd() - 0.5) * dp;
				sp[count]->weight = Fw > 0. ? 1 : -1;
				sp[count]->time = 0;

				sp[count]->status = true;

				count++;
			}
		}
	}

	if (count == 0) {

		printf("There are no particles in the system.\n");
		exit(0);
	}
	
	printf("Initial number of signed particles = %d\n", count);

	return count;
}