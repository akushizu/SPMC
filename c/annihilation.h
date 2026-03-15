#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "sp.h"
#include "rng.h"

void density(struct Sp* sp[], const int number, const int Nx, const int Np, int* Fw[]) {

	int i, j;

	for (int n = 0; n < number; n++) {

		i = sp[n]->xcell;
		j = sp[n]->pcell;

		if ((i >= 0 && i < Nx) && (j >= 0 && j < Np)) {

			Fw[i][j] += sp[n]->weight;
		}
	}

	return;
}

int comp(const void* sp1p, const void* sp2p) {

	if (((*(struct Sp**)sp1p)->status == true) && ((*(struct Sp**)sp2p)->status == false)) {

		return -1;
	}	
	else if (((*(struct Sp**)sp1p)->status == false) && ((*(struct Sp**)sp2p)->status == false)) {

		return 0;
	}
	else if (((*(struct Sp**)sp1p)->status == true) && ((*(struct Sp**)sp2p)->status == true)) {

		if ((*(struct Sp**)sp1p)->xcell < (*(struct Sp**)sp2p)->xcell) {

			return -1;
		}
		else if ((*(struct Sp**)sp1p)->xcell > (*(struct Sp**)sp2p)->xcell) {

			return 1;
		}
		else if ((*(struct Sp**)sp1p)->pcell < (*(struct Sp**)sp2p)->pcell) {

			return -1;
		}
		else if ((*(struct Sp**)sp1p)->pcell > (*(struct Sp**)sp2p)->pcell) {

			return 1;
		}
		else {

			return 0;
		}
	}
	else {

		return 1;
	}
}

void check_range(struct Sp* sp[], const int number, int Nx, int Np) {

	for (int n = 0; n < number; n++) {

		if (sp[n]->status == false) {

			continue;
		}
		else {

			if (!(sp[n]->xcell >= 0 && sp[n]->xcell < Nx) || !(sp[n]->pcell >= 0 && sp[n]->pcell < Np)) {

				sp[n]->status = false;
			}
		}
	}
}

int sparse_matrix(struct Sp* sp[], const int number, int* sparse[]) {

	int count = 0;

	if (sp[0]->status != false) {

		sparse[count] = (int*)calloc(3, sizeof(int));
		sparse[count][0] = sp[0]->xcell;
		sparse[count][1] = sp[0]->pcell;
		sparse[count][2] = sp[0]->weight;
		count++;
	}

	for (unsigned int i = 1; i < number; i++) {

		if (sp[i]->status == false) {

			continue;
		}
		else if (sp[i]->xcell != sp[i - 1]->xcell) {

			if (sparse[count - 1][2] == 0) {

				free(sparse[count - 1]);
				sparse[count - 1] = NULL;
				count--;
			}

			sparse[count] = (int*)calloc(3, sizeof(int));
			sparse[count][0] = sp[i]->xcell;
			sparse[count][1] = sp[i]->pcell;
			sparse[count][2] = sp[i]->weight;
			count++;
		}
		else if (sp[i]->pcell != sp[i - 1]->pcell) {

			if (sparse[count - 1][2] == 0) {

				free(sparse[count - 1]);
				sparse[count - 1] = NULL;
				count--;
			}

			sparse[count] = (int*)calloc(3, sizeof(int));
			sparse[count][0] = sp[i]->xcell;
			sparse[count][1] = sp[i]->pcell;
			sparse[count][2] = sp[i]->weight;
			count++;
		}
		else {

			sparse[count - 1][2] += sp[i]->weight;
		}
	}

	return count;
}

int annihilation(struct Sp* sp[], const int number, const double x[], const int Nx, const double p[], const double Np) {
	// grid, O(n)

	printf("\nAnnihilating particles...\n");
	
	double dx = x[1] - x[0];
	double dp = p[1] - p[0];

	printf("# of particles before annihilation = %d\n", number);

	int** Fw = (int**)calloc(Nx, sizeof(int*));

	for (int i = 0; i < Nx; i++) {

		Fw[i] = (int*)calloc(Np, sizeof(int));
	}
	
	density(sp, number, Nx, Np, Fw);

	int count = 0;

	for (int i = 0; i < Nx; i++) {

		for (int j = 0; j < Np; j++) {

			for (int k = 0; k < abs(Fw[i][j]); k++) {

				sp[count]->number = count;

				sp[count]->xcell = i;
				sp[count]->x = x[i] + (rnd() - 0.5) * dx;
				sp[count]->pcell = j;
				sp[count]->p = p[j] + (rnd() - 0.5) * dp;
				sp[count]->weight = Fw[i][j] > 0. ? 1 : -1;

				sp[count]->status = true;

				count++;
			}
		}
	}
					
	for (int i = count; i < number; i++) {

		free(sp[i]);
		sp[i] = NULL;
	}

	free(Fw);

	printf("# of particles after annihilation = %d\n\n", count);

	return count;
}

int annihilation_a(struct Sp* sp[], const int number, const double x[], const int Nx, const double p[], const double Np) {
	// sort, O(n log n)

	printf("\nAnnihilating particles...\n");

	double dx = x[1] - x[0];
	double dp = p[1] - p[0];

	printf("# of particles before annihilation = %d\n", number);

	check_range(sp, number, Nx, Np);

	qsort(sp, number, sizeof(sp[0]), comp);

	int** sparse = (int**)calloc(number, sizeof(int*));

	int sparnum = sparse_matrix(sp, number, sparse);

	int count = 0;

	double m = sp[0]->m;
	double t = sp[0]->time;

	for (int i = 0; i < number; i++) {

		free(sp[i]);
		sp[i] = NULL;
	}

	for (int i = 0; i < sparnum; i++) {

		for (int n = 0; n < abs(sparse[i][2]); n++) {

			sp[count] = (struct Sp*)malloc(sizeof(struct Sp));;
			sp[count]->number = count;

			sp[count]->m = m;
			sp[count]->xcell = sparse[i][0];
			sp[count]->x = x[sparse[i][0]] + (rnd() - 0.5) * dx;
			sp[count]->pcell = sparse[i][1];
			sp[count]->p = p[sparse[i][1]] + (rnd() - 0.5) * dp;
			sp[count]->weight = sparse[i][2] > 0. ? 1 : -1;
			sp[count]->time = t;

			sp[count]->status = true;

			count++;
		}
	}

	for (int i = count; i < number; i++) {

		free(sp[i]);
		sp[i] = NULL;
	}

	free(sparse);

	printf("# of particles after annihilation = %d\n\n", count);

	return count;
}