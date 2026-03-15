#include <stdio.h>
#include "sp.h"

void saveRho(struct Sp* sp[], const int number, int rho[], const double t[], const int Nt, int ti, const double x[], const int Nx, int Np) {

	int i, j;

	for (int n = 0; n < number; n++) {

		i = sp[n]->xcell;
		j = sp[n]->pcell;

		if ((i >= 0 && i < Nx) && (j >= 0 && j < Np)) {

			rho[i] += sp[n]->weight;
		}
	}

	FILE *den = NULL;

	if (ti == 0) {

		errno_t err = fopen_s(&den, "rho.dat", "w");
		fprintf(den, "x\trho\n");
		fclose(den);
	}

	errno_t err = fopen_s(&den, "rho.dat", "a");

	for (i = 0; i < 20; i++) {

		fprintf(den, "-");
	}

	fprintf(den, "\n");

	fprintf(den, "# = %d of %d\n", ti, Nt - 1);
	fprintf(den, "t = %g\n", t[ti]);
	fprintf(den, "\n");

	for (i = 0; i < Nx; i++) {

		fprintf(den, "%.3f\t%d\n", x[i], rho[i]);

	}
	
	fprintf(den, "\n");
	fflush(den);
	fclose(den);

	return;
}

void saveFw(struct Sp* sp[], const int number, int* Fw[], const double t[], const int Nt, int ti, const double x[], const int Nx, const double p[], int Np) {

	int i, j;

	for (int n = 0; n < number; n++) {

		i = sp[n]->xcell;
		j = sp[n]->pcell;

		if ((i >= 0 && i < Nx) && (j >= 0 && j < Np)) {

			Fw[i][j] += sp[n]->weight;
		}
	}

	FILE* fw = NULL;
	
	if (ti == 0) {

		errno_t err = fopen_s(&fw, "Fw.dat", "w");
		fprintf(fw, "x\n");

		for (i = 0; i < Nx; i++) {

			fprintf(fw, "%.3f\t", x[i]);
		}

		fprintf(fw, "\n");
		fprintf(fw, "p\n");

		for (i = 0; i < Np; i++) {

			fprintf(fw, "%.4f\t", p[i]);
		}

		fprintf(fw, "\n");
		fprintf(fw, "Fw\n");
		fclose(fw);
	}

	errno_t err = fopen_s(&fw, "Fw.dat", "a");
	
	for (i = 0; i < 20; i++) {

		fprintf(fw, "-");
	}

	fprintf(fw, "\n");

	fprintf(fw, "# = %d of %d\n", ti, Nt - 1);
	fprintf(fw, "t = %g\n", t[ti]);
	fprintf(fw, "\n");

	for (j = Np - 1; j > -1; j--) {

		for (i = 0; i < Nx; i++) {

			fprintf(fw, "%d\t", Fw[i][j]);
		}

		fprintf(fw, "\n");
	}

	fprintf(fw, "\n\n");
	fflush(fw);
	fclose(fw);

	return;
}

void save(struct Sp* sp[], const int number, const double t[], const int Nt, int ti, const double x[], const int Nx, const double p[], const int Np) {

	int* rho = (int*)calloc(Nx, sizeof(int));
	int** Fw = (int**)calloc(Nx, sizeof(int*));

	for (int i = 0; i < Nx; i++) {

		Fw[i] = (int*)calloc(Np, sizeof(int));
	}

	saveRho(sp, number, rho, t, Nt, ti, x, Nx, Np);
	saveFw(sp, number, Fw, t, Nt, ti, x, Nx, p, Np);

	free(rho);
	free(Fw);

	return;
}