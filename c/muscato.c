#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <math.h>
#include "sp.h"
#include "init.h"
#include "save.h"
#include "wmc_muscato.h"
#include "annihilation.h"

#define pi 3.14159265358979323846

#define Nmax 2000000

#define hbar 6.62607015 / 1.602176634 / 2 / pi
double m = 0.32;
double e = 1;

#define Nx 200
#define Np 400
#define Nt 500

double xrange[] = { -30, 30 };
double prange[] = { -10 * hbar, 10 * hbar };

int initialParticles = 160000;

//int fanni = 10;

#define dx (xrange[1] - xrange[0]) / Nx
#define dp (prange[1] - prange[0]) / Np
double dt = 0.05;

double x0 = -15;
double s0 = 2.852;
double p0 = 0.7 * hbar;
double xc = 0;
double a = 0.3;
double sigma = 1;

static double wFunction(double x, double p) {

	return 1. / pi / hbar * exp(-0.5 * pow(((x - x0) / s0), 2)) * exp(-2. * pow((s0 * (p - p0) / hbar), 2));
}

static double potential(double x) {

	return e * a * exp(-0.5 * pow(((x - xc) / sigma), 2));
}

static double wKernel(double x, double p) {

	return 2. * a * sigma * sqrt(2 * pi) / pi / pow(hbar, 2) * exp(-2. * pow((sigma * p / hbar), 2)) * sin(2. * p * x / hbar);
}

double cutoff = 8 * hbar;
#define vwh 2 * a * sigma * sqrt(2 * pi) / pi / pow(hbar, 2)
#define Gamma vwh * cutoff

int main()
{
	srand(12345);
	
	struct Sp* sp[Nmax] = {NULL};

	double x[Nx];
	double p[Np];
	double t[Nt + 1];

	for (int i = 0; i < Nx; i++) {

		x[i] = xrange[0] + dx * ((double)i + 1. / 2.);
	}

	for (int i = 0; i < Np; i++) {

		p[i] = prange[0] + dp * ((double)i + 1. / 2.);
	}

	for (int i = 0; i <= Nt; i++) {

		t[i] = 0. + (double)i * dt;
	}

	printf("Maximum number of particles allowed = %d\n\n", Nmax);

	struct timespec t0;
	timespec_get(&t0, TIME_UTC);

	int numberOfParticles = init(sp, Nmax, m, x, Nx, p, Np, wFunction, initialParticles);

	save(sp, numberOfParticles, t, Nt, 0, x, Nx, p, Np);

	struct timespec t1 = t0;
	
	struct timespec t2;
	timespec_get(&t2, TIME_UTC);

	double tta = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) / 1000000000.0;
	double ttb = (double)(t2.tv_sec - t0.tv_sec) + (double)(t2.tv_nsec - t0.tv_nsec) / 1000000000.0;

	printf("Elapsed time is %f seconds.\n", tta);
	printf("Total elapsed time is %f seconds.\n", ttb);

	double particleTime = 0;

	for (int i = 1; i <= Nt; i++) {

		particleTime += dt;

		printf("\nStep %d of %d --- time = %g\n\n", i, Nt, particleTime);

		numberOfParticles = wmc(sp, numberOfParticles, Nmax, x, Nx, p, Np, dt, wKernel, cutoff, vwh, Gamma);

		if (numberOfParticles > 480000) {

			numberOfParticles = annihilation(sp, numberOfParticles, x, Nx, p, Np);
		}

		save(sp, numberOfParticles, t, Nt, i, x, Nx, p, Np);

		t1 = t2;
		timespec_get(&t2, TIME_UTC);
		tta = (double)(t2.tv_sec - t1.tv_sec) + (double)(t2.tv_nsec - t1.tv_nsec) / 1000000000.0;
		ttb = (double)(t2.tv_sec - t0.tv_sec) + (double)(t2.tv_nsec - t0.tv_nsec) / 1000000000.0;

		printf("Elapsed time is %f seconds.\n", tta);
		printf("Total elapsed time is %f seconds.\n", ttb);
	}

	printf("\nOutput files saved\n\n");
}
