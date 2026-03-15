#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "sp.h"
#include "rng.h"

void propagate(struct Sp* sp, const double x0, const double dx, const double dt) {

	sp->x = sp->x + sp->p / sp->m * dt;
	sp->xcell = (int)floor((sp->x - (x0 - dx / 2.)) / dx);
	sp->time = sp->time + dt;
}

bool range_check(struct Sp* sp, const int Nx, const int Np) {

	if ((sp->xcell >= 0 && sp->xcell < Nx) && (sp->pcell >= 0 && sp->pcell < Np)) {

		return true;
	}
	else {

		return false;
	}
}

void create_particles(struct Sp* sp[], const int n, int* currentParticles, int* newParticles, int Nmax, double dt, double p0, double dp, double (*Vw)(double, double), double cutoff, double vwh, double Gamma) {

	if (rnd() < (Gamma * dt)) {

		double ps = (2 * rnd() - 1) * cutoff;

		if (rnd() < fabs(Vw(sp[n]->x, ps)) / vwh) {

			*currentParticles += 2;
			*newParticles += 2;

			if (*currentParticles > Nmax) {

				printf("The number of particles has reached the limit.\n");
				exit(0);
			}

			sp[*currentParticles - 2] = (struct Sp*)malloc(sizeof(struct Sp));
			sp[*currentParticles - 2]->number = *currentParticles - 2;
			sp[*currentParticles - 2]->m = sp[n]->m;
			sp[*currentParticles - 2]->x = sp[n]->x;
			sp[*currentParticles - 2]->xcell = sp[n]->xcell;
			sp[*currentParticles - 2]->p = sp[n]->p + ps;
			sp[*currentParticles - 2]->pcell = (int)floor((sp[*currentParticles - 2]->p - (p0 - dp / 2.)) / dp);
			sp[*currentParticles - 2]->weight = sp[n]->weight * (Vw(sp[n]->x, ps) > 0. ? 1 : -1);
			sp[*currentParticles - 2]->time = sp[n]->time;
			sp[*currentParticles - 2]->status = true;

			sp[*currentParticles - 1] = (struct Sp*)malloc(sizeof(struct Sp));
			sp[*currentParticles - 1]->number = *currentParticles - 1;
			sp[*currentParticles - 1]->m = sp[n]->m;
			sp[*currentParticles - 1]->x = sp[n]->x;
			sp[*currentParticles - 1]->xcell = sp[n]->xcell;
			sp[*currentParticles - 1]->p = sp[n]->p - ps;
			sp[*currentParticles - 1]->pcell = (int)floor((sp[*currentParticles - 1]->p - (p0 - dp / 2.)) / dp);
			sp[*currentParticles - 1]->weight = -sp[n]->weight * (Vw(sp[n]->x, ps) > 0. ? 1 : -1);
			sp[*currentParticles - 1]->time = sp[n]->time;
			sp[*currentParticles - 1]->status = true;
		}
	}
}

int wmc(struct Sp* sp[], const int number, int Nmax, const double x[], const int Nx, const double p[], const int Np, double dt, double (*Vw)(double, double), double cutoff, double vwh, double Gamma) {

	printf("Evolving particles...\n");

	double x0 = x[0];
	double dx = x[1] - x[0];
	double p0 = p[0];
	double dp = p[1] - p[0];

	int originalParticles = number;

	int currentParticles = originalParticles;
	int newParticles = 0, outsideParticles = 0;

	for (int n = 0; n < originalParticles; n++) {

		if (sp[n]->status == false) {

			outsideParticles++;
			continue;
		}

		propagate(sp[n], x0, dx, dt);

		if (range_check(sp[n], Nx, Np) == false) {

			sp[n]->status = false;
			outsideParticles++;
			continue;
		}

		create_particles(sp, n, &currentParticles, &newParticles, Nmax, dt, p0, dp, Vw, cutoff, vwh, Gamma);
	}

	printf("Number of created particles = %d\n", newParticles);
	printf("Total number of particles = %d\n", currentParticles);

	for (int n = originalParticles; n < currentParticles; n++) {

		if (range_check(sp[n], Nx, Np) == false) {

			sp[n]->status = false;
			outsideParticles++;
		}
	}

	printf("Number of outside particles = %d\n", outsideParticles);

	return currentParticles;
}
