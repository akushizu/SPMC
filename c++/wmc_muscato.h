#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include "sp.h"

using namespace std;

void propagate(Sp* sp, const double& x0, const double& dx, const double& dt) {

	sp->setx(sp->getx() + sp->getp() / sp->getm() * dt);
	sp->setXcell((int)floor((sp->getx() - (x0 - dx / 2.)) / dx));
	sp->setTime(sp->getTime() + dt);
}

bool range_check(Sp* sp, const int& Nx, const int& Np) {

	if ((sp->getXcell() >= 0 && sp->getXcell() < Nx) && (sp->getPcell() >= 0 && sp->getPcell() < Np)) {

		return true;
	}
	else {

		return false;
	}
}

void create_particles(vector< Sp* >& sp, const int n, int& currentParticles, int& newParticles, int& Nmax, double& dt, double& p0, double& dp, function< double(double, double) >& Vw, double& cutoff, double& vwh, double& Gamma) {

	uniform_real_distribution<double> rnd(0, 1);

	if (rnd(gen) < (Gamma * dt)) {

		double ps = (2 * rnd(gen) - 1) * cutoff;

		if (rnd(gen) < fabs(Vw(sp[n]->getx(), ps)) / vwh) {

			currentParticles += 2;
			newParticles += 2;

			if (currentParticles > Nmax) {

				cout << "The number of particles has reached the limit." << endl;
				exit(0);
			}

			sp.push_back(nullptr);
			sp[currentParticles - 2] = new Sp(currentParticles - 2);
			sp[currentParticles - 2]->setm(sp[n]->getm());
			sp[currentParticles - 2]->setx(sp[n]->getx());
			sp[currentParticles - 2]->setXcell(sp[n]->getXcell());
			sp[currentParticles - 2]->setp(sp[n]->getp() + ps);
			sp[currentParticles - 2]->setPcell((int)floor((sp[currentParticles - 2]->getp() - (p0 - dp / 2.)) / dp));
			sp[currentParticles - 2]->setWeight(sp[n]->getWeight() * (Vw(sp[n]->getx(), ps) > 0. ? 1 : -1));
			sp[currentParticles - 2]->setTime(sp[n]->getTime());
			sp[currentParticles - 2]->setStatus(true);

			sp.push_back(nullptr);
			sp[currentParticles - 1] = new Sp(currentParticles - 1);
			sp[currentParticles - 1]->setm(sp[n]->getm());
			sp[currentParticles - 1]->setx(sp[n]->getx());
			sp[currentParticles - 1]->setXcell(sp[n]->getXcell());
			sp[currentParticles - 1]->setp(sp[n]->getp() - ps);
			sp[currentParticles - 1]->setPcell((int)floor((sp[currentParticles - 1]->getp() - (p0 - dp / 2.)) / dp));
			sp[currentParticles - 1]->setWeight(-sp[n]->getWeight() * (Vw(sp[n]->getx(), ps) > 0. ? 1 : -1));
			sp[currentParticles - 1]->setTime(sp[n]->getTime());
			sp[currentParticles - 1]->setStatus(true);
		}
	}
}

int wmc(vector< Sp* >& sp, int Nmax, vector< double >& x, vector< double >& p, double dt, function< double(double, double) > Vw, double& cutoff, double& vwh, double& Gamma) {

	cout << "Evolving particles..." << "\n";

	int Nx = (int)x.size();
	int Np = (int)p.size();
	double x0 = x[0];
	double dx = x[1] - x[0];
	double p0 = p[0];
	double dp = p[1] - p[0];

	int originalParticles = (int)sp.size();

	int currentParticles = originalParticles;
	int newParticles = 0, outsideParticles = 0;

	for (int n = 0; n < originalParticles; n++) {

		if (sp[n]->getStatus() == false) {

			outsideParticles++;
			continue;
		}

		propagate(sp[n], x0, dx, dt);

		if (range_check(sp[n], Nx, Np) == false) {

			sp[n]->setStatus(false);
			outsideParticles++;
			continue;
		}

		create_particles(sp, n, currentParticles, newParticles, Nmax, dt, p0, dp, Vw, cutoff, vwh, Gamma);
	}

	cout << "Number of created particles = " << newParticles << "\n";
	cout << "Total number of particles = " << currentParticles << "\n";

	for (int n = originalParticles; n < currentParticles; n++) {

		if (range_check(sp[n], Nx, Np) == false) {

			sp[n]->setStatus(false);
			outsideParticles++;
		}
	}

	cout << "Number of outside particles = " << outsideParticles << "\n";

	return currentParticles;
}
