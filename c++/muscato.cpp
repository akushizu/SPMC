#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <chrono>
#include "sp.h"
#include "rng.h"
#include "init.h"
#include "save.h"
#include "wmc_muscato.h"
#include "annihilation.h"

using namespace std;

namespace rrr
{
	mt19937_64 gen;
}

uniform_real_distribution<double> rnd(0, 1);

const double pi = 3.14159265358979323846;

const int Nmax = 2000000;

double hbar = 6.62607015 / 1.602176634 / 2 / pi;
double m = 0.32;
double e = 1;

const int Nx = 200;
const int Np = 400;
const int Nt = 500;

vector< double > xrange{ -30, 30 };
vector< double > prange{ -10 * hbar, 10 * hbar };

int initialParticles = 160000;

//int fanni = 10;

double dx = (xrange[1] - xrange[0]) / Nx;
double dp = (prange[1] - prange[0]) / Np;
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
double vwh = 2 * a * sigma * sqrt(2 * pi) / pi / pow(hbar, 2);
double Gamma = vwh * cutoff;

int main()
{
	gen.seed(12345);
	
	vector< Sp* > sp;

	vector< double > x(Nx);
	vector< double > p(Np);
	vector< double > t(Nt + 1);

	for (int i = 0; i < Nx; i++) {

		x[i] = xrange[0] + dx * ((double)i + 1. / 2.);
	}

	for (int i = 0; i < Np; i++) {

		p[i] = prange[0] + dp * ((double)i + 1. / 2.);
	}

	for (int i = 0; i <= Nt; i++) {

		t[i] = 0. + (double)i * dt;
	}

	cout << "Maximum number of particles allowed = " << Nmax << "\n";

	auto t0 = chrono::high_resolution_clock::now();

	int numberOfParticles = init(sp, Nmax, m, x, p, wFunction, initialParticles);

	save(sp, t, 0, x, p);

	auto t1 = t0;
	auto t2 = chrono::high_resolution_clock::now();
	std::chrono::duration< double > ta = t2 - t1;
	std::chrono::duration< double > tb = t2 - t0;
	cout << "Elapsed time is " << ta.count() << " seconds." << "\n";
	cout << "Total elapsed time is " << tb.count() << " seconds." << endl;

	double particleTime = 0;

	for (int i = 1; i <= Nt; i++) {

		particleTime += dt;

		cout << "\nStep " << i << " of " << Nt << " --- time = " << particleTime << "\n\n";

		numberOfParticles = wmc(sp, Nmax, x, p, dt, wKernel, cutoff, vwh, Gamma);

		if (numberOfParticles > 480000) {

			annihilation_a(sp, x, p);
		}

		save(sp, t, i, x, p);

		t1 = t2;
		t2 = chrono::high_resolution_clock::now();
		ta = t2 - t1;
		tb = t2 - t0;
		cout << "Elapsed time is " << ta.count() << " seconds." << "\n";
		cout << "Total elapsed time is " << tb.count() << " seconds." << endl;
	}

	cout << "\nOutput files saved\n" << endl;
}
