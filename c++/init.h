#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include "sp.h"
#include "rng.h"

using namespace std;
using namespace rrr;

int init(vector< Sp* >& sp, int Nmax, double m, vector< double >& x, vector< double >& p, function< double(double, double) > fw, int number) {

	cout << "Initializing fw..." << "\n";

	int Nx = (int)x.size();
	int Np = (int)p.size();
	double dx = x[1] - x[0];
	double dp = p[1] - p[0];

	int Fw;
	int actualNumber = 0;

	int count = 0;
	uniform_real_distribution<double> rnd(0, 1);
	
	for (int i = 0; i < Nx; i++) {

		for (int j = 0; j < Np; j++) {

			Fw = (int)round((double)number * fw(x[i], p[j]) * dx * dp);

			actualNumber += (int)fabs(Fw);

			if (actualNumber > Nmax) {

				cout << "The number of particles has reached the limit." << endl;
				exit(0);
			}

			for (int k = 0; k < fabs(Fw); k++) {

				sp.push_back(nullptr);
				sp[count] = new Sp(count);
				
				sp[count]->setm(m);
				sp[count]->setXcell(i);
				sp[count]->setx(x[i] + (rnd(gen) - 0.5) * dx);
				sp[count]->setPcell(j);
				sp[count]->setp(p[j] + (rnd(gen) - 0.5) * dp);
				sp[count]->setWeight(Fw > 0. ? 1 : -1);
				sp[count]->setTime(0);

				sp[count]->setStatus(true);

				count++;
			}
		}
	}

	if (count == 0) {

		cout << "There are no particles in the system." << endl;
		exit(0);
	}
	
	cout << "Initial number of signed particles = " << count << "\n";

	return count;
}