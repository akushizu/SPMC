#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include "sp.h"
#include "rng.h"

using namespace std;

void density(vector< Sp* >& sp, const int Nx, const int Np, vector< vector< int > >& Fw) {

	int i, j;

	for (int n = 0; n < sp.size(); n++) {

		i = sp[n]->getXcell();
		j = sp[n]->getPcell();

		if ((i >= 0 && i < Nx) && (j >= 0 && j < Np)) {

			Fw[i][j] += sp[n]->getWeight();
		}
	}

	return;
}

bool comp(Sp* sp1, Sp* sp2) {

	if ((sp1->getStatus() == true) && (sp2->getStatus() == false)) {

		return true;
	}	
	else if ((sp1->getStatus() == false) && (sp2->getStatus() == false)) {

		return false;
	}
	else if ((sp1->getStatus() == true) && (sp2->getStatus() == true)) {

		if (sp1->getXcell() < sp2->getXcell()) {

			return true;
		}
		else if (sp1->getXcell() > sp2->getXcell()) {

			return false;
		}
		else if (sp1->getPcell() < sp2->getPcell()) {

			return true;
		}
		else if (sp1->getPcell() > sp2->getPcell()) {

			return false;
		}
		else {

			return false;
		}
	}
	else {

		return false;
	}
}

void check_range(vector< Sp* >& sp, int Nx, int Np) {

	for (int n = 0; n < sp.size(); n++) {

		if (sp[n]->getStatus() == false) {

			continue;
		}
		else {

			if (!(sp[n]->getXcell() >= 0 && sp[n]->getXcell() < Nx) || !(sp[n]->getPcell() >= 0 && sp[n]->getPcell() < Np)) {

				sp[n]->setStatus(false);
			}
		}
	}
}

void sparse_matrix(vector< Sp* >& sp, vector< vector < int > >& sparse) {

	if (sp[0]->getStatus() != false) {

		sparse.push_back(vector< int >{sp[0]->getXcell(), sp[0]->getPcell(), sp[0]->getWeight()});
	}

	for (unsigned int i = 1; i < sp.size(); i++) {

		if (sp[i]->getStatus() == false) {

			continue;
		}
		else if (sp[i]->getXcell() != sp[i - 1]->getXcell()) {

			if (sparse.back()[2] == 0) {

				sparse.pop_back();
			}

			sparse.push_back(vector< int >{sp[i]->getXcell(), sp[i]->getPcell(), sp[i]->getWeight()});
		}
		else if (sp[i]->getPcell() != sp[i - 1]->getPcell()) {

			if (sparse.back()[2] == 0) {

				sparse.pop_back();
			}

			sparse.push_back(vector< int >{sp[i]->getXcell(), sp[i]->getPcell(), sp[i]->getWeight()});
		}
		else {

			sparse[sparse.size() - 1][2] += sp[i]->getWeight();
		}
	}
}

void annihilation(vector< Sp* >& sp, vector< double >& x, vector< double >& p) {
	// grid, O(n)

	cout << "\nAnnihilating particles..." << "\n";
	
	int Nx = (int)x.size();
	int Np = (int)p.size();
	double dx = x[1] - x[0];
	double dp = p[1] - p[0];

	cout << "# of particles before annihilation = " << sp.size() << "\n";

	vector< vector< int > > Fw(Nx, vector< int >(Np, 0));
	
	density(sp, Nx, Np, Fw);

	int count = 0;
	uniform_real_distribution<double> rnd(0, 1);

	for (int i = 0; i < Nx; i++) {

		for (int j = 0; j < Np; j++) {

			for (int k = 0; k < fabs(Fw[i][j]); k++) {

				sp[count]->setNumber(count);

				sp[count]->setXcell(i);
				sp[count]->setx(x[i] + (rnd(gen) - 0.5) * dx);
				sp[count]->setPcell(j);
				sp[count]->setp(p[j] + (rnd(gen) - 0.5) * dp);
				sp[count]->setWeight(Fw[i][j] > 0. ? 1 : -1);

				sp[count]->setStatus(true);

				count++;
			}
		}
	}
					
	for (int i = count; i < sp.size(); i++) {

		delete sp[i];
		sp[i] = nullptr;
	}

	while (sp.back() == nullptr) {

		sp.pop_back();
	}

	cout << "# of particles after annihilation = " << count << "\n\n";

	return;
}

void annihilation_a(vector< Sp* >& sp, vector < double >& x, vector < double >& p) {
	// sort, O(n log n)

	cout << "\nAnnihilating particles...\n";

	int Nx = (int)x.size();
	int Np = (int)p.size();
	double dx = x[1] - x[0];
	double dp = p[1] - p[0];

	cout << "# of particles before annihilation = " << sp.size() << "\n";

	check_range(sp, Nx, Np);

	stable_sort(sp.begin(), sp.end(), comp);

	vector< vector< int > > sparse;

	sparse_matrix(sp, sparse);

	int count = 0;

	uniform_real_distribution<double> rnd(0, 1);

	double m = sp[0]->getm();
	double t = sp[0]->getTime();

	for (int i = 0; i < sp.size(); i++) {

		delete sp[i];
		sp[i] = nullptr;
	}

	for (int i = 0; i < (int)sparse.size(); i++) {

		for (int n = 0; n < fabs(sparse[i][2]); n++) {

			sp[count] = new Sp(count);

			sp[count]->setm(m);
			sp[count]->setXcell(sparse[i][0]);
			sp[count]->setx(x[sparse[i][0]] + (rnd(gen) - 0.5) * dx);
			sp[count]->setPcell(sparse[i][1]);
			sp[count]->setp(p[sparse[i][1]] + (rnd(gen) - 0.5) * dp);
			sp[count]->setWeight(sparse[i][2] > 0. ? 1 : -1);
			sp[count]->setTime(t);

			sp[count]->setStatus(true);

			count++;
		}
	}

	for (int i = count; i < sp.size(); i++) {

		delete sp[i];
		sp[i] = nullptr;
	}

	while (sp.back() == nullptr) {

		sp.pop_back();
	}

	cout << "# of particles after annihilation = " << sp.size() << "\n\n";

	return;
}