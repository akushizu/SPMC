#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>

using namespace std;

void saveRho(vector< Sp* >& sp, vector< int >& rho, vector< double >& t, int& ti, vector< double >& x, int Np) {

	int Nx = (int)x.size();
	int Nt = (int)t.size();
	
	int i, j;

	for (int n = 0; n < (int)sp.size(); n++) {

		i = sp[n]->getXcell();
		j = sp[n]->getPcell();

		if ((i >= 0 && i < Nx) && (j >= 0 && j < Np)) {

			rho[i] += sp[n]->getWeight();
		}
	}

	if (ti == 0) {

		ofstream den{ "rho.dat", ios::out };
		den << "x\trho" << "\n";
		den.close();
	}

	ofstream den{ "rho.dat", ios::app };

	for (i = 0; i < 20; i++) {

		den << "-";
	}

	den << "\n";

	den << "# = " << ti << " of " << Nt - 1 << "\n";
	den << "t = " << t[ti] << "\n";
	den << "\n";

	for (i = 0; i < Nx; i++) {

		den << fixed << setprecision(3) << x[i] << "\t" << rho[i] << "\n";
	}
	
	den << endl;
	den.close();

	return;
}

void saveFw(vector< Sp* >& sp, vector< vector< int > >& Fw, vector< double >& t, int& ti, vector< double >& x, vector< double >& p) {

	int Nx = (int)x.size();
	int Np = (int)p.size();
	int Nt = (int)t.size();

	int i, j;

	for (int n = 0; n < (int)sp.size(); n++) {

		i = sp[n]->getXcell();
		j = sp[n]->getPcell();

		if ((i >= 0 && i < Nx) && (j >= 0 && j < Np)) {

			Fw[i][j] += sp[n]->getWeight();
		}
	}
	
	if (ti == 0) {

		ofstream fw{ "Fw.dat", ios::out };
		fw << "x" << "\n";

		for (i = 0; i < Nx; i++) {

			fw << fixed << setprecision(3) << x[i] << "\t";
		}

		fw << "\n";
		fw << "p" << "\n";

		for (i = 0; i < Np; i++) {
			fw << fixed << setprecision(4) << p[i] << "\t";
		}

		fw << "\n";
		fw << "Fw" << "\n";
		fw.close();
	}

	ofstream fw{ "Fw.dat", ios::app };
	
	for (i = 0; i < 20; i++) {

		fw << "-";
	}

	fw << "\n";

	fw << "# = " << ti << " of " << Nt - 1 << "\n";
	fw << "t = " << t[ti] << "\n";
	fw << "\n";

	for (j = Np - 1; j > -1; j--) {

		for (i = 0; i < Nx; i++) {

			fw << Fw[i][j] << "\t";
		}

		fw << "\n";
	}

	fw << "\n" << endl;
	fw.close();

	return;
}

void save(vector< Sp* >& sp, vector< double >& t, int ti, vector< double >& x, vector< double >& p) {

	vector< int > rho(x.size(), 0);
	vector< vector< int > > Fw(x.size(), vector< int >(p.size(), 0));

	saveRho(sp, rho, t, ti, x, (int)p.size());
	saveFw(sp, Fw, t, ti, x, p);

	return;
}