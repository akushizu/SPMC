#ifndef SP_H
#define SP_H

#include <vector>

using namespace std;

class Sp
{
protected:
	int number;
	int id = -1;

	double m;
	double x;
	int xcell;
	double p;
	int pcell;
	int weight;

	double time;

	bool status;

public:
	Sp(int n = -1, double M = -999.999, double X = -999.999, int xc = -999, double P = -999.999, int pc = -999, int w = 0, double t = -999.999, bool sta = false)
		: number(n), m(M), x(X), xcell(xc), p(P), pcell(pc), weight(w), time(t), status(sta)
	{

	}

	void setNumber(int n)
	{
		number = n;
	}

	int getNumber() const
	{
		return number;
	}

	void setm(double n)
	{
		m = n;
	}

	double getm() const
	{
		return m;
	}

	void setx(double r)
	{
		x = r;
	}

	double getx() const
	{
		return x;
	}

	void setXcell(int xn)
	{
		xcell = xn;
	}

	int getXcell() const
	{
		return xcell;
	}

	void setp(double k)
	{
		p = k;
	}

	double getp() const
	{
		return p;
	}

	void setPcell(int pn)
	{
		pcell = pn;
	}

	int getPcell() const
	{
		return pcell;
	}

	void setWeight(int n)
	{
		weight = n;
	}

	int getWeight() const
	{
		return weight;
	}

	void setTime(double t)
	{
		time = t;
	}

	double getTime() const
	{
		return time;
	}

	void setStatus(bool n)
	{
		status = n;
	}

	bool getStatus() const
	{
		return status;
	}
};

class Spn : public Sp
{
private:
	int dimension;

	vector< double > xn;
	vector< int > xncell;
	vector< double > pn;
	vector< int > pncell;

public:
	Spn(int n = -1, double M = -999.999, int w = 0, double t = -999.999, bool sta = false)
		: dimension(2), Sp(n, M, -999.999, -999, -999.999, -999, w, t, sta)
	{
		for (int i = 1; i <= dimension; i++) {

			xn.push_back(-999.999);
			xncell.push_back(-999);
			pn.push_back(-999.999);
			pncell.push_back(-999);
		}
	}

	void setx(vector< double > r)
	{
		xn = r;
	}

	vector< double > getx() const
	{
		return xn;
	}

	void setXcell(vector< int > xn)
	{
		xncell = xn;
	}

	vector< int > getXcell() const
	{
		return xncell;
	}

	void setp(vector< double > k)
	{
		pn = k;
	}

	vector< double > getp() const
	{
		return pn;
	}

	void setPcell(vector< int > pn)
	{
		pncell = pn;
	}

	vector< int > getPcell() const
	{
		return pncell;
	}
};

#endif