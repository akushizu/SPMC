#ifndef SP_H
#define SP_H

#include <stdbool.h>

struct Sp
{
	int number;
	int id;

	double m;
	double x;
	int xcell;
	double p;
	int pcell;
	int weight;

	double time;

	bool status;
};

struct Spn
{
	int dimension;

	double* xn;
	int* xncell;
	double* pn;
	int* pncell;

	int number;
	int id;

	double m;
	double x;
	int xcell;
	double p;
	int pcell;
	int weight;

	double time;

	bool status;
};

#endif