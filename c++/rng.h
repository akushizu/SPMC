#ifndef RNG_H
#define RNG_H

#include <random>

using namespace std;

namespace rrr
{
	extern mt19937_64 gen;
	extern uniform_real_distribution<double> rnd(0, 1);
}

#endif 
