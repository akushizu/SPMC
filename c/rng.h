#ifndef RNG_H
#define RNG_H

#include <stdlib.h>

double rnd() {

	return (double)rand() / (double)RAND_MAX;
}

#endif 
