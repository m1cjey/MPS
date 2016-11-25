#ifndef HYPERELASTIC
#define HYPERELASTIC

#include <vector>


class hyperelastic
{
public:

	unsigned ID;
	double lambda;
	double half_p[DIMENSION];
	double stress[DIMENSION][DIMENSION];
	double differential_p[DIMENSION];
	double p[DIMENSION];
	int NEI[200];
	int N;
};

class hyperelastic2
{
public:
	unsigned ID;
	double DgDq[DIMENSION];
};

#endif