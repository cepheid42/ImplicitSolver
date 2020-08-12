#ifndef REGIMPLICIT_B_FIELD_H
#define REGIMPLICIT_B_FIELD_H

#include "constants.cuh"

class Bfield {
public:
	Bfield() :
			Hx(new float[nz * ny * nx]{}),
			Hy(new float[nz * ny * nx]{}),
			Hz(new float[nz * ny * nx]{})
	{}

	~Bfield() {
		delete[] Hx;
		delete[] Hy;
		delete[] Hz;
	}

public:
	float* Hx;
	float* Hy;
	float* Hz;
};

#endif //REGIMPLICIT_B_FIELD_H
