#ifndef REGIMPLICIT_B_FIELD_H
#define REGIMPLICIT_B_FIELD_H

#include "constants.cuh"

class Bfield {
public:
	Bfield() :
		Bx(new float[nz * ny * nx]{}),
		By(new float[nz * ny * nx]{}),
		Bz(new float[nz * ny * nx]{})
	{}

	~Bfield() {
		delete[] Bx;
		delete[] By;
		delete[] Bz;
	}

public:
	float* Bx;
	float* By;
	float* Bz;
};

#endif //REGIMPLICIT_B_FIELD_H
