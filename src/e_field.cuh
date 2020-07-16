#ifndef REGIMPLICIT_E_FIELD_H
#define REGIMPLICIT_E_FIELD_H

#include "constants.cuh"

class Efield {
public:
	Efield() :
//		Ez    (new float[nx]{}),
//		ez    (new float[nx]{}),
//		ez_rhs(new float[nx]{})
		Ex    (new float[nz * ny * nx]{}),
		ex    (new float[nz * ny * nx]{}),
		ex_rhs(new float[nz * ny * nx]{}),

		Ey    (new float[nz * ny * nx]{}),
		ey    (new float[nz * ny * nx]{}),
		ey_rhs(new float[nz * ny * nx]{}),

		Ez    (new float[nz * ny * nx]{}),
		ez    (new float[nz * ny * nx]{}),
		ez_rhs(new float[nz * ny * nx]{})
	{}

	~Efield() {
		delete[] Ex;
		delete[] ex;
		delete[] ex_rhs;

		delete[] Ey;
		delete[] ey;
		delete[] ey_rhs;

		delete[] Ez;
		delete[] ez;
		delete[] ez_rhs;
	}

public:
	float* Ex;
	float* ex;
	float* ex_rhs;

	float* Ey;
	float* ey;
	float* ey_rhs;

	float* Ez;
	float* ez;
	float* ez_rhs;
};

#endif //REGIMPLICIT_E_FIELD_H
