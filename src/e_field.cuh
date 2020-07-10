#ifndef REGIMPLICIT_E_FIELD_H
#define REGIMPLICIT_E_FIELD_H

#include "constants.cuh"

class Efield {
public:
	Efield() :
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

	~Efield() = default;

public:
	std::unique_ptr<float[]> Ex;
	std::unique_ptr<float[]> ex;
	std::unique_ptr<float[]> ex_rhs;

	std::unique_ptr<float[]> Ey;
	std::unique_ptr<float[]> ey;
	std::unique_ptr<float[]> ey_rhs;

	std::unique_ptr<float[]> Ez;
	std::unique_ptr<float[]> ez;
	std::unique_ptr<float[]> ez_rhs;
};

#endif //REGIMPLICIT_E_FIELD_H
