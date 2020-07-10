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

	~Bfield() = default;

public:
	std::unique_ptr<float[]> Bx;
	std::unique_ptr<float[]> By;
	std::unique_ptr<float[]> Bz;
};

#endif //REGIMPLICIT_B_FIELD_H
