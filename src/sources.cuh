#ifndef REGIMPLICIT_SOURCES_H
#define REGIMPLICIT_SOURCES_H

#include "constants.cuh"

class Source {
public:
	Source() :
		Jx(new float[nz * ny * nx]{}),
		Jy(new float[nz * ny * nx]{}),
		Jz(new float[nz * ny * nx]{})
	{}

	~Source() {
		delete[] Jx;
		delete[] Jy;
		delete[] Jz;
	}

public:
	float* Jx;
	float* Jy;
	float* Jz;
};

void inc_ey(float* Ey, int time_step) {
	for (int j = 0; j < ny; j++) {
		auto ind = get_index(0, j, int(nz / 2));
		float time = float(time_step) * dt;
		float time_dep_c = std::cos(2.0f * pi * freq * time) * std::exp(-1.0f * std::pow(time - t0, n0) / (2.0f * std::pow(sig0, n0)));

		Ey[ind] = ez0 * time_dep_c * std::exp(-1.0f * std::pow(float(j) * dy, n0) / (2.0f * std::pow(4.0f * lambda, n0)));
	}
}

void inc_ez(float* Ez, int time_step) {
	for (int k = 0; k < nz; k++) {
		auto ind = get_index(0, int(ny / 2), k);
		float time = float(time_step) * dt;
		float time_dep_s = std::sin(2.0f * pi * freq * time) * std::exp(-1.0f * std::pow(time - t0, n0) / (2.0f * std::pow(sig0, n0)));

		Ez[ind] = ez0 * time_dep_s * std::exp(-1.0f * std::pow(float(k) * dy, n0) / (2.0f * std::pow(4.0f * lambda, n0)));
	}
}


#endif //REGIMPLICIT_SOURCES_H
