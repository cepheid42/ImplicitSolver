#ifndef REGIMPLICIT_SOURCES_H
#define REGIMPLICIT_SOURCES_H

#include "constants.cuh"

class Source {
public:
	Source() :
		Jx(new float[nz * ny * nx]{}),
		Jy(new float[nz * ny * nx]{}),
		Jz(new float[nz * ny * nx]{})
//		Jz(new float[nx]{})
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

// Verified
void inc_ey(float* Jy, int time_step) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			auto ind = get_index(0, j, k);
			float time = float(time_step) * dt;
			float time_dep_c = std::cos(2.0f * pi * freq * time) * std::exp(-1.0f * std::pow(time - t0, n0) / (2.0f * std::pow(sig0, n0)));

			Jy[ind] = time_dep_c * std::exp(-1.0f * std::pow(float(j) * dy, n0) / (2.0f * std::pow(4.0f * lambda, n0)));
		}
	}
}

// Verified
void inc_ez(float* Jz, int time_step) {
	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			auto ind = get_index(0, j, k);
			float time = float(time_step) * dt;
			float time_dep_s = std::sin(2.0f * pi * freq * time) * std::exp(-1.0f * std::pow(time - t0, n0) / (2.0f * std::pow(sig0, n0)));

			Jz[ind] = time_dep_s * std::exp(-1.0f * std::pow(float(j) * dy, n0) / (2.0f * std::pow(4.0f * lambda, n0)));
		}
	}
}

//void inc_ez_wave(float *Jz, int q) {
//	auto ind = get_index(0, int(ny / 2) + 1, int(nz / 2) + 1);
//	Jz[ind] = std::sin(2.0f * pi * freq * (float(q) * dt));
//}

#endif //REGIMPLICIT_SOURCES_H
