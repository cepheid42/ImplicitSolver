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
//void inc_ey(float* Jy, int time_step) {
//	for (int j = 0; j < ny; j++) {
//		auto ind = get_index(0, j, int(nz / 2));
//		float time = float(time_step) * dt;
//		float time_dep_c = std::cos(2.0f * pi * freq * time) * std::exp(-1.0f * std::pow(time - t0, n0) / (2.0f * std::pow(sig0, n0)));
//
//		Jy[ind] = ez0 * time_dep_c * std::exp(-1.0f * std::pow(float(j) * dy, n0) / (2.0f * std::pow(4.0f * lambda, n0)));
//	}
//}
//
//// Verified
//void inc_ez(float* Jz, int time_step) {
//	for (int k = 0; k < nz; k++) {
//		auto ind = get_index(0, int(ny / 2), k);
//		float time = float(time_step) * dt;
//		float time_dep_s = std::sin(2.0f * pi * freq * time) * std::exp(-1.0f * std::pow(time - t0, n0) / (2.0f * std::pow(sig0, n0)));
//
//		Jz[ind] = ez0 * time_dep_s * std::exp(-1.0f * std::pow(float(k) * dy, n0) / (2.0f * std::pow(4.0f * lambda, n0)));
//	}
//}

void verify_jz(float* Jz, int q) {
	auto t = float(q) * dt;
	auto chi = (t - t0) / tau;

	auto ind = get_index(int(nx / 2), int(ny / 2), int(nz / 2));
	Jz[ind] = chi * std::exp(-1.0f * std::pow(chi, 2.0f));
}

#endif //REGIMPLICIT_SOURCES_H
