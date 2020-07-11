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

	~Source() = default;

	// Verified
	void inc_ez(int time_step, int j, int k) const {
		int ind = get_index(0, j, k);

		float time = float(time_step) * dt;
		float time_dep_s = sinf(2.0f * pi * freq * time) * expf(-1.0f * powf(time - t0, n0) / (2.0f * powf(sig0, n0)));

		Jz[ind] = ez0 * time_dep_s * expf(-1.0f * powf(float(j) * dy, n0) / (2.0f * powf(4.0f * lambda, n0)));
	}

	// Verified
	void inc_ey(int time_step, int j, int k) const {
		int ind = get_index(0, j, k);
		float time = float(time_step) * dt;
		float time_dep_c = cosf(2.0f * pi * freq * time) * expf(-1.0f * powf(time - t0, n0) / (2.0f * powf(sig0, n0)));

		Jy[ind] = ez0 * time_dep_c * expf(-1.0f * powf(float(j) * dy, n0) / (2.0f * powf(4.0f * lambda, n0)));
	}

public:
	std::unique_ptr<float[]> Jx;
	std::unique_ptr<float[]> Jy;
	std::unique_ptr<float[]> Jz;
};

#endif //REGIMPLICIT_SOURCES_H
