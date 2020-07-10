#ifndef REGIMPLICIT_SOURCES_H
#define REGIMPLICIT_SOURCES_H

#include "constants.cuh"

class Source {
public:
	Source() {
		Jx = new float[nz * ny * nx]();
		Jy = new float[nz * ny * nx]();
		Jz = new float[nz * ny * nx]();
	}

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
inline float inc_ez(int time_step, int y_loc) {
	float time = float(time_step) * dt;
	float time_dep_s = sinf(2.0f * pi * freq * time) * expf(-1.0f * powf(time - t0, n0) / (2.0f * powf(sig0, n0)));

	return ez0 * time_dep_s * expf(-1.0f * powf(float(y_loc), n0) / (2.0f * powf(4.0f * lambda, n0)));
}

// Verified
inline float inc_ey(int time_step, int y_loc) {
	float time = float(time_step) * dt;
	float time_dep_c = cosf(2.0f * pi * freq * time) * expf(-1.0f * powf(time - t0, n0) / (2.0f * powf(sig0, n0)));

	return ez0 * time_dep_c * expf(-1.0f * powf(float(y_loc), n0) / (2.0f * powf(4.0f * lambda, n0)));
}


#endif //REGIMPLICIT_SOURCES_H
