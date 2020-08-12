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

void update_sources(Source& s, float t) {
	float sigX = 10.0f * dx;
	float sigY = 10.0f * dy;
	float sigZ = 10.0f * dz;
	float X0 = dx * float(nx - 1) / 2.0f;
	float Y0 = dy * float(ny - 1) / 2.0f;
	float Z0 = dz * float(nz - 1) / 2.0f;
	float J0 = -140.0f * 0.007945f * (x_resolution / 96.0f);

	float time_dep_s = std::sin(2.0f * pi * freq * t) * std::exp(-1.0f * sqr(t - t0) / (2.0f * sqr(sig0)));

	for (int k = 0; k < nz; k++) {
		for (int j = 0; j < ny; j++) {
			for (int i = 0; i < nx; i++) {
				auto ind = get_index(i, j, k);

				float r2 = sqr(i * dx - X0) / (2.0f * sqr(sigX)) + sqr(j * dy - Y0) / (2.0f * sqr(sigY)) + sqr(k * dz - Z0) / (2.0f * sqr(sigZ));
				float sDep = std::exp(-r2);

				s.Jx[ind] = 0.0f;
				s.Jy[ind] = 0.0f;
				s.Jz[ind] = J0 * sDep * time_dep_s;
			}
		}
	}
}



#endif //REGIMPLICIT_SOURCES_H