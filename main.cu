#include <iostream>

#include "constants.cuh"

class Efield {
public:
	Efield() : Ez(nullptr) {}
	~Efield() = default;

	void create() {
		checkErr(cudaMallocManaged(&Ez, nx * sizeof(float)))
		checkErr(cudaDeviceSynchronize())
	}

	void destroy() {
		checkErr(cudaDeviceSynchronize())
		checkErr(cudaFree(Ez))
	}

	void zero() {
		for (int i = 0; i < nx; i++) {
			Ez[i] = 0.0f;
		}
	}

public:
//	float *Ex;
//	float *Ey;
	float *Ez;
};

class Bfield {
public:
	Bfield() : Bx(nullptr) {}
	~Bfield() = default;

	void create() {
		checkErr(cudaMallocManaged(&Bx, nx * sizeof(float)))
		checkErr(cudaDeviceSynchronize())
	}

	void destroy() {
		checkErr(cudaDeviceSynchronize())
		checkErr(cudaFree(Bx))
	}

	void zero() {
		for (int i = 0; i < nx; i++) {
			Bx[i] = 0.0f;
		}
	}
public:
	float* Bx;
//	float* By;
//	float* Bz;
};


int main() {

	Efield e;
	e.create();
	e.zero();

	Bfield b;
	b.create();
	b.zero();

	return 0;
}
