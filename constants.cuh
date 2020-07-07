#ifndef IMPLICITSOLVER_CONSTANTS_CUH
#define IMPLICITSOLVER_CONSTANTS_CUH

const float c0 = 299792458.0;  // speed of light, m/s
const float eps0 = 8.85418782e-12;
const float mu0 = 1.25663706e-6;
const float pi = 3.1415926535;

const float freq = 14.0e6;
const float lamb = c0 / freq;
const float T = 1.0f / freq;
const float t0 = 6.0f * T;
const float sig0 = 2.0f * T;
const float n0 = 2.0f;

const int x_resolution = 64;
const int num_wavelengths = 20;

const int nx = num_wavelengths * x_resolution + 1;
const float dx = (num_wavelengths * lamb - 0.0f) / float(nx - 1);

const int ny = nx;
const float dy = dx;

const int nz = nx;
const float dz = dx;


// courant number
const float C = 1.0;

const float dt = C * dx / c0;
const int nt = int(num_wavelengths * lamb / c0 / dt) + 1;

#define checkErr(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true) {
	if (code != cudaSuccess) {
		fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
		if (abort) exit(code);
	}
}

#endif //IMPLICITSOLVER_CONSTANTS_CUH
