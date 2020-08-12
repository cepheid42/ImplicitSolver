#include "constants.cuh"
#include "e_field.cuh"
#include "b_field.cuh"
#include "sources.cuh"
#include "update_functions.cuh"
#include "tdma.cuh"
#include "file_io.cuh"

using namespace std;

const int step = 10;

void run_loop(Efield& e, Bfield& b, Source& s) {
	Timer update_loop_timer;
	update_loop_timer.start();

	Tridiagonal td_x_half(nx, dy); // ex n+1/2
	Tridiagonal td_y_half(ny, dz); // ey n+1/2
	Tridiagonal td_z_half(nz, dx); // ez n+1/2

	Tridiagonal td_x_one(nx, dz); // ex n+1
	Tridiagonal td_y_one(ny, dx); // ey n+1
	Tridiagonal td_z_one(nz, dy); // ez n+1

	// Begin time loop
	for (int q = 0; q < nt; q++) {
		auto t = float(q) * dt;
		update_sources(s, t);

		/* n -> n + 1/2 */
		// Implicit e update
		implicit_ex_half(e, b, s);
		x_solve(td_x_half, e.ex_rhs, e.ex);

		implicit_ey_half(e, b, s);
		y_solve(td_y_half, e.ey_rhs, e.ey);

		implicit_ez_half(e, b, s);
		z_solve(td_z_half, e.ez_rhs, e.ez);

		// Explicit E update
		explicit_E(e.Ex, e.ex);
		explicit_E(e.Ey, e.ey);
		explicit_E(e.Ez, e.ez);

		// Explicit H update
		explicit_Hx_half(b, e);
		explicit_Hy_half(b, e);
		explicit_Hz_half(b, e);

		/* n + 1/2 -> n + 1 */
		// Implicit e update
		implicit_ex_one(e, b);
		x_solve(td_x_one, e.ex_rhs, e.ex);

		implicit_ey_one(e, b);
		y_solve(td_y_one, e.ey_rhs, e.ey);

		implicit_ez_one(e, b);
		z_solve(td_z_one, e.ez_rhs, e.ez);

		// Explicit E update
		explicit_E(e.Ex, e.ex);
		explicit_E(e.Ey, e.ey);
		explicit_E(e.Ez, e.ez);

		// Explicit H update
		explicit_Hx_one(b, e);
		explicit_Hy_one(b, e);
		explicit_Hz_one(b, e);

		if (q % step == 0) {
			cout << q << "/" << nt << ": " << update_loop_timer.split() << endl;
			snapshot(q, e, b);
		}
	}

	cout << "Update loop completed. Processing arrays." << endl;
	snapshot(nt, e, b);
	update_loop_timer.stop();
	cout << "Processing completed.\nTotal time: " << update_loop_timer.total << endl;
}

int main() {
	save_params(step);

	Efield e;
	Bfield b;
	Source s;

	run_loop(e, b, s);
	return 0;
}
