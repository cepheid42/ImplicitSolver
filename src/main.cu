#include "constants.cuh"
#include "e_field.cuh"
#include "b_field.cuh"
#include "sources.cuh"
#include "update_functions.cuh"
#include "tdma.cuh"
#include "file_io.cuh"

using namespace std;

const int step = 10;


void run_loop(Efield& e, Bfield& b, Source& s, Tridiagonal& tdm) {
	Timer update_loop_timer;
	update_loop_timer.start();

	// Begin time loop
	for (int q = 0; q < nt; q++) {
		// Source
		for (int k = 0; k < nz; k++) {
			for (int j = 0; j < ny; j++) {
				auto ind = get_index(0, j, k);
				s.Jy[ind] = inc_ey(q, j);
				s.Jz[ind] = inc_ez(q, j);
			}
		}

		// c1 = dt / (2 * eps0)
		// c2 = dt / (2 * mu0)

		/* N -> N + 1/2 */
		// Implicit update
		implicit_e_half(e.ex_rhs, e.Ex, b.Bz, s.Jx, ddy); // ex = Ex + c1 * ddy * Bz - c1 * Jx
		implicit_e_half(e.ey_rhs, e.Ey, b.Bx, s.Jy, ddz); // ey = Ey + c1 * ddz * Bx - c1 * Jy
		implicit_e_half(e.ez_rhs, e.Ez, b.By, s.Jz, ddx); // ez = Ez + c1 * ddx * By - c1 * Jz

		tdm.TDMAsolver(e.ex_rhs, e.ex);
		tdm.TDMAsolver(e.ey_rhs, e.ey);
		tdm.TDMAsolver(e.ez_rhs, e.ez);

		// Explicit update
		explicit_e(e.Ex, e.ex);
		explicit_e(e.Ey, e.ey);
		explicit_e(e.Ez, e.ez);

		explicit_b_half(b.Bx, e.ey, ddz); // Bx = Bx + c2 * ddz * ey
		explicit_b_half(b.By, e.ez, ddx); // By = By + c2 * ddx * ez
		explicit_b_half(b.Bz, e.ex, ddy); // Bz = Bz + c2 * ddy * ex

		/* N + 1/2 -> N + 1 */
		// Implicit update
		implicit_e_one(e.ex_rhs, e.Ex, b.By, ddz); // ex = Ex - c1 * ddz * By
		implicit_e_one(e.ey_rhs, e.Ey, b.Bz, ddx); // ey = Ey - c1 * ddx * Bz
		implicit_e_one(e.ez_rhs, e.Ez, b.Bx, ddy); // ez = Ez - c1 * ddy * Bx

		// Explicit update
		explicit_e(e.Ex, e.ex);
		explicit_e(e.Ey, e.ey);
		explicit_e(e.Ez, e.ez);

		explicit_b_one(b.Bx, e.ez, ddy); // Bx = Bx - c2 * ddy * ez
		explicit_b_one(b.By, e.ex, ddz); // By = By - c2 * ddz * ex
		explicit_b_one(b.Bz, e.ey, ddx); // Bz = Bz - c2 * ddx * ey

		if (q % step == 0) {
			update_loop_timer.split();
			cout << q << "/" << nt << ": " << setprecision(3) << update_loop_timer.elapsed() << "s" << endl;
			snapshot(q, e, b);
		}
	}
	update_loop_timer.stop();
	cout << "Total loop time: " << setprecision(3) << update_loop_timer.time() << "s" << endl;
}

int main() {

	// todo:
	//  1.) Figure out why outputs are always 0
	//  2.) Plot slices
	//  3.) Output processing
	//  4.) Verify update functions work correctly


	save_params(step);

	Efield e;
	Bfield b;
	Source s;

	Tridiagonal tdm;
	tdm.init();

	run_loop(e, b, s, tdm);
	return 0;
}
