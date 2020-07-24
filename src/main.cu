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

	Tridiagonal td_x_half(ny, ddy);
	Tridiagonal td_y_half(nz, ddz);
	Tridiagonal td_z_half(nx, ddx);

	Tridiagonal td_x_one(nz, ddz);
	Tridiagonal td_y_one(nx, ddx);
	Tridiagonal td_z_one(ny, ddy);

	// Begin time loop
	for (int q = 0; q < nt; q++) {
		// Sources
//		inc_ey(e.Ey, q);
//		inc_ez(e.Ez, q);
		auto ind = true_middle();
		auto a = ((float(q) * dt) - t0) / tau;
		s.Jz[ind] = a * exp(-1.0f * (a * a));

		// c1 = dt / (2 * eps0)
		// c2 = dt / (2 * mu0)

		/* N -> N + 1/2 */
		// Implicit update
		implicit_ex_half(e.ex_rhs, e.Ex, b.Bz, s.Jx); // ex = Ex + c1 * ddy * Bz - c1 * Jx
		implicit_ey_half(e.ey_rhs, e.Ey, b.Bx, s.Jy); // ey = Ey + c1 * ddz * Bx - c1 * Jy
		implicit_ez_half(e.ez_rhs, e.Ez, b.By, s.Jz); // ez = Ez + c1 * ddx * By - c1 * Jz

		ddy_solve(td_x_half, e.ex_rhs, e.ex);
		ddz_solve(td_y_half, e.ey_rhs, e.ey);
		ddx_solve(td_z_half, e.ez_rhs, e.ez);


		// Explicit update
		explicit_E(e.Ex, e.ex);
		explicit_E(e.Ey, e.ey);
		explicit_E(e.Ez, e.ez);

		explicit_bx_half(b.Bx, e.ey); // Bx = Bx + c2 * ddz * ey
		explicit_by_half(b.By, e.ez); // By = By + c2 * ddx * ez
		explicit_bz_half(b.Bz, e.ex); // Bz = Bz + c2 * ddy * ex


		/* N + 1/2 -> N + 1 */
		// Implicit update
		implicit_ex_one(e.ex_rhs, e.Ex, b.By); // ex = Ex - c1 * ddz * By
		implicit_ey_one(e.ey_rhs, e.Ey, b.Bz); // ey = Ey - c1 * ddx * Bz
		implicit_ez_one(e.ez_rhs, e.Ez, b.Bx); // ez = Ez - c1 * ddy * Bx

		ddz_solve(td_x_one, e.ex_rhs, e.ex);
		ddx_solve(td_y_one, e.ey_rhs, e.ey);
		ddy_solve(td_z_one, e.ez_rhs, e.ez);

		// Explicit update
		explicit_E(e.Ex, e.ex);
		explicit_E(e.Ey, e.ey);
		explicit_E(e.Ez, e.ez);

		explicit_bx_one(b.Bx, e.ez); // Bx = Bx - c2 * ddy * ez
		explicit_by_one(b.By, e.ex); // By = By - c2 * ddz * ex
		explicit_bz_one(b.Bz, e.ey); // Bz = Bz - c2 * ddx * ey

		if (q % step == 0) {
			cout << q << "/" << nt << ": snapshot taken (" << update_loop_timer.split() << "s)" << endl;
			snapshot(q, e, b, s);
		}
	}
	update_loop_timer.stop();
	cout << "Total loop time: " << setprecision(3) << update_loop_timer.total << "s" << endl;
}

int main() {
	save_params(step);

	Efield e;
	Bfield b;
	Source s;

	run_loop(e, b, s);
	return 0;
}
