#include <vector>                                   // header for vector type
#include <fstream>                                  // header for csv output

int main() {
    double rho_max = 0.1;                           // max rho, cars / ft
    double v_max = 50.0;                            // max v, 50 ft / s

    double dx = 50.0;                               // cell size, ft
    int N_cars = 100;                               // # of cars, vehicles
    double t_f = 60;                                // simulation duration, s

    int N_x_stopped = N_cars / (dx * rho_max);      // # of cells for cars stopped on red
    int N_x_maxdx = (v_max * t_f) / dx;             // # of cells for max possible distance
    int N_x_bc = 1;                                 // # of cells for boundary conditons
    int N_x = N_x_stopped + N_x_maxdx + N_x_bc + 15;// # of cells needed + extra

    std::vector<double> S_f = {0.3, 0.6, 0.9};      // safety factors

    for (int j = 0; j < S_f.size(); j++) {   // output all flux and S_f combos in one loop
        double t = 0.0;                                 // simulation time, s
        double t_out = 20.0;                            // output time, s
        std::vector<double> rho(N_x, 0.0);              // cell densities
        std::vector<double> a(N_x, 0.0);                // cell characteristic speeds
        std::vector<double> flux(N_x, 0.0);             // cell fluxes

        for (int i = 1; i <= N_x_stopped; i++) {        // initial condition setup
            rho[i] = rho_max;
        }

        while (t < t_f) {                               // scheme
            double a_max = 0.0;

            for (int i = 0; i < N_x; i++) {             // characteristic speed calculation
                a[i] = v_max * (1.0 - (2.0 * rho[i]) / rho_max);
                if (std::abs(a[i]) > a_max) { a_max = std::abs(a[i]); }
            }

            double dt = S_f[j] * (dx / a_max);          // time step calculation, s

            for (int i = 0; i < N_x - 1; i++) {         // flux calculation
                double rho_L = rho[i];
                double rho_R = rho[i + 1];
                double a_L = a[i];
                double a_R = a[i + 1];

                if (a_L <= 0.0 && a_R >= 0.0) { flux[i] = rho_max * v_max / 4.0; }
                else if (a_L + a_R >= 0.0) { flux[i] = rho_L * v_max * (1.0 - rho_L / rho_max); }
                else { flux[i] = rho_R * v_max * (1.0 - rho_R / rho_max); }

                double rho_LW = 0.5 * (rho_L + rho_R) - 0.5 * (dt / dx) *
                (rho_R * v_max * (1.0 - rho_R / rho_max) - rho_L * v_max * (1.0 - rho_L / rho_max));
                flux[i] = rho_LW * v_max * (1.0 - rho_LW / rho_max);

                if (a_L + a_R >= 0.0) { flux[i] = rho_L * v_max * (1.0 - rho_L / rho_max); }
                else { flux[i] = rho_R * v_max * (1.0 - rho_R / rho_max); }

                if (a_L <= 0.0 && a_R >= 0.0) { flux[i] = 0.5 *
                (rho_L * v_max * (1.0 - rho_L / rho_max) + rho_R * v_max * (1.0 - rho_R / rho_max)
                - 0.5 * (a_R - a_L) * (rho_R - rho_L)); }
                else if (a_L + a_R >= 0.0) { flux[i] = rho_L * v_max * (1.0 - rho_L / rho_max); }
                else { flux[i] = rho_R * v_max * (1.0 - rho_R / rho_max); }
            }
            flux[N_x - 1] = rho[N_x - 1] * v_max * (1.0 - rho[N_x - 1] / rho_max);      // flux right BC

            for (int i = 1; i < N_x; i++) {             // rho update
                rho[i] = rho[i] - dt / dx * (flux[i] - flux[i - 1]);
            }
            t += dt;

            if (t > t_out) {                            // csv output
                std::ofstream output(std::format("rhovx_{}_0{}.csv", int(t), int(S_f[j] * 10.0)));
                for (int i = 0; i < N_x; i++) {
                    std::println(output, "{}, {}", dx * i + (dx / 2.0) - dx * (N_x_stopped + N_x_bc), rho[i]);
                }
                output.close();
                t_out += 20.0;
            }
        }
    }
}