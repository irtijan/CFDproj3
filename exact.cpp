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

    std::vector<double> t = {20.0, 40.0, 60.0};     // simulation times, s

    for (int j = 0; j < t.size(); j++) {            // output all times in one loop
        std::vector<double> rho(N_x, 0.0);              // cell densities
        std::vector<double> a(N_x, 0.0);                // cell characteristic speeds
        std::vector<double> flux(N_x, 0.0);             // cell fluxes

        for (int i = 1; i <= N_x; i++) {                 // initial condition setup
            double x = dx * i + (dx / 2.0) - dx * (N_x_stopped + N_x_bc);
            if (x / t[j] < -v_max) { rho[i] = rho_max; }
            else if (v_max < x / t[j]) { rho[i] = 0; }
            else { rho[i] = (rho_max / 2.0) * (1.0 - x / (v_max * t[j])); }
        }
                          
        std::ofstream output(std::format("rhovx_exact_{}.csv", int(t[j])));  // csv output
        for (int i = 0; i < N_x; i++) {
            std::println(output, "{}, {}", dx * i + (dx / 2.0) - dx * (N_x_stopped + N_x_bc), rho[i]);
        }
        output.close();
    }
}