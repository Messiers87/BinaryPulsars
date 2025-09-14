#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <cmath>
#include <vector>
#include <complex>
#include <algorithm>
#include <chrono>
namespace fs = std::filesystem;


// ------- Physical constants
constexpr long double G = 6.67430e-11L;         // m^3 kg^-1 s^-2
constexpr long double MSUN = 1.98847e30L;          // kg
constexpr long double DAY = 86400.0L;             // seconds
constexpr long double PI = 3.141592653589793238462643383279502884L;
constexpr long double C = 299792458.0L;          // m/s

// ------- Angle helpers
inline long double deg2rad(long double d) { return d * PI / 180.0L; }
inline long double rad2deg(long double r) { return r * 180.0L / PI; }

//normalising function
inline long double norm2pi(long double x) {
    long double y = fmod(x, 2.0L * PI);
    if (y < 0) y += 2.0L * PI;
    return y;
}

//solving kepler equation, E - e*sin(E) = M
long double solve_kepler_E(long double M, long double e) {
    M = norm2pi(M);
    if (e == 0.0L) return M;

    long double E = (e < 0.8L) ? M : PI;
    for (int it = 0; it < 60; ++it) {
        long double f = E - e * sinl(E) - M;
        long double fp = 1.0L - e * cosl(E);
        long double dE = -f / fp;
        E += dE;
        if (fabsl(dE) < 1e-14L) break;
    }
    return norm2pi(E);
}

//function to get true anomaly from eccentric anomaly
long double E_to_true_f(long double E, long double e) {
    long double s = sqrtl(1.0L + e) * sinl(E / 2.0L);
    long double c = sqrtl(1.0L - e) * cosl(E / 2.0L);
    long double f = 2.0L * atan2l(s, c);
    return norm2pi(f);
}

struct OrbitParams {
    long double Mp_sun;   // mass of pulsar [M_sun]
    long double Mc_sun;   // mass of companion [M_sun]
    long double Po_day;   // orbital period [day]
    long double inc_deg;  // inclination i [deg]
    long double e;        // eccentricity
    long double Tp_s;     // epoch of periastron [s]
    long double varpi_deg;// longitude of periastron ̟ [deg]
};

// Projected semi-major axis a_p' (Eq. 25)
long double projected_semi_major_axis(const OrbitParams& p) {
    long double Mp = p.Mp_sun * MSUN;
    long double Mc = p.Mc_sun * MSUN;
    long double Po = p.Po_day * DAY;
    long double sin_i = sinl(deg2rad(p.inc_deg));

    // equation 25
    long double ap_prime = cbrtl((powl(Po / (2.0L * PI), 2.0L) * G * (Mp + Mc))) * (Mc / (Mp + Mc)) * sin_i; // meters
    return ap_prime;
}

// LOS position at time t (equation 22)
long double r_l(const OrbitParams& p, long double t) {
    long double Po = p.Po_day * DAY;
    long double omega_o = 2.0L * PI / Po;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);

    // Kepler: M = ω_o (t - T_p), solve for E, then true anomaly f
    long double M = omega_o * (t - p.Tp_s);
    long double E = solve_kepler_E(M, e);
    long double f = E_to_true_f(E, e);
    long double varpi = deg2rad(p.varpi_deg);

    long double rl = ap * sinl(f + varpi) * (1.0L - e * e) / (1.0L + e * cosl(f));
    return rl;
}

// LOS velocity v_l(t) (equation 24)
long double v_l(const OrbitParams& p, long double t) {
    long double Po = p.Po_day * DAY;
    long double omega_o = 2.0L * PI / Po;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);

    // Kepler: M = ω_o (t - T_p), solve for E, then true anomaly f
    long double M = omega_o * (t - p.Tp_s);
    long double E = solve_kepler_E(M, e);
    long double f = E_to_true_f(E, e);
    long double varpi = deg2rad(p.varpi_deg);
    long double pref = (2.0L * PI / Po) * ap / sqrtl(1.0L - e * e);

    //v_l(t)
    long double vl = pref * (cosl(f + varpi) + e * cosl(varpi)); // m/s
    return vl;
}

// LOS acceleration a_l(f) as a function of true anomaly (Eq. 26)
long double a_l_of_f(const OrbitParams& p, long double f) {
    long double Po = p.Po_day * DAY;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);
    long double varpi = deg2rad(p.varpi_deg);
    long double pref = -powl(2.0L * PI / Po, 2.0L) * ap / powl(1.0L - e * e, 2.0L);

    long double al = pref * sinl(f + varpi) * powl(1.0L + e * cosl(f), 2.0L); // m/s^2
    return al;
}

// ++++++++++ gamma function (Eq. 17) 
long double gamma1(const OrbitParams& baseParams, long double f, long double T_obs, int m, long double Pp_s) {

    OrbitParams p = baseParams; 

    long double Po = p.Po_day * DAY;
    long double omega_o = 2.0L * PI / Po;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);

    // compute Tp so that f(0) = f0
    long double f0 = deg2rad(f);
    long double tan_halfE = tanl(f0 / 2.0L) * sqrtl((1.0L - e) / (1.0L + e));
    long double E0 = 2.0L * atanl(tan_halfE);

    E0 = norm2pi(E0);
    long double M0 = E0 - e * sinl(E0);
    p.Tp_s = -M0 / omega_o; 

    // time sampling
    const int N = 3000; // reduced resolution (was 5000)
    std::vector<long double> t_grid(N + 1);
    std::vector<long double> r_grid(N + 1);
    for (int i = 0; i <= N; i++) {
        t_grid[i] = T_obs * ((long double)i) / ((long double)N);
        r_grid[i] = r_l(p, t_grid[i]);
    }

    long double r0 = r_grid[0];

    // constants
    long double omega_p = 2.0L * PI / Pp_s;
    long double K = m * omega_p / C;

    // velocity bounds
    long double vmin = 1.0e10L;
    long double vmax = -1.0e10L;
    for (int i = 1; i <= N; i++) {
        long double v = (r_grid[i] - r_grid[i - 1]) / (t_grid[i] - t_grid[i - 1]);
        if (v < vmin) vmin = v;
        if (v > vmax) vmax = v;
    }
    long double margin = 1e3L; // add 1 km/s margin

    // helper lambda for trapezoidal integration
    auto compute_gamma = [&](long double alpha) {
        std::complex<long double> sum = 0.0L;
        for (int j = 0; j < N; j++) {
            long double phase1 = K * (r_grid[j] - r0 - alpha * t_grid[j]);
            long double phase2 = K * (r_grid[j + 1] - r0 - alpha * t_grid[j + 1]);
            std::complex<long double> e1 = std::polar(1.0L, phase1);
            std::complex<long double> e2 = std::polar(1.0L, phase2);
            sum += 0.5L * (e1 + e2) * (t_grid[j + 1] - t_grid[j]);
        }
        return std::abs(sum) / T_obs;
        };

    // === Stage 1: coarse scan ===
    int NA_coarse = 100;
    long double alpha_best = 0.0L;
    long double gamma_best = -1.0L;

    for (int i = 0; i < NA_coarse; i++) {
        long double alpha = (vmin - margin) +
            (vmax - vmin + 2.0L * margin) * ((long double)i) / ((long double)NA_coarse - 1);
        long double gamma = compute_gamma(alpha);
        if (gamma > gamma_best) {
            gamma_best = gamma;
            alpha_best = alpha;
        }
    }

    // === Stage 2: refine scan around best alpha ===
    int NA_refine = 100;
    long double refine_width = (vmax - vmin) / 50.0L; // narrower window
    long double alpha_ref_best = alpha_best;
    long double gamma_ref_best = gamma_best;

    for (int i = 0; i < NA_refine; i++) {
        long double alpha = alpha_best - refine_width / 2.0L +
            refine_width * ((long double)i) / ((long double)NA_refine - 1);
        long double gamma = compute_gamma(alpha);
        if (gamma > gamma_ref_best) {
            gamma_ref_best = gamma;
            alpha_ref_best = alpha;
        }
    }

    return gamma_ref_best;
}



int main() {
    auto t1 = std::chrono::high_resolution_clock::now();
    // 1. ======== Figure 2
    
    //Mp = 1.4 , Mc = 0.3, Po = 0.5 day, i = 60,e = 0.5, Tp = 0, varpi = variable
    OrbitParams baseFig2{ 1.4L, 0.3L, 0.5L, 60.0L, 0.5L, 0.0L, 0.0L };

    std::vector<long double> varpi_upper = { 0.0L, 90.0L, 180.0L, 270.0L };
    std::vector<long double> varpi_lower = { 60.0L, 120.0L, 240.0L, 300.0L };

    //Sampling over full orbital Period
    const int N = 1000; // resolution
    long double Po_sec = baseFig2.Po_day * DAY;

    //generate csv file
    fs::create_directories("data");
    std::ofstream f2("data/fig2_velocity.csv");
    f2.setf(std::ios::fixed); f2 << std::setprecision(6);

    // Header
    f2 << "t_sec";
    //appending columns (upper panel)
    for (auto w : varpi_upper) f2 << ",v_varpi_" << (int)w << "deg";
    //appending columns (lower panel)
    for (auto w : varpi_lower) f2 << ",v_varpi_" << (int)w << "deg";
    //appending column (circular case)
    f2 << ",v_circular" << '\n';

    // Prepare circular case: e=0, varpi arbitrary (set 0)
    OrbitParams circ = baseFig2; circ.e = 0.0L; circ.varpi_deg = 0.0L;

    for (int i = 0; i <= N; i++) {
        long double t = Po_sec * ((long double)i) / ((long double)N);
        f2 << t;

        //upper panel
        for (auto w : varpi_upper) {
            auto p = baseFig2; p.varpi_deg = w;
            f2 << "," << v_l(p, t);

        }

        //lower panel
        for (auto w : varpi_lower) {
            auto p = baseFig2;
            p.varpi_deg = w;
            f2 << "," << v_l(p, t);
        }

        //circular case
        f2 << "," << v_l(circ, t);
        f2 << '\n';
    }
    f2.close();

    // 2. ======== Figure 3 
    
    // Mp = 1.4 M⊙, Mc = 0.3 M⊙,Po = 0.1 day i = 60◦, e = 0.5, Tp=0, varpi = 0◦  
    OrbitParams fig3{ 1.4L, 0.3L, 0.1L, 60.0L, 0.5L, 0.0L, 0.0L };

    std::ofstream f3("data/fig3_acc.csv");
    f3.setf(std::ios::fixed); f3 << std::setprecision(6);
    f3 << "f_deg,a_l__m/s2" << '\n';

    const int NF = 360; // 1-degree steps
    for (int k = 0; k <= NF; k++) {
        long double fdeg = (long double)k;
        long double f = deg2rad(fdeg);
        f3 << fdeg << "," << a_l_of_f(fig3, f) << '\n';
    }

    f3.close();

	// 3. ========= gammma function

	OrbitParams gammaParams{ 1.4L, 0.3L, 0.1L, 60.0L, 0.5L, 0.0L, 0.0L };

    long double T_obs = 500.0L;
	int m = 4;
	long double Pp_s = 0.01L;

    // list of starting true anomalies
    std::vector<long double> f0_list = { 0.0L, 10.0L, 40.0L, 50.0L, 60.0L, 70.0L, 170.0L, 180.0L, 190.0L, 350.0L };

	// ensure the output directory exists
	fs::create_directories("data");

	// open output file
    std::ofstream table1("data/table1_results.csv");
	table1.setf(std::ios::fixed); table1 << std::setprecision(6);
	table1 << "f0_deg,gamma1,w" << '\n';



    auto compute_w = [&](long double f0_deg, long double e) {
        long double f0 = deg2rad(f0_deg);
        long double num = powl(1.0L + e * cosl(f0), -2.0L);
        long double denom = powl(1.0L + e * cosl(PI), -2.0L);
        return num / denom;
        };

	// compute gamma1 for each f0
    for (auto f0 : f0_list) {
        std::cout << "loop 2: main func gamma loop start" << std::endl;
        std::cout << "-------------" << std::endl;
		std::cout << "Computing gamma1 for f0 = " << f0 << std::endl;
        long double gamma = gamma1(gammaParams, f0, T_obs, m, Pp_s);

        long double w = compute_w(f0, gammaParams.e);

        table1 << f0 << "," << gamma << w << '\n';
        std::cout << "f0=" << f0 << "deg -> gamma1 =" << gamma << "w = "<< w << std::endl;

    }

    table1.close();

    // 4. ====== figure 4 data

	//Mp, Mc, Po, i, e, Tp, varpi
    OrbitParams fig4{ 1.4L, 0.3L, 0.01L, 60.0L, 0.5L, 0.0L, 0.0L };

	fs::create_directories("data");
	std::ofstream f4("data/fig4_data.csv");
	f4.setf(std::ios::fixed); f4 << std::setprecision(6);
	f4 << "Po_day,Pp_s,gamma1" << '\n';

    std::vector<long double> Po_days = { 0.1L, 0.2L, 0.5L, 1.0L, 10.0L, 100.0L};
    std::vector<long double> Pp_s_vec = { 0.001L, 0.005L, 0.01L, 0.05L, 0.1L, 1.0L };

    long double T_obs_fig4 = 1000.0L;
    int m_fig4 = 4;

    for (auto Po_d : Po_days) {
		std::cout << "-------------" << std::endl;
        std::cout << "loop 3: figure 4 data loop" << std::endl;
            for (auto Pp_s_val : Pp_s_vec) {
            fig4.Po_day = Po_d;

            long double sum_num = 0.0L;
            long double sum_den = 0.0L;
            for (int fdeg = 0; fdeg < 360; fdeg += 2) {   // greater than 360 points
                long double gamma_f = gamma1(fig4, (long double)fdeg, T_obs_fig4, m_fig4, Pp_s_val);
                long double w = powl(1.0L + fig4.e * cosl(deg2rad((long double)fdeg)), -2.0L);
                sum_num += gamma_f * w;
                sum_den += w;
            }
            long double gamma_avg = sum_num / sum_den;

            f4 << Po_d << "," << Pp_s_val << "," << gamma_avg << '\n';
            
        }; 
            std::cout << "===== end of loop 3 ===== " << std::endl;
    };

    f4.close();

    //after computing gamma1 
	std::cerr << "wrote data/fig4_data.csv with gamma1 for f0=180deg !\n";

    auto t2 = std::chrono::high_resolution_clock::now();
    std::cerr << "gamma1 call look" << std::chrono::duration<double>(t2 - t1).count() << "s\n";
    std::cout << "end of program" << std::endl;
    std::cin.get();

    return 0;


}


// 1094 sec