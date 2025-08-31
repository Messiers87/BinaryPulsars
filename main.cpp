#include <iostream>
#include <fstream>
#include <filesystem>
#include <iomanip>
#include <cmath>
#include <vector>
#include <algorithm>
namespace fs = std::filesystem;


// ------- Physical constants
constexpr long double G = 6.67430e-11L;         // m^3 kg^-1 s^-2
constexpr long double MSUN = 1.98847e30L;          // kg
constexpr long double DAY = 86400.0L;             // seconds
constexpr long double PI = 3.141592653589793238462643383279502884L;

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

// Line-of-sight velocity v_l(t) (Eq. 24)
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

// Line-of-sight acceleration a_l(f) as a function of true anomaly (Eq. 26)
long double a_l_of_f(const OrbitParams& p, long double f) {
    long double Po = p.Po_day * DAY;
    long double e = p.e;
    long double ap = projected_semi_major_axis(p);
    long double varpi = deg2rad(p.varpi_deg);
    long double pref = -powl(2.0L * PI / Po, 2.0L) * ap / powl(1.0L - e * e, 2.0L);

    long double al = pref * sinl(f + varpi) * powl(1.0L + e * cosl(f), 2.0L); // m/s^2
    return al;
}


int main() {
    // --- Figure 2  ------------
    //Mp = 1.4 M⊙, Mc = 0.3 M⊙, Po = 0.5 day, i = 60◦,e = 0.5, Tp = 0, varpi = variable
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

    // ---- Figure 3 ------------
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

    std::cerr << "wrote fig2_velcoity.csv and fig3_acc.csv";
    return 0;


}