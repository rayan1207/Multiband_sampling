#include "mini_ami.hpp"
using namespace std;
double mband::beta(const std::complex<double>& A, const std::complex<double>& B, const std::complex<double>& C,
    const std::complex<double>& g, const std::complex<double>& n, const std::complex<double>& s) {
    std::complex<double > e1= pow(2.0 * pow(A, 3) - 3.0 * pow(A, 2) * B - 3.0 * pow(A, 2) * C + sqrt(4.0 * pow(-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0), 3.0) +
        pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * (pow(A, 2.0)) * C - 3.0 * A * (pow(B, 2.0)) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) -
            3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 2.0))
        - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) -
        18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 1.0 / 3.0) / (3.0 * pow(2.0, 1.0 / 3.0)) - (pow(2.0, 1 / 3.0) * (-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) +
            B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0))) / (3.0 * pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C + sqrt(4.0 * pow(-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) -
                3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0), 3.0) + pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) +
                    9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) +
                    54.0 * pow(n, 3.0) * pow(s, 2.0), 2.0)) - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) +
                9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 1.0 / 3.0)) + (1.0 / 3.0) * (A + B + C);
    return std::real(e1);
}
double mband::alpha(const std::complex<double>& A, const std::complex<double>& B, const std::complex<double>& C,
    const std::complex<double>& g, const std::complex<double>& n, const std::complex<double>& s) {
    std::complex<double > E2 = -(complex<double>(1.0, -sqrt(3.0)) * pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C + sqrt(4.0 * pow(-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0), 3.0)
        + pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0)
            + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 2.0)) - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0)
        - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) +
        9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 1.0 / 3.0)) / (6.0 * pow(2.0, 1.0 / 3.0)) + (complex<double>(1.0, sqrt(3.0)) * (-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0))) / (3.0 * pow(2.0, 2.0 / 3.0) * pow(2.0 * pow(A, 3.0)
            - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C + sqrt(4.0 * pow(-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0), 3.0) + pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C
                - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0)
                + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3) * pow(s, 2), 2.0)) - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) -
            18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 1.0 / 3.0)) + (1.0 / 3.0) * (A + B + C);
    return std::real(E2);
}
double mband::gamma(const std::complex<double>& A, const std::complex<double>& B, const std::complex<double>& C,
    const std::complex<double>& g, const std::complex<double>& n, const std::complex<double>& s) {
    std::complex<double > E2 = -(complex<double>(1.0, sqrt(3.0)) * pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C + sqrt(4.0 * pow(-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0), 3.0)
        + pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0)
            + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 2.0)) - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0)
        - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) +
        9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 1.0 / 3.0)) / (6.0 * pow(2.0, 1.0 / 3.0)) + (complex<double>(1.0, -sqrt(3.0)) * (-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0))) / (3.0 * pow(2.0, 2.0 / 3.0) * pow(2.0 * pow(A, 3.0)
            - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C + sqrt(4.0 * pow(-pow(A, 2.0) + A * B + A * C - pow(B, 2.0) + B * C - pow(C, 2.0) - 3.0 * pow(g, 2.0) - 6.0 * pow(n, 2.0) * pow(s, 2.0) - 3.0 * pow(n, 2.0), 3.0) + pow(2.0 * pow(A, 3.0) - 3.0 * pow(A, 2.0) * B - 3.0 * pow(A, 2.0) * C - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C
                - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) - 18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0)
                + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3) * pow(s, 2), 2.0)) - 3.0 * A * pow(B, 2.0) + 12.0 * A * B * C - 3.0 * A * pow(C, 2.0) + 9.0 * A * pow(g, 2.0) - 9.0 * A * pow(n, 2.0) * pow(s, 2.0) + 9.0 * A * pow(n, 2.0) + 2.0 * pow(B, 3.0) - 3.0 * pow(B, 2.0) * C - 3.0 * B * pow(C, 2.0) + 9.0 * B * pow(g, 2.0) + 18.0 * B * pow(n, 2.0) * pow(s, 2.0) -
            18.0 * B * pow(n, 2.0) + 2.0 * pow(C, 3.0) - 18.0 * C * pow(g, 2.0) - 9.0 * C * pow(n, 2.0) * pow(s, 2.0) + 9.0 * C * pow(n, 2.0) + 54.0 * pow(n, 3.0) * pow(s, 2.0), 1.0 / 3.0)) + (1.0 / 3.0) * (A + B + C);
    return std::real(E2);
}

double mband::SRO_Hubbard_Energy(NewAmiCalc::ext_vars ext, std::vector<double> momenta, int species, mband::params_param param) {
    const double mu_a = 1.0;
    const double mu_b = 1.1;
    const double t = 1.0;
    const double t_t = 0.1;
    const double tp = 0.8;
    const double tpp = 0.3;
    const double t_orb = param.t_orb;
    const double s = 1;
    const double SOC = param.SOC;

    const double kx = momenta[0];
    const double ky = momenta[1];

    const double cos_kx = cos(kx);
    const double cos_ky = cos(ky);
    const double sin_kx = sin(kx);
    const double sin_ky = sin(ky);

    const double A = -2 * (t * cos_kx + t_t * cos_ky) - mu_a; // dxz
    const double B = -2 * (t_t * cos_kx + t * cos_ky) - mu_a; // dyz
    const double C = -2 * tp * (cos_kx + cos_ky) - 4 * tpp * (cos_kx * cos_ky) - mu_b; // dxy
    const double g = -4 * t_orb * sin_kx * sin_ky;

   
    if (species == 1 || species == 2) {
            return mband::gamma(A, B, C, g, SOC, s);  
        } else if (species == 3 || species == 4) {
            return mband::beta(A, B, C, g, SOC, s); 
        } else if (species == 5 || species == 6) {    
            return mband::alpha(A, B, C, g, SOC, s); 
        }
    else {
        std::cerr << "Species number should be 1-6 for Tri-layer Hubbard problem" << std::endl;
        return 0.0;
    }
}