// Checks whether the phase covarian dynamics in the Heisenberg picture is CP divisible in the Schrodinger picture
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <unsupported/Eigen/KroneckerProduct>
#include <complex>
#include <vector>
#include <cstring>
#include <fstream>
#include <ctime>
#include <cstdlib>
#include <iostream>

using namespace std;
using namespace Eigen;

static complex<double> I(0,1), one(1,0), zero(0,0);
static Matrix2cd sigma_x {{1.,0.},{0.,-1.}}, sigma_y {{0,-I},{I,0}}, sigma_z {{0.,1.},{1.,0.}}, sigma_p {{1.,-1.},{1.,-1.}}, sigma_m {{1.,1.},{-1.,-1.}}, id {{1,0},{0,1}};

static Vector2cd plus_state {{1.,0.}}, minus_state {{0.,1.}};

// Parameters
double dt = .005, tmax = 2., threshold_neg = 1e-1;

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B-B*A;}
Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B+B*A;}
Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

double eps (double t) {
    return 4.*(1. + .5*cos(2.*t));
}

Matrix2cd H (double t) {
    return eps(t)*sigma_x;
}

double gamma_p (double t) {
    return .1*(1. + .8*cos(t));
}

double gamma_m (double t) {
    return .1*(1. - .8*cos(t));
}

Matrix2cd J (const Matrix2cd &rho, double t) {
  return gamma_p(t)*sigma_m*rho*sigma_p + gamma_m(t)*sigma_p*rho*sigma_m;
}

Matrix2cd Gamma (double t) {
    return gamma_p(t)*sigma_m*sigma_p + gamma_m(t)*sigma_p*sigma_m;
}

Matrix2cd L (const Matrix2cd &rho, double t) {
    return I*comm(H(t), rho) + J(rho, t) - .5*anticomm(Gamma(t), rho);
}

int main () {
    vector<Matrix2cd> sigmas(4), sigmas_t(4);
    sigmas[0] = id;
    sigmas[1] = sigma_x;
    sigmas[2] = sigma_y;
    sigmas[3] = sigma_z;
    sigmas_t[0] = id;
    sigmas_t[1] = sigma_x;
    sigmas_t[2] = sigma_y;
    sigmas_t[3] = sigma_z;

    ofstream out_tmax("tmax.txt"), out("rates_H.txt");
    out_tmax << tmax;

    // Computing the Bloch matrix \Lambda = tr[\sigma_i \Phi_t[\sigma_j]]
    Matrix4d Lambda = Matrix4d::Zero(4,4);
    for (double t = 0.; t < tmax; t += dt) {
        Matrix4d Lambda_old = Lambda, Lambda_old_S = Lambda_old.transpose();
        for (int j = 0; j < 4; ++j) {
            sigmas_t[j] += dt*L(sigmas_t[j], t);
            for (int i = 0; i < 4; ++i) {
                Lambda(i,j) = .5*(sigmas[i]*sigmas_t[j]).trace().real();
            }
        }

        if (t != 0.) {
        //if (true) {
            // Computing the infinitesimal time evolution 
            Matrix4d Lambda_S = Lambda.transpose(), dLambda_dt_S = (Lambda_S - Lambda_old_S)/dt, Lambda_inv_S = Lambda_old_S.inverse();
            Matrix4d infinitesimal_S = Lambda_S*Lambda_inv_S;

            // Computing the images of the Pauli matrices in Schrodinger picture
            // under the infinitesimal evolution
            Matrix2cd id_t_S = .5*(infinitesimal_S(0,0)*id + infinitesimal_S(1,0)*sigma_x + infinitesimal_S(2,0)*sigma_y + infinitesimal_S(3,0)*sigma_z),
                        sigma_x_t_S = .5*(infinitesimal_S(0,1)*id + infinitesimal_S(1,1)*sigma_x + infinitesimal_S(2,1)*sigma_y + infinitesimal_S(3,1)*sigma_z),
                        sigma_y_t_S = .5*(infinitesimal_S(0,2)*id + infinitesimal_S(1,2)*sigma_x + infinitesimal_S(2,2)*sigma_y + infinitesimal_S(3,2)*sigma_z),
                        sigma_z_t_S = .5*(infinitesimal_S(0,3)*id + infinitesimal_S(1,3)*sigma_x + infinitesimal_S(2,3)*sigma_y + infinitesimal_S(3,3)*sigma_z);
            // From here, the images of |i><j|
            Matrix2cd E00_t = .5*(id_t_S + sigma_z_t_S), E01_t = .5*(sigma_x_t_S - I*sigma_y_t_S),
                        E10_t = .5*(sigma_x_t_S + I*sigma_y_t_S), E11_t = .5*(id_t_S - sigma_z_t_S);
            Matrix2cd E00 = .5*(id + sigma_z), E01 = .5*(sigma_x - I*sigma_y),
                        E10 = .5*(sigma_x + I*sigma_y), E11 = .5*(id - sigma_z);

            // Using the infinitesimal Schrodinger evolution to get the Choi state
            Matrix4cd Choi = KroneckerProduct(E00_t, E00) + KroneckerProduct(E01_t, E01) +
                                KroneckerProduct(E10_t, E10) + KroneckerProduct(E11_t, E11);
            SelfAdjointEigenSolver<Matrix4cd> es(Choi);
            Vector4d eigvals = es.eigenvalues();
            out << eigvals[0] << " " << eigvals[1] << " " << eigvals[2] << " " << eigvals[3] << " ";

            // Now for P divisibility: checking the maximum value of Bloch radius
            double max_r = 0.;
            Matrix3d R = infinitesimal_S.block<3,3>(1,1);
            Vector3d v = {infinitesimal_S(1,0), infinitesimal_S(2,0), infinitesimal_S(3,0)};
            for (double theta = 0.; theta < 2.*M_PI; theta += M_PI/20.) {
                for (double phi = 0.; phi < M_PI; phi += M_PI/20.) {
                    Vector3d n = {sin(phi)*cos(theta), sin(phi)*sin(theta), cos(phi)};
                    double r = (R*n + v).norm();
                    if (r > max_r) max_r = r;
                }
            }

            out << max_r << endl;
        }
    }

    return 0;
}