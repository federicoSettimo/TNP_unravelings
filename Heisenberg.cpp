#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
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
static Eigen::Matrix2cd sigma_x {{0,1},{1,0}};
static Eigen::Matrix2cd sigma_y {{0,-I},{I,0}};
static Eigen::Matrix2cd sigma_z {{1,0},{0,-1}};
static Eigen::Matrix2cd sigma_p {{0,1},{0,0}};
static Eigen::Matrix2cd sigma_m {{0,0},{1,0}};
static Eigen::Matrix2cd id {{1,0},{0,1}};

static Eigen::Vector2cd ground_state {{0.,1.}};
static Eigen::Vector2cd excited_state {{1.,0.}};
static Eigen::Vector2cd plus_state {{1./sqrt(2.),1./sqrt(2.)}};
static Eigen::Vector2cd minus_state {{1./sqrt(2.),-1./sqrt(2.)}};

Matrix2cd comm (const Matrix2cd &A, const Matrix2cd &B) {return A*B-B*A;}
Matrix2cd anticomm (const Matrix2cd &A, const Matrix2cd &B) {return A*B+B*A;}
Matrix2cd projector (const Vector2cd &psi) {return psi*psi.adjoint();}

double eps (double t) {
    return 1. + .5*cos(2.*t);
}

Matrix2cd H (double t) {
    return eps(t)*sigma_x;
}

double gamma_p (double t) {
    return .5*(1. + .8*cos(t));
}

double gamma_m (double t) {
    return .5*(1. - .8*cos(t));
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

double observable_x (const Matrix2cd &rho) {return real((rho*sigma_x).trace());}
double observable_z (const Matrix2cd &rho) {return real((rho*sigma_z).trace());}

double rand01 () {return rand()/(double)RAND_MAX;}

int main () {
    srand(time(NULL));

    double dt = .005, tmax = 2.;
    int Nensemble = 10000;

    Vector2cd psi0 = (excited_state + .4*plus_state).normalized();
    Matrix2cd rho = projector(psi0);

    ofstream out_rates("rates.txt"), out_tr("trace.txt"), out_obs_ex("obs.txt") , out_obs_avg("obs_avg.txt"), out_tmax("tmax.txt");
    out_tmax << tmax;

    vector<Vector2cd> psi;
    vector<bool> destroyed;
    for (int i = 0; i < Nensemble; ++i) {
        psi.push_back(psi0);
        destroyed.push_back(false);
    }

    for (double t = 0.; t < tmax; t += dt) {
        Matrix2cd K0 = -H(t) -.5*I*Gamma(t), rho_avg = Matrix2cd::Zero();
        int Nt = 0;

        //cout << t << ", " << psi.size() << endl;

        for (int i = 0; i < psi.size(); ++i) {
            // If destroyed do nothing
            if (!destroyed[i]) {
                rho_avg += projector(psi[i])/(double)Nensemble;
                Nt++;

                //double a = real((psi[i])[1]), sq = sqrt(1.-a*a),
                double phi1 = (-3.*gamma_m(t) + gamma_p(t))/(4.*sqrt(2.)),
                    phi0 = (gamma_m(t) - 3.*gamma_p(t))/(4.*sqrt(2.));
                Vector2cd Phi = phi0*ground_state + phi1*excited_state;
                Phi = 0.*ground_state;
                Matrix2cd K = K0 - .5*I*Phi*(psi[i].adjoint());
                //cout << Phi << endl;
                Matrix2cd R = J(projector(psi[i]), t) + .5*(psi[i]*(Phi.adjoint()) + Phi*(psi[i].adjoint()));
                // Compute all probabiities
                double pdet = ((id - I*K*dt)*psi[i]).squaredNorm(), pj = real(R.trace())*dt;
                double p_create = pdet +pj - 1., p_destroy = 1. - pdet - pj;

                // First of all: should we destroy it?
                if (p_destroy > 0) {
                    double z = rand01();
                    if (z > 1. - p_destroy)
                        destroyed[i] = true;
                    else if (z < pj) {
                        ComplexEigenSolver<Matrix2cd> eigs;
                        eigs.compute(R);
                        Vector2cd eigval = eigs.eigenvalues();
                        Matrix2cd eigvec = eigs.eigenvectors();
                        if (z < real(eigval[0])*dt)
                            psi[i] = eigvec.col(0);
                        else psi[i] = eigvec.col(1);
                    }
                    else psi[i] = (id - I*K*dt)*psi[i];
                }
                else { // Or should we create a new copy?
                    double z = rand01();
                    if (z > 1. - p_create) {
                        psi.push_back(psi[i]);
                        destroyed.push_back(false);
                    }
                    // Now renormalize the probabilities
                    //double normalizing_factor = pdet + pj + pj;
                    //R /= normalizing_factor;
                    //pdet /= normalizing_factor;
                    z = rand01();
                    if (z < real(R.trace())*dt) {
                        ComplexEigenSolver<Matrix2cd> eigs;
                        eigs.compute(R);
                        Vector2cd eigval = eigs.eigenvalues();
                        Matrix2cd eigvec = eigs.eigenvectors();
                        if (z < real(eigval[0])*dt)
                            psi[i] = eigvec.col(0);
                        else psi[i] = eigvec.col(1);
                    }
                    else psi[i] = (id - I*K*dt)*psi[i];
                }
                psi[i] *= exp(-I*arg(psi[i][1]));
                psi[i].normalize();
            }
        }
        out_rates << gamma_m(t) << " " << gamma_p(t) << " " << eps(t) << endl;
        out_tr << real(rho.trace()) << " " << real(rho_avg.trace()) << " " << (double)Nt/(double)Nensemble << endl;
        out_obs_ex << observable_x(rho) << " " << observable_z(rho) << endl;
        out_obs_avg << observable_x(rho_avg) << " " << observable_z(rho_avg) << endl;

        rho += L(rho, t)*dt;
    }

    return 0;
}