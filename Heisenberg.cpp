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
static Eigen::Matrix2cd sigma_x {{1.,0.},{0.,-1.}};
static Eigen::Matrix2cd sigma_z {{0.,1.},{1.,0.}};
static Eigen::Matrix2cd sigma_p {{1.,-1.},{1.,-1.}};
static Eigen::Matrix2cd sigma_m {{1.,1.},{-1.,-1.}};
static Eigen::Matrix2cd id {{1,0},{0,1}};

static Eigen::Vector2cd plus_state {{1.,0.}};
static Eigen::Vector2cd minus_state {{0.,1.}};

// Parameters
double dt = .005, tmax = 2., threshold_neg = 1e-1;
int Nensemble = 10000, Ntraj = 6;

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

Matrix2cd Gamma_L (double t) {
    return gamma_m(t)*sigma_m*sigma_p + gamma_p(t)*sigma_p*sigma_m;
}

Matrix2cd L (const Matrix2cd &rho, double t) {
    return I*comm(H(t), rho) + J(rho, t) - .5*anticomm(Gamma(t), rho);
}

double observable_x (const Matrix2cd &rho) {return real((rho*sigma_x).trace());}
double observable_z (const Matrix2cd &rho) {return real((rho*sigma_z).trace());}

double rand01 () {return rand()/(double)RAND_MAX;}
double max (double a, double b) {return a >= b ? a : b;}

Vector2cd jump (const Matrix2cd &R, double z) {
    ComplexEigenSolver<Matrix2cd> eigs;
    eigs.compute(R);
    Vector2cd eigval = eigs.eigenvalues();
    Matrix2cd eigvec = eigs.eigenvectors();

    if (real(eigval[0]) < -threshold_neg || real(eigval[1]) < -threshold_neg) {
        cerr << "Error: negative eigenvalue" << endl;
        cerr << "\tEigenvalues: " << real(eigval[0]) << " " << real(eigval[1]) << endl;
    }
    if (z < real(eigval[0])*dt)
        return eigvec.col(0);
    else 
        return eigvec.col(1);
}

int main () {
    srand(time(NULL));

    Vector2cd psi0 = (plus_state + .3*minus_state).normalized();
    Matrix2cd rho = projector(psi0);

    ofstream out_rates("rates.txt"), out_tr("trace.txt"), out_obs_ex("obs.txt") , out_obs_avg("obs_avg.txt"), out_tmax("tmax.txt"), out_traj_x("traj_x.txt"), out_traj_z("traj_z.txt");
    out_tmax << tmax << endl << Ntraj;

    vector<Vector2cd> psi;
    vector<bool> destroyed;
    for (int i = 0; i < Nensemble; ++i) {
        psi.push_back(psi0);
        destroyed.push_back(false);
    }

    for (double t = 0.; t < tmax; t += dt) {
        Matrix2cd K0 = -H(t) -.5*I*Gamma(t), rho_avg = Matrix2cd::Zero();
        int Nt = 0;
        double gp = gamma_p(t), gm = gamma_m(t), ee = eps(t);

        cout << t << ", " << psi.size() << endl;

        for (int i = 0; i < psi.size(); ++i) {
            if (i < Ntraj) {
                if (!destroyed[i]) {
                    out_traj_x << observable_x(projector(psi[i])) << " ";
                    out_traj_z << observable_z(projector(psi[i])) << " ";
                }
                else { // Some value for destroyed trajectories
                    out_traj_x << -10. << " ";
                    out_traj_z << -10. << " ";
                }
            }
            // If destroyed do nothing
            if (!destroyed[i]) {
                rho_avg += projector(psi[i])/(double)Nensemble;
                Nt++;

                Vector2cd Phi;
                if (abs(psi[i][0]) > .99) // |+>
                    Phi = (2.*gp - 2.*gm)*minus_state - (gm + gp)*plus_state;
                else if (abs(psi[i][1]) > .99) // |->
                    Phi = (2.*gp - 2.*gm)*plus_state - (gm + gp)*minus_state;
                else {
                    double theta = arg(psi[i][1]), a = abs(psi[i][1]), sq = sqrt(1. - a*a),
                        beta = -2.*(gm - 2.*cos(theta)*a*sq*(gm-gp) + gp)/(2.*cos(2.*theta)*a);
                    complex<double> phip = (2.*((a*(-1 + pow(a,2))*exp(3.*I*theta)*(1. + exp(2.*I*theta))*(gm - gp) + 
                    sqrt(1 - pow(a,2))*exp(4.*I*theta)*(gm + gp))/
                  (1. + exp(4.*I*theta)) + 
                 a*(exp(I*theta)*(-gm + gp) + 
                    a*sqrt(1 - pow(a,2))*(1. + exp(2.*I*theta))*(gm + gp))))/pow(a,2);
                    complex<double> phim = (2*a*sqrt(1 - pow(a,2))*(1. + exp(2.*I*theta))*(gm - gp) - 
                    2.*exp(1.*I*theta)*(gm + gp))/(a*(1. + exp(4.*I*theta)));
                    Phi = phip*plus_state + phim*minus_state;
                }
                //Phi = 0.*plus_state;
                Matrix2cd K = K0 - .5*I*Phi*(psi[i].adjoint());
                Matrix2cd R = J(projector(psi[i]), t) + .5*(psi[i]*(Phi.adjoint()) + Phi*(psi[i].adjoint()));
                
                // Compute all probabilities
                double pj = real(R.trace())*dt,
                    p_d = max(0, real(((Gamma(t) - Gamma_L(t))*projector(psi[i])).trace()) * dt),
                    p_c = max(0, -real(((Gamma(t) - Gamma_L(t))*projector(psi[i])).trace()) * dt),
                    pdet = 1. - pj - p_d - p_c;
                double z = rand01();
                if (z < pj) // Jump
                    psi[i] = jump(R, z);
                if (z >= pj && z < p_c + pj) { // Create a copy
                    psi[i] = (id - I*K*dt)*psi[i];
                    psi.push_back(psi[i]);
                    destroyed.push_back(false);
                }
                else if (z >= p_c + pj && z < p_c + p_d + pj) // Destroy it
                    destroyed[i] = true;
                else // Det evo
                    psi[i] = (id - I*K*dt)*psi[i];
                psi[i] *= exp(-I*arg(psi[i][0]));
                psi[i].normalize();
            }
        }
        out_rates << gamma_m(t) << " " << gamma_p(t) << " " << eps(t) << endl;
        out_tr << real(rho.trace()) << " " << real(rho_avg.trace()) << " " << (double)Nt/(double)Nensemble << endl;
        out_obs_ex << observable_x(rho) << " " << observable_z(rho) << endl;
        out_obs_avg << observable_x(rho_avg) << " " << observable_z(rho_avg) << endl;
        out_traj_x << endl;
        out_traj_z << endl;

        rho += L(rho, t)*dt;
    }

    return 0;
}