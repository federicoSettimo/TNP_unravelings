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

MatrixXcd comm (const MatrixXcd &A, const MatrixXcd &B) {return A*B-B*A;}
MatrixXcd anticomm (const MatrixXcd &A, const MatrixXcd &B) {return A*B+B*A;}
MatrixXcd projector (const VectorXcd &psi) {return psi*psi.adjoint();}

// Parameters
//string label = "0"; // xi = 0
//string label = "+"; // xi = .1
string label = "-"; // xi = -1
int dim_max = 2, Nensemble = 10000;
double gamma = 1., Omega = 1., phi = .2, nbar = .5, dt = .0005, tmax = 2.;

MatrixXd a () {
    MatrixXd aa = MatrixXd::Zero(dim_max, dim_max);
    for (int n = 1; n < dim_max; ++n) {
        aa(n - 1, n) = std::sqrt(n);
    }
    return aa;
}

MatrixXd adag () {
    return a().transpose();
}

MatrixXcd id = MatrixXcd::Identity(dim_max, dim_max);
MatrixXcd H = .5*Omega*(a()*exp(-2.*I*phi) + adag()*exp(2.*I*phi));

MatrixXcd J (const MatrixXcd &rho, double t) {
    double xi;
    if (label == "0") xi = 0.;
    else if (label == "+") xi = .01;
    else xi = -1.;
    return gamma*(nbar + 1.)*a()*rho*adag() + gamma*nbar*adag()*rho*a() + xi*(nbar + 1.)*a()*rho*adag();
}

MatrixXcd Gamma (double t) {
    return gamma*(nbar + 1.)*adag()*a() + gamma*nbar*a()*adag();
}

MatrixXcd L (const MatrixXcd &rho, double t) {
    return -I*comm(H, rho) + J(rho, t) - .5*anticomm(Gamma(t), rho);
}

double observable_0 (const MatrixXcd &rho) {return real(rho(0,0));}
double observable_1 (const MatrixXcd &rho) {return real(rho(1,1));}

double rand01 () {return rand()/(double)RAND_MAX;}

int main () {
    srand(time(NULL));

    VectorXd psi0 = VectorXd::Zero(dim_max);
    psi0(0) = 1.;
    psi0(1) = 1.;
    psi0.normalize();
    MatrixXcd rho = projector(psi0);

    ofstream out_tr("trace_"+label+".txt"), out_obs_ex("obs_"+label+".txt") , out_obs_avg("obs_avg_"+label+".txt"), out_tmax("tmax_"+label+".txt");
    out_tmax << tmax;

    vector<VectorXcd> psi;
    vector<bool> destroyed;
    for (int i = 0; i < Nensemble; ++i) {
        psi.push_back(psi0);
        destroyed.push_back(false);
    }

    for (double t = 0.; t < tmax; t += dt) {
        cout << t << ", " << psi.size() << endl;
        MatrixXcd K = H -.5*I*Gamma(t), rho_avg = MatrixXcd::Zero(dim_max, dim_max);
        int Nt = 0;

        //cout << t << ", " << psi.size() << endl;

        for (int i = 0; i < psi.size(); ++i) {
            // If destroyed do nothing
            if (!destroyed[i]) {
                rho_avg += projector(psi[i])/(double)Nensemble;
                Nt++;

                MatrixXcd R = J(projector(psi[i]), t);
                // Compute all probabiities
                double pdet = ((id - I*K*dt)*psi[i]).squaredNorm(), pj = real(R.trace())*dt;
                double p_create = pdet +pj - 1., p_destroy = 1. - pdet - pj;

                // First of all: should we destroy it?
                if (p_destroy > 0) {
                    double z = rand01();
                    if (z > 1. - p_destroy)
                        destroyed[i] = true;
                    else if (z < pj) {
                        ComplexEigenSolver<MatrixXcd> eigs;
                        eigs.compute(R);
                        VectorXcd eigval = eigs.eigenvalues();
                        MatrixXcd eigvec = eigs.eigenvectors();
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
                    z = rand01();
                    if (z < real(R.trace())*dt) {
                        ComplexEigenSolver<MatrixXcd> eigs;
                        eigs.compute(R);
                        VectorXcd eigval = eigs.eigenvalues();
                        MatrixXcd eigvec = eigs.eigenvectors();
                        if (z < real(eigval[0])*dt)
                            psi[i] = eigvec.col(0);
                        else psi[i] = eigvec.col(1);
                    }
                    else psi[i] = (id - I*K*dt)*psi[i];
                }
                psi[i].normalize();
            }
        }
        out_tr << real(rho.trace()) << " " << real(rho_avg.trace()) << " " << (double)Nt/(double)Nensemble << endl;
        out_obs_ex << observable_0(rho) << " " << observable_1(rho) << endl;
        out_obs_avg << observable_0(rho_avg) << " " << observable_1(rho_avg) << endl;

        rho += L(rho, t)*dt;
    }

    return 0;
}