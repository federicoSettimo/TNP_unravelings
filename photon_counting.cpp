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

// Run it with all 3 labels to have all data for plots
//string label = "0"; // xi = 0
//string label = "+"; // xi = .1
string label = "-"; // xi = -1

// Parameters
int dim_max = 10, Nensemble = 1000;
double gamma = 1., Omega = 1., phi = .2, nbar = .5, dt = .001, tmax = 2.;
// Possible values of xi
double xi_0 = 0., xi_p = 0.01, xi_m = -1.;

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
    if (label == "0") xi = xi_0;
    else if (label == "+") xi = xi_p;
    else xi = xi_m;
    return gamma*(nbar + 1.)*a()*rho*adag() + gamma*nbar*adag()*rho*a() + xi*(nbar + 1.)*a()*rho*adag();
}

MatrixXcd Gamma (double t) {
    return gamma*(nbar + 1.)*adag()*a() + gamma*nbar*a()*adag();
}

MatrixXcd Gamma_L (double t) {
    double xi;
    if (label == "0") xi = xi_0;
    else if (label == "+") xi = xi_p;
    else xi = xi_m;
    return gamma*(nbar + 1.)*adag()*a() + gamma*nbar*a()*adag() + xi*(nbar + 1.)*adag()*a();
}

MatrixXcd L (const MatrixXcd &rho, double t) {
    return -I*comm(H, rho) + J(rho, t) - .5*anticomm(Gamma(t), rho);
}

double observable_0 (const MatrixXcd &rho) {return real(rho(0,0));}
double observable_1 (const MatrixXcd &rho) {return real(rho(1,1));}

double rand01 () {return rand()/(double)RAND_MAX;}
double max (double a, double b) {return a >= b ? a : b;}

// Chooses where the jump happens
VectorXcd jump (const MatrixXcd &R, double z) {
    ComplexEigenSolver<MatrixXcd> eigs;
    eigs.compute(R);
    VectorXcd eigval = eigs.eigenvalues();
    MatrixXcd eigvec = eigs.eigenvectors();

    double sum_previous_eigs = 0.;
    for (int j = 0; j < dim_max; ++j) {
        if (z >= sum_previous_eigs*dt && z < (sum_previous_eigs + real(eigval(j)))*dt)
            return eigvec.col(j)*exp(-arg(eigvec.col(j)[0]));
        sum_previous_eigs += real(eigval(j));
    }
    return eigvec.col(dim_max-1)*exp(-arg(eigvec.col(dim_max-1)[0]));
}

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

    int index_time = 0;
    for (double t = 0.; t < tmax; t += dt) {
        if (index_time % 100 == 0)
            cout << t << ", " << psi.size() << endl;
        index_time++;

        MatrixXcd K = H -.5*I*Gamma(t), rho_avg = MatrixXcd::Zero(dim_max, dim_max);
        int Nt = 0;

        //cout << t << ", " << psi.size() << endl;

        for (int i = 0; i < psi.size(); ++i) {
            // If destroyed do nothing
            if (!destroyed[i]) {
                rho_avg += projector(psi[i])/(double)Nensemble;
                Nt++;

                MatrixXcd R = J(projector(psi[i]), t);

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