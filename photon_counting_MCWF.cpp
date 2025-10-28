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
double gamma = 1., Omega = 1., phi = .2, nbar = .5, dt = .01, tmax = 1., dz = 0.01;
int dim_max = 10, Nensemble = 10000, Nsmooth = (int)tmax/dt/50;

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

MatrixXcd J (const MatrixXcd &rho, double t, double z = 0.) {
    return gamma*(nbar + 1.)*a()*rho*adag() + gamma*nbar*adag()*rho*a() + z*(nbar + 1.)*a()*rho*adag();
}

MatrixXcd J_extra (const MatrixXcd &rho, double t) { // J in the new MEs
    return (nbar + 1.)*a()*rho*adag();
}

MatrixXcd Gamma (double t) {
    return gamma*(nbar + 1.)*adag()*a() + gamma*nbar*a()*adag();
}

MatrixXcd Gamma_L (double t, double z) {
    return gamma*(nbar + 1.)*adag()*a() + gamma*nbar*a()*adag() + z*(nbar + 1.)*adag()*a();
}

MatrixXcd L (const MatrixXcd &rho, double t, double z = 0.) {
    return -I*comm(H, rho) + J(rho, t, z) - .5*anticomm(Gamma(t), rho);
}

double observable_0 (const MatrixXcd &rho) {return real(rho(0,0));}
double observable_1 (const MatrixXcd &rho) {return real(rho(1,1));}

double rand01 () {return rand()/(double)RAND_MAX;}
double max (double a, double b) {return a >= b ? a : b;}

vector<MatrixXcd> unravel (const VectorXcd &psi0, int k, const vector<MatrixXcd> &rhoKm1);
vector<double> evolve_ex (const VectorXcd &psi0, const double z, int index);

vector<MatrixXcd> smooth (const vector<MatrixXcd> &rho); // Smoothens the MC evolution

void unravel_tr (const VectorXcd &psi0, const double z, int index);

int main () {
    srand(42);

    VectorXd psi0 = VectorXd::Zero(dim_max);
    psi0(0) = 1.;
    psi0(1) = 1.;
    psi0.normalize();

    vector<double> rho0_ex, rho1_ex, rho2_ex, rho3_ex, rho4_ex;

    vector<double> zs {-2.*dz, -dz, 0., dz, 2.*dz};

    ofstream out_tmax("tmax.txt"), out_mom("moments.txt"), out_mom_ex("moments_ex.txt");
    out_tmax << tmax;
    for (double z : zs) out_tmax << "\n" << z;
    out_tmax.close();

    // Getting the exact solution + derivatives
    cout << "Unraveling the trace and getting the exact solutions...\n";
    int i = 0;
    cout << "\tz = " << zs[i] << "\n";
    unravel_tr(psi0, zs[i], i);
    rho0_ex = evolve_ex(psi0, zs[i], i); ++i;
    cout << "\tz = " << zs[i] << "\n";
    unravel_tr(psi0, zs[i], i);
    rho1_ex = evolve_ex(psi0, zs[i], i); ++i;
    cout << "\tz = " << zs[i] << "\n";
    unravel_tr(psi0, zs[i], i);
    rho2_ex = evolve_ex(psi0, zs[i], i); ++i;
    cout << "\tz = " << zs[i] << "\n";
    unravel_tr(psi0, zs[i], i);
    rho3_ex = evolve_ex(psi0, zs[i], i); ++i;
    cout << "\tz = " << zs[i] << "\n";
    unravel_tr(psi0, zs[i], i);
    rho4_ex = evolve_ex(psi0, zs[i], i);

    vector<MatrixXcd> rho0, rho1, rho2, rho3, rho4;
    // Compute the moments etc...
    cout << "\nUnraveling the moments...\n";
    cout << "\tUnraveling tau_0...\n";
    rho0 = smooth(unravel(psi0, 0, rho0));
    cout << "\tUnraveling tau_1...\n";
    rho1 = smooth(unravel(psi0, 1, rho0));
    cout << "\tUnraveling tau_2...\n";
    rho2 = smooth(unravel(psi0, 2, rho1));
    cout << "\tUnraveling tau_3...\n";
    rho3 = smooth(unravel(psi0, 3, rho2));
    cout << "\tUnraveling tau_4...\n";
    rho4 = smooth(unravel(psi0, 4, rho3));


    for (int i = 0; i < rho0.size(); ++i) {
        double m1_ex = -.5*((rho1_ex[i] - rho2_ex[i])/dz + (rho0_ex[i] - rho1_ex[i])/dz);
        double m2_ex = ((rho0_ex[i] - 2.*rho1_ex[i] + rho2_ex[i])/(dz*dz)
            + (rho1_ex[i] - 2.*rho2_ex[i] + rho3_ex[i])/(dz*dz)
            + (rho2_ex[i] - 2.*rho3_ex[i] + rho4_ex[i])/(dz*dz))/3.;
        double m3_ex = -(rho0_ex[i]-2*rho1_ex[i]+2*rho3_ex[i]-rho4_ex[i])/(2.*dz*dz*dz);
        double m4_ex = (rho0_ex[i]-4*rho1_ex[i]+6*rho2_ex[i]-4*rho3_ex[i]+rho4_ex[i])/(dz*dz*dz*dz);

        out_mom_ex << m1_ex << " " << m2_ex << " " << m3_ex << " " << m4_ex << endl;
        out_mom << real(rho1[i].trace()) << " "
            << real(rho2[i].trace()) << " "
            << real(rho3[i].trace()) << " "
            << real(rho4[i].trace()) << endl;
    }

    return 0;
}

vector<double> evolve_ex (const VectorXcd &psi0, const double z, int index) {
    vector<double> m;
    MatrixXcd rho = projector(psi0);
    ofstream out_tr("trace_ex_"+to_string(index)+".txt");
    int it = 0;
    for (double t = 0.; t < tmax; t += dt) {
        rho += L(rho, t, z) * dt;
        double trr = real(rho.trace());
        m.push_back(trr);
        it++;
        out_tr << 1. - trr << endl;
    }
    return m;
}

vector<MatrixXcd> smooth (const vector<MatrixXcd> &rho) {
    int n = rho.size(), half = Nsmooth / 2;
    vector<MatrixXcd> rho_smooth;

    for (int i = 0; i < n; ++i) {
        int start = max(0, i - half);
        int end   = min(n - 1, i + half);
        MatrixXcd sum = MatrixXcd::Zero(dim_max, dim_max);
        for (int j = start; j <= end; ++j)
            sum += rho[j];
        rho_smooth.push_back( sum / (end - start + 1) );
    }

    return rho_smooth;
}

vector<MatrixXcd> unravel (const VectorXcd &psi0, int k, const vector<MatrixXcd> &rhoKm1) {
    vector<VectorXcd> psi;
    vector<bool> isZero;
    for (int i = 0; i < Nensemble; ++i) {
        psi.push_back(psi0);
        isZero.push_back(k>=1); // If k >= 1, the initial condition is actually a zero vector 
    }

    vector<MatrixXcd> m;
    int it = 0;
    for (double t = 0.; t < tmax; t += dt) {
        MatrixXcd K = H -.5*I*Gamma(t),
            rho_avg = MatrixXcd::Zero(dim_max, dim_max),
            At = MatrixXcd::Zero(dim_max, dim_max);
        if (k > 0)
            At = k * J_extra(rhoKm1[it], t);
        
        // Getting the eigendecomposition of A for the extra jumps
        ComplexEigenSolver<MatrixXcd> eigs;
        eigs.compute(At);
        VectorXcd eigval = eigs.eigenvalues();
        MatrixXcd eigvec = eigs.eigenvectors();

        for (int i = 0; i < psi.size(); ++i) {
            double z = rand01();
            // If it's zero: easy, just create copies
            if (isZero[i]) {
                bool stop = false;
                for (int j = 0; j < dim_max && !stop; ++j) {
                    double ll = real(eigval(j)) * dt;
                    if (z < ll) {
                        stop = true;
                        isZero[i] = false;
                        psi[i] = eigvec.col(j).normalized();
                    }
                    z -= ll;
                }
            }
            else {
                rho_avg += projector(psi[i])/(double)Nensemble;
                double pj_m = gamma*(nbar + 1.)*(a()*psi[i]).squaredNorm()*dt,
                        pj_p = gamma*nbar*(adag()*psi[i]).squaredNorm()*dt,
                        p_copy = real(At.trace()) * dt;
                if (z < p_copy) {
                    if (i < Nensemble) { // Is it right? Maybe
                        for (int j = 0; j < dim_max; ++j) {
                            double ll = real(eigval(j)) * dt;
                            if (z < ll) {
                                psi.push_back(eigvec.col(j).normalized());
                                isZero.push_back(false);
                                break;
                            }
                            z -= ll;
                        }
                    }
                }
                else if (z < p_copy + pj_m) // Jump minus
                    psi[i] = a()*psi[i];
                
                else if (z < p_copy + pj_m + pj_p) // Jump plus
                    psi[i] = adag()*psi[i];
                else
                    psi[i] = (id - I*K*dt)*psi[i];
                psi[i].normalize();
            }
        }
        m.push_back(rho_avg);
        it++;
    }
    return m;
}

void unravel_tr (const VectorXcd &psi0, const double z, int index) {
    vector<VectorXcd> psi;
    vector<bool> destroyed;
    for (int i = 0; i < Nensemble; ++i) {
        psi.push_back(psi0);
        destroyed.push_back(false);
    }

    ofstream out_tr("trace_"+to_string(index)+".txt");

    int it = 0;
    for (double t = 0.; t < tmax; t += dt) {
        MatrixXcd K = H -.5*I*Gamma(t);
        MatrixXcd rho_avg = MatrixXcd::Zero(dim_max, dim_max);

        for (int i = 0; i < psi.size(); ++i) {
            // If destroyed do nothing
            if (!destroyed[i]) {
                rho_avg += projector(psi[i])/(double)Nensemble;

                // Compute all probabilities
                double pj_m = gamma*(nbar + 1.)*(1. + z)*(a()*psi[i]).squaredNorm()*dt,
                    pj_p = gamma*nbar*(adag()*psi[i]).squaredNorm()*dt,
                    pj = pj_m + pj_p,
                    p_d = max(0, real(((Gamma(t) - Gamma_L(t,z))*projector(psi[i])).trace()) * dt),
                    p_c = max(0, -real(((Gamma(t) - Gamma_L(t,z))*projector(psi[i])).trace()) * dt),
                    pdet = 1. - pj_m - pj_p - p_d - p_c;
                double z = rand01();
                if (z < pj) {// Jump
                    if (z < pj_m) psi[i] = a()*psi[i];
                    else psi[i] = adag()*psi[i];
                }
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
        double trr = real(rho_avg.trace());
        out_tr << 1. - trr << endl;

        it++;
    }
    return;
}