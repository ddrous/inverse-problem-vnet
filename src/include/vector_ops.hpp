#ifndef INCLUDED_VECTOR_ARITHMETICS
#define INCLUDED_VECTOR_ARITHMETICS

#include <cmath>
#include <vector>
#include <iostream>
// #include <stdexcept>

/*****************************************
 * Pour les operations arithmetiques avec le type vecteur
 */

using namespace std;


// Check that the vectors have the same size
void check_size(vector<double> const& u, vector<double> const& v){
    if (u.size() != v.size() || u.size() != 2) {
            cout << "u = " << u.size() << "v" << v.size() << endl;
        throw std::string("ERREUR: Les tailles des vecteurs ne correspondent pas");
    }
}


double l2_norm(vector<double> const& u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}


vector<double> prod(double a, vector<double> const& u) {
    vector<double> product(u.size());
    for (int i = 0; i < u.size(); ++i) {
        product[i] = a * u[i];
    }
    return product;
}


vector<double> add(vector<double> const& u, vector<double> const& v) {
    check_size(u, v);

    vector<double> sum(u.size());
    for (int i = 0; i < u.size(); ++i) {
        sum[i] = u[i] + v[i];
    }
    return sum;
}


double dot(vector<double> const& u, vector<double> const& v) {
    check_size(u, v);
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
        accum += u[i] * v[i];
    }
    return accum;
}

#endif
