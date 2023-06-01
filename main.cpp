#include <iostream>
#include <cassert>
#include <vector>
#include <cmath>
#include "matrix.hpp"

// Y = Ax^3 + Cx + D

// Sum i from 0 to n-1 (xi^powx * yi^powy)
double coef(const std::vector<double> &x, double powx, const std::vector<double> &y, double powy) {
    double sum = 0;

    for (unsigned i = 0; i < x.size(); ++i) {
        sum += std::pow(y[i], powy) * std::pow(x[i], powx);
    }

    return sum;
}

std::vector<double> solve(const Matrix<double> &mat, const std::vector<double> &y) {
    assert(mat.N == y.size());

    double d = mat.det();

    assert(d != 0);

    std::vector<double> x;
    x.resize(y.size());

    for (unsigned j = 0; j < y.size(); ++j) {
        Matrix<double> matj = mat;

        for (int i = 0; i < y.size(); ++i) {
            matj.at(i, j) = y[i];
        }

        double dj = matj.det();

        x[j] = dj/d;
    }

    return x;
}

double avg(const std::vector<double> &arr) {
    assert(!arr.empty());

    double sum = 0;

    for (unsigned i = 0; i < arr.size(); ++i) {
        sum += arr[i];
    }

    return sum/arr.size();
}

double quadraticErr(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &a) {
    assert(a.size() == 3);

    double sum = 0;

    for (unsigned i = 0; i < x.size(); ++i) {
        double f = a[0] + a[1]*x[i] + a[2]*x[i]*x[i]*x[i];

        sum += (f - y[i])*(f  - y[i]);
    }

    return sum/x.size();
}

double avgErrLinear(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &a) {
    assert(a.size() == 2);

    double sum = 0;

    for (unsigned i = 0; i < x.size(); ++i) {
        double f = a[0] + a[1]*x[i];

        sum += (f - y[i]);
    }

    return sum/x.size();
}

double avgErrFunc(const std::vector<double> &x, const std::vector<double> &y, const std::vector<double> &a) {
    assert(a.size() == 3);

    double sum = 0;

    for (unsigned i = 0; i < x.size(); ++i) {
        double f = a[0] + a[1]*x[i] + a[2]*x[i]*x[i]*x[i];

        sum += (f - y[i]);
    }

    return sum/x.size();
}

double correlationCoef(const std::vector<double> &x, const std::vector<double> &y) {
    assert(x.size() == y.size());

    double mx = avg(x);
    double my = avg(y);

    double u = 0;

    double ux = 0, uy = 0;

    for (unsigned i = 0; i < x.size(); ++i) {
        u += (x[i] - mx)*(y[i] - my);
        ux += (x[i] - mx)*(x[i] - mx);
        uy += (y[i] - my)*(y[i] - my);
    }

    return u/std::sqrt(ux*uy);
}

int main() {
    unsigned n;
    std::vector<double> x, y;

    std::cin>>n;

    x.resize(n);
    for (unsigned i = 0; i < n; ++i) {
        std::cin>>x[i];
    }

    y.resize(n);
    for (unsigned i = 0; i < n; ++i) {
        std::cin>>y[i];
    }

    {
        std::vector<double> right;
        Matrix<double> system(2);

        system.at(0, 0) = n;
        system.at(0, 1) = coef(x, 1, y, 0);

        system.at(1, 0) = coef(x, 1, y, 0);
        system.at(1, 1) = coef(x, 2, y, 0);

        right.push_back(coef(x, 0, y, 1));
        right.push_back(coef(x, 1, y, 1));

        std::vector<double> a = solve(system, right);

        std::cout<<"Regression line: "<<a[0]<<" + "<<a[1]<<"*x"<<std::endl;

        std::cout<<"Correlation coeficent: "<<correlationCoef(x, y)<<std::endl;

        std::cout<<"Average error (linear): "<<avgErrLinear(x, y, a)<<std::endl;
    }

    {
        std::vector<double> right;
        Matrix<double> system(3);

        system.at(0, 0) = n;
        system.at(0, 1) = coef(x, 1, y, 0);
        system.at(0, 2) = coef(x, 3, y, 0);

        system.at(1, 0) = coef(x, 1, y, 0);
        system.at(1, 1) = coef(x, 2, y, 0);
        system.at(1, 2) = coef(x, 4, y, 0);

        system.at(2, 0) = coef(x, 3, y, 0);
        system.at(2, 1) = coef(x, 4, y, 0);
        system.at(2, 2) = coef(x, 6, y, 0);

        right.push_back(coef(x, 0, y, 1));
        right.push_back(coef(x, 1, y, 1));
        right.push_back(coef(x, 3, y, 1));

        std::vector<double> a = solve(system, right);

        std::cout<<"Function: ";
        std::cout<<a[0]<<" + "<<a[1]<<"*x + "<<a[2]<<"*x**3"<<std::endl;

        std::cout<<"Quadratic error: "<<quadraticErr(x, y, a)<<std::endl;

        std::cout<<"Average error (function): "<<avgErrFunc(x, y, a)<<std::endl;
    }


    return EXIT_SUCCESS;
}
