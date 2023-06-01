#ifndef MATRIX_HPP
#define MATRIX_HPP

#include <iostream>
#include <iomanip>
#include <cstring>
#include <cassert>

template<class T>
class Matrix {
    T *data;

public:
    const unsigned N;
    const unsigned M;

    Matrix(unsigned n, unsigned m): N(n), M(m) {
        data = new T[n*m];
    }

    Matrix(unsigned n): Matrix(n, n) {}

    Matrix(const Matrix<T> &other): N(other.N), M(other.M) {
        data = new T[N*M];

        memcpy(data, other.data, N*M*sizeof(T));
    }

    Matrix<T> operator+(const Matrix<T> &other) const {
        assert(other.N == N && other.M == M);

        Matrix<T> result(N, M);

        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                result.at(i, j) = at(i, j) + other.at(i, j);
            }
        }

        return result;
    }

    Matrix<T> operator*(const Matrix<T> &other) const {
        assert(M == other.N);

        Matrix<T> result(N, other.M);

        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < other.M; ++j) {
                result.at(i, j) = 0;

                for (unsigned k = 0; k < M; ++k) {
                    result.at(i, j) += at(i, k) * other.at(k, j);
                }
            }
        }

        return result;
    }

    Matrix<T> operator*(T num) const {
        Matrix<T> result(N, M);

        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                result.at(i, j) = at(i, j)*num;
            }
        }

        return result;
    }

    Matrix<T> operator/(T num) const {
        Matrix<T> result(N, M);

        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                result.at(i, j) = at(i, j)/num;
            }
        }

        return result;
    }

    Matrix<T> trans() const {
        assert(M == N);

        Matrix<T> t(M, N);

        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j <= i; ++j) {
                t.at(j, i) = at(i, j);
                t.at(i, j) = at(j, i);
            }
        }

        return t;
    }

    T &at(unsigned i, unsigned j) {
        return data[i*M + j];
    }

    const T &at(unsigned i, unsigned j) const {
        return data[i*M + j];
    }

    const T det() const {
        assert(N == M);

        if (N == 0) { return 0; }
        if (N == 1) { return *data; }
        if (N == 2) { return at(0, 0)*at(1, 1) - at(1, 0)*at(0, 1); }

        T result = 0;

        for (unsigned j = 0; j < N; ++j) {
            Matrix<T> minor(N-1, N-1);

            for (unsigned mi = 0, copyi = 1; mi < minor.N; ++mi, ++copyi) {
                for (unsigned mj = 0, copyj = 0; mj < minor.M; ++mj, ++copyj) {
                    if (copyj == j) ++copyj;

                    minor.at(mi, mj) = at(copyi, copyj);
                }
            }

            result += at(0, j)*minor.det()*(j % 2 == 0 ? T(1) : T(-1));
        }

        return result;
    }

    Matrix<T> inverse() const {
        assert(N == M);

        Matrix<T> result(N);

        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < N; ++j) {
                Matrix<T> minor(N-1, N-1);

                for (unsigned mi = 0, copyi = 0; mi < minor.N; ++mi, ++copyi) {
                    if (copyi == i) ++copyi;

                    for (unsigned mj = 0, copyj = 0; mj < minor.M; ++mj, ++copyj) {
                        if (copyj == j) ++copyj;

                        minor.at(mi, mj) = at(copyi, copyj);
                    }
                }

                result.at(i, j) = minor.det() * ((i+j)%2 == 0 ? T(1) : T(-1));
            }
        }

        return result.trans() / det();
    }

    void readFromStream(std::istream &in) {
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                in>>data[i*M + j];
            }
        }
    }

    void writeToStream(std::ostream &out) const {
        for (unsigned i = 0; i < N; ++i) {
            for (unsigned j = 0; j < M; ++j) {
                out<<std::setw(12)<<data[i*M + j];
            }

            out<<std::endl;
        }
    }

    ~Matrix() {
        delete[] data;
    }
};

#endif
