//
// Created by John Lambert-Admin on 4/22/19.
//

#pragma once

#include <string>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>

#include "csparse_utils.h"
#include <Eigen/Dense>

using VectorXd = Eigen::VectorXd; // equiv. to Matrix<double,Dynamic,1>
using MatrixXd = Eigen::MatrixXd; // equiv. to Matrix<double,Dynamic,Dynamic>

class SparseLinSolver{

    public:
        SparseLinSolver(MatrixXd & A, VectorXd & b, VectorXd & x) : A_(A), b_(b), x_(x)
        {}

    void solve_cplusplus(std::string solver_type)
    {
            size_t n = A_.cols();
            size_t m = A_.rows();
            size_t idx_1d = 0;
            double *A_arr = (double *) malloc( sizeof(double) *m*n);
            double *b_arr = (double *) malloc( sizeof(double) *m);
            double *x_arr = (double *) malloc( sizeof(double) *m);
            size_t i = 0;
            size_t j = 0;

            for (i = 0; i < A_.rows(); i++)
            {
                    for (j = 0; j < A_.cols() ; j++)
                    {
                        // n entries in each row, plus up to current column
                        idx_1d = i*n + j;
                        A_arr[idx_1d] = A_(i,j);
                    }
                    b_arr[i] = b_(i,0);
            }

            solve_c(A_arr, b_arr, x_arr, A_.rows(), A_.cols() );

            free(A_arr);
            free(b_arr);
            free(x_arr);
    }


private:
    MatrixXd A_;
    VectorXd b_;
    VectorXd x_;
};




