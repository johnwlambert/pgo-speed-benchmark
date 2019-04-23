//
// Created by John Lambert-Admin on 4/19/19.
//

#pragma once

#include <string>
#include <iostream>
#include <numeric>

#include "se2.h"
#include "pose_graph.h"
#include "manifold_constraints.h"
#include "timing_utils.h"

#include "sparse_lin_solver.h"


class GaussNewtonOptimizer
{
public:
    GaussNewtonOptimizer(PoseGraph2D & g, std::string dataset_name);
    void optimize();

    private:
        MatrixXd H_;
        VectorXd b_;
        VectorXd dx_;
        PoseGraph2D g_;
        std::string dataset_name_;
        static constexpr size_t max_iters_{100};

        // pre-allocate space for the error vector and Jacobians
        MatrixXd e_;
        MatrixXd A_; // Jacobian w.r.t. x1
        MatrixXd B_; // Jacobian w.r.t. x2

        void linearize();
        void solve_system();
        void linearize_pose_landmark_constraint(Vec3 x_i, Vec2 l, Vec2 z);
        void linearize_pose_pose_constraint(Vec3 v_i, Vec3 v_j, Vec3 z_ij);
};


