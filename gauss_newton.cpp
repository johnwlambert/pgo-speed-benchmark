//
// Created by John Lambert-Admin on 4/20/19.
//

#include "gauss_newton.h"
#include "sparse_lin_solver.h"

/*
 * Constructor -- initialize the graph object.
 * Allocate enough space for dense H and b.
 */
SparseGaussNewtonOptimizer::SparseGaussNewtonOptimizer(PoseGraph2D & g, std::string dataset_name):
                            g_(g), dataset_name_(dataset_name)
{
    // could allocate the sparse H and the vector b here
    // get number of rows -- x already accounts for 3 values per vertex
    size_t state_vec_len = g.x_.rows();
    H_.resize(state_vec_len,state_vec_len);
    b_.resize(state_vec_len,1);
    dx_.resize(state_vec_len,1);

    bool analyze_nnz = false;
    if (analyze_nnz==true)
    {
        //        nnz = nnz_of_graph(g)
        //        print('H has shape ', H.shape, ' #elements in H is ', H.shape[0]*H.shape[1] ,' with ', nnz,' nnz')
        //        nnz_entry_percent = nnz / (H.shape[0]*H.shape[1]) * 100
        //        print('NNZ % = ', nnz_entry_percent)
        //        print('Not NNZ % = ', 100. - nnz_entry_percent)
    }
}


/*
 * Run the iterations for Gauss Newton algorithm.
 * Recall that g contains pose variables in a state vector "x" of shape (K,1).
 * We add delta x (aka "dx") to the state vector.
 */
void SparseGaussNewtonOptimizer::optimize()
{
    bool visualize=false;
    //// plot the initial state of the graph
    //plot_graph(fig, g, 0)

    //print('Initial error %f\n', compute_global_error(g))
    // maximum allowed dx
    double EPSILON = std::pow(10.0,-4);
    double error = 0;

    // carry out the iterations
    for (size_t iter=0; iter<max_iters_; iter++){
        std::cout << "Performing iteration " << iter << std::endl;

        linearize();
        solve_system();

        // Apply the solution to the state vector g.x
        g_.x_ += dx_;

        g_ = normalize_angles(g_);
        if (visualize==true){
            // plot the current state of the graph
            //plot_graph(fig, g, iter)
            // err = compute_global_error(g)
        }
        std::cout << "Current error = " << error;
        // termination criterion
        if (dx_.norm() < std::pow(1.0,-10) ){
            break;
        }
    }
    std::cout << "Final error =" << error << std::endl;
}


/*
 * We solve the linear system, solution stored in dx_
 * Instead of inverting H explicitly, we use a sparse solver.
 * Compute dx_, a matrix of shape (K,1) representing changed pose variables("delta x")
 */
void SparseGaussNewtonOptimizer::solve_system()
{
    std::cout << "\tSystem size: " << H_.rows() << " x " << H_.cols() << std::endl;
    std::cout << "\tSolving system (may take some time) ..." << std::endl;
    std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();


    SparseLinSolver sls(H_,b_,dx_);
    std::string solver_type = "cholesky";
    // solve for dx_ with sparse matrices
    sls.solve_cplusplus (solver_type);

    double duration = compute_elapsed_wall_clock_time(start_time, std::chrono::steady_clock::now() );
    std::cout << "Linear system solve complete -- took " << duration << " sec.";
}



/*
 * Each constraint is linearized and added to the Hessian.
 * Compute the addend term to H and b for each of our constraints.
        e the error vector
        A Jacobian wrt x1
        B Jacobian wrt x2
 */
void SparseGaussNewtonOptimizer::linearize()
{
    std::cout << "\tLinearize and build system" << std::endl;
    for (auto & edge : g_.edges_)
    {
        size_t i_v_id = edge.from_v_id_;
        size_t j_v_id = edge.to_v_id_;

        size_t i = g_.vertex_map_[i_v_id]->x_offset_idx_;
        size_t j = g_.vertex_map_[j_v_id]->x_offset_idx_;

        // height of Jacobian blocks
        size_t i_dim = 3; // Jacobian block w.r.t. i always 3x3
        size_t j_dim = 0;

        // pose-pose constraint
        if (edge.edge_type_=="P") {

            // edge.fromIdx and edge.toIdx describe the location of
            // the first element of the pose in the state vector
            j_dim = 3; // Jacobian block w.r.t. j 3x3 here

            // get (3x1) column vector out for 1st and 2nd robot poses
            Vec3 x1 = g_.x_.block(i, 0, 3, 1);
            Vec3 x2 = g_.x_.block(j, 0, 3, 1);

            // Computing the error and the Jacobians, using edge measurement "z"
            // e the error vector. A Jacobian wrt x1. B Jacobian wrt x2
            linearize_pose_pose_constraint(x1, x2, edge.measurement_);
            // Update H matrix and vector b, using Omega (the information matrix)
            // compute the blocks of H^k
            Vec3 b_i = -1* A_.transpose() * edge.information_ * e_;
            Vec3 b_j = -1* B_.transpose() * edge.information_ * e_;
            Mat3x3 H_ii = A_.transpose() *  edge.information_ * A_;
            Mat3x3 H_ij = A_.transpose() *  edge.information_ * B_;
            Mat3x3 H_jj = B_.transpose() *  edge.information_ * B_;
            // accumulate the blocks in H and b
            H_.block(i,i,i_dim,i_dim) += H_ii;
            H_.block(j,j,j_dim,j_dim) += H_jj;
            H_.block(i,j,i_dim,j_dim) += H_ij;
            H_.block(j,i,j_dim,i_dim) += H_ij.transpose();

            b_.block(i,0,i_dim,1) += b_i;
            b_.block(j,0,j_dim,1) += b_j;

        } else if (edge.edge_type_=="L") {
            // pose-landmark constraint
            j_dim = 2; // Jacobian block w.r.t. j 2x3 here

            // edge.fromIdx and edge.toIdx describe the location of
            // the first element of the pose and the landmark in the state vector
            // You should use also this index when updating the elements
            // of the H matrix and the vector b.

            Vec3 x1 = g_.x_.block(i, 0, 3, 1); // the robot pose
            Vec2 x2 = g_.x_.block(j, 0, 2, 1); // the landmark

            // Computing the error and the Jacobians, using measurement "z"
            linearize_pose_landmark_constraint(x1, x2, edge.measurement_);
            // Update H matrix and vector b, using Omega (the information matrix)
            // compute the blocks of H^k
            Vec3 b_i = -1* A_.transpose() * edge.information_ * e_;
            Vec2 b_j = -1* B_.transpose() * edge.information_ * e_;
            Mat3x3 H_ii = A_.transpose() *  edge.information_ * A_;
            Mat3x2 H_ij = A_.transpose() *  edge.information_ * B_;
            Mat2x2 H_jj = B_.transpose() *  edge.information_ * B_;
            // accumulate the blocks in H and b
            H_.block(i,i,i_dim,i_dim) += H_ii;
            H_.block(j,j,j_dim,j_dim) += H_jj;
            H_.block(i,j,i_dim,j_dim) += H_ij;
            H_.block(j,i,j_dim,i_dim) += H_ij.transpose();
            b_.block(i,0,i_dim,1) += b_i;
            b_.block(j,0,j_dim,1) += b_j;
        }
    }

    // Add the prior for one pose of the edges
    // fixes one node to remain at its current location
    // note that the system (H b) is obtained only from
    // relative constraints. H is not full rank.
    // we solve the problem by anchoring the position of the first vertex.
    // this can be expressed by adding the equation
    //   deltax(1:3,1)=0;

    // set upper-left 3x3 block of matrix
    H_.block(0,0,3,3) += Mat3x3::Identity(3,3);

    // save_matrix_image( H.copy(), iter, dataset_name )
}








/*
 * Compute the error of a pose-pose constraint. We
    compute homogeneous matrix transforms from the pose
    vectors (x,y,theta).

    We extract the poses of the vertices and the mean of this k'th edge.
    Computes the taylor expansion of the error function of this k'th edge.

    Args:
    -	v_i: 3x1 vector (x,y,theta) of the first robot pose
    -	v_j: 3x1 vector (x,y,theta) of the second robot pose
    -	z 3x1 vector (x,y,theta) of the measurement

            Returns:
    -	e: 3x1 error of the constraint -  e_k(x)
    -	A: 3x3 Jacobian wrt x1 - d e_k(x) / d(x_i)
    -	B: 3x3 Jacobian wrt x2 -  d e_k(x) / d(x_j)
 */
void SparseGaussNewtonOptimizer::linearize_pose_pose_constraint(Vec3 v_i, Vec3 v_j, Vec3 z_ij)
{
    e_.resize(3,1);
    A_.resize(3,3);
    B_.resize(3,3);

    // compute the homogeneous transforms of the previous solutions
    SE2 Z_ij(z_ij);
    SE2 w_T_i(v_i);
    SE2 w_T_j(v_j);

    // compute the displacement between x_i and x_j
    Mat3x3 f_ij = w_T_i.inverse() * w_T_j.mat_3x3_;

    // this below is too long to explain, to understand it derive it by hand
    Vec2 dt_ij = w_T_j.t_ - w_T_i.t_;

    A_ = -1 * Mat3x3::Identity(3,3);
    // left 2x2 block
    A_.block(0,0,2,2) = -w_T_i.R_.transpose();
    // far-right column
    A_.block(0,2,2,1) = w_T_i.dRT_dtheta_ * dt_ij;

    B_ = Mat3x3::Identity(3,3);
    // left 2x2 block
    B_.block(0,0,2,2) = w_T_i.R_.transpose();
    // far-right column
    B_.block(0,2,2,1) = Vec2::Zero();

    Mat3x3 ztinv = Z_ij.inverse();
    e_ = SE2_mat(ztinv * f_ij).as_pose_vector();
    // far-right column
    ztinv.block(0,2,2,1) = Vec2::Zero();
    A_ = ztinv * A_;
    B_ = ztinv * B_;

}

/*
 *  Compute the error of a pose-landmark constraint
 *  and the Jacobians of the error.
 *
 *  Args:
 *  -	x: 3x1 vector (x,y,theta) of the robot pose
 *  -	l: 2x1 vector (x,y) of the landmark
 *  -	z: 2x1 vector (x,y) of the measurement, the position of the landmark in
 *  the coordinate frame of the robot given by the vector x
 *
 *          Output
 *  -	e: 2x1 error of the constraint
 *  -	A: 2x3 Jacobian wrt x
 *  -	B: 2x2 Jacobian wrt l
 */
void SparseGaussNewtonOptimizer::linearize_pose_landmark_constraint(Vec3 x_i, Vec2 l, Vec2 z)
{
    e_.resize(2,1);
    A_.resize(2,3);
    B_.resize(2,2);

    SE2 w_T_i(x_i);
    e_ = w_T_i.R_.transpose() * (l - w_T_i.t_) - z;

    // left 2x2 block
    A_.block(0,0,2,2) = -w_T_i.R_.transpose();
    // far right column
    A_.block(0,2,2,1) = w_T_i.dRT_dtheta_ * (l - w_T_i.t_);
    B_ = w_T_i.R_.transpose();
}
