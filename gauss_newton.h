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

/*
 * Args:
    -	g: graph
    -	numIterations: the number of iterations of Gauss-Newton
 */
void run_lsslam(PoseGraph2D g, std::string dataset_name)
{
    bool visualize=false;
    size_t max_iters = 100;
    //// plot the initial state of the graph
    //plot_graph(fig, g, 0)

    //print('Initial error %f\n', compute_global_error(g))
    // maximum allowed dx
    double EPSILON = std::pow(10.0,-4);

    // Error
    double err = 0;

    // carry out the iterations
    for (size_t iter=0; iter<max_iters; iter++){

        std::cout << "Performing iteration " << i << std::endl;
        std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();
        VectorXd dx = linearize_and_solve(g, iter, dataset_name);
        double duration = compute_elapsed_wall_clock_time(start_time, std::chrono::steady_clock::now() );
        std::cout << "Iter " << iter << " took " << duration << " sec.";

        // Apply the solution to the state vector g.x
        g.x_ += dx;

        PoseGraph2D g = normalize_angles(g);
        if (visualize==true){
            // plot the current state of the graph
            //plot_graph(fig, g, iter)
            // err = compute_global_error(g)
        }

        std::cout << "Current error = " << err;

        // termination criterion
        if (dx.norm() < std::pow(1.0,-10) ){
            break;
        }
    }


    std::cout << "Final error =" << err << std::endl;

}



