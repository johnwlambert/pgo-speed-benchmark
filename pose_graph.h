//
// Created by John Lambert-Admin on 4/19/19.
//

#pragma once


#include <string>
#include <cmath>
#include <map>
#include <vector>

#include <Eigen/Dense>
#include "string_utils.h"


using VectorXd = Eigen::VectorXd; // equiv. to Matrix<double,Dynamic,1>
using MatrixXd = Eigen::MatrixXd; // equiv. to Matrix<double,Dynamic,Dynamic>


/*
 * Pose Graph Vertex
 * v_id is "vertex ID"
 */
class VertexPGO {
    public:
        size_t v_id_;
        size_t x_offset_idx_;
        size_t dim_;

        VertexPGO(size_t v_id, size_t x_offset_idx, size_t dim);
};

/*
 * Pose Graph Edge
 * "from_v_id" is "from" vertex ID
 * "to_v_id" is "to" vertex ID
 */
class EdgePGO{

    public:
        std::string edge_type_;
        size_t from_v_id_;
        size_t to_v_id_;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> measurement_;
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> information_;

    EdgePGO(std::string edge_type, size_t from_v_id, size_t to_v_id,
            VectorXd measurement, MatrixXd information);
};


class PoseGraph2D {
    public:
        PoseGraph2D(std::string dataset_name);
        void get_poses_landmarks();

        std::string dataset_name_;
        std::vector<EdgePGO> edges_;
        std::map<size_t, VertexPGO*> vertex_map_;

        // state vector (concatenated pose vectors)
        VectorXd x_;
        VectorXd pose_start_idxs_;
        VectorXd landmark_start_idxs_;

        private:
                void read_vertex_data(std::string vertices_fpath);
                void read_edge_data(std::string edges_fpath);
                void read_initial_state_vector(std::string initial_state_fpath);
};
