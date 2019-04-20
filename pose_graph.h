//
// Created by John Lambert-Admin on 4/19/19.
//

#pragma once

#include <string>
#include <cmath>
#include <map>
#include <vector>

#include <Eigen/Dense>

/*
 * Pose Graph Vertex
 * v_id is "vertex ID"
 */
class VertexPGO {

public:
    VertexPGO(size_t v_id, size_t x_offset_idx,size_t dim);

    size_t v_id_;
    size_t x_offset_idx_;
    size_t dim_;
};

/*
 * Pose Graph Edge
 * "from_v_id" is "from" vertex ID
 * "to_v_id" is "to" vertex ID
 */
class EdgePGO{

    EdgePGO(std::string edge_type, size_t from_v_id, size_t to_v_id,
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> measurement,
            Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> information);

    public:
        std::string edge_type_;
        size_t from_v_id_;
        size_t to_v_id_;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> measurement_;
        Eigen::Matrix<int, Eigen::Dynamic, Eigen::Dynamic> information_;
};

