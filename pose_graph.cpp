//
// Created by John Lambert-Admin on 4/19/19.
//

#include "pose_graph.h"
#include <fstream>
#include <vector>
#include <string>
#include <iostream>


/*
 * Constructor loads the graph from disk.
 */
PoseGraph2D::PoseGraph2D(std::string dataset_name):
        dataset_name_(dataset_name)
{
    std::string edges_fpath = "datasets/" + dataset_name + "/" + dataset_name + "_edges.txt";
    std::string vertices_fpath = "datasets/" + dataset_name + "/" + dataset_name + "_vertices.txt";
    std::string initial_state_fpath = "datasets/" + dataset_name + "/" + dataset_name + "_initial_state.txt";

    read_edge_data(edges_fpath);
    read_vertex_data(vertices_fpath);
    read_initial_state_vector(initial_state_fpath);
}


/*
 * Args:
    -   edges_fpath

    Returns:
    -   edges: list of EdgePGO objects
 */
void PoseGraph2D::read_edge_data(std::string edges_fpath)
{
    std::ifstream ifp(edges_fpath);
    std::string line;
    if (ifp.is_open())
    {
        while ( std::getline (ifp,line) )
        {
            auto edge_info = split_comma_delimited_string(line);
            std::string edge_type = edge_info[0];
            std::size_t from_v_id = str2size_t(edge_info[1]);
            std::size_t to_v_id = str2size_t(edge_info[2]);
            size_t dim;
            if (edge_type == "L"){
                dim = 2;
            } else if (edge_type == "P"){
                dim = 3;
            }
            VectorXd measurement;
            measurement.resize(dim,1);
            for (size_t i=0; i<dim;i++)
            {
                measurement(i,1) = std::stod(edge_info[i+3]);
            }

            MatrixXd information;
            information.resize(dim,dim);
            for (size_t i=0; i<dim*dim;i++)
            {
                auto row = i/3;
                auto col = i%3;
                information(row,col) = std::stod(edge_info[3+dim+i]);
            }
            EdgePGO edge(edge_type, from_v_id, to_v_id, measurement, information);
            edges_.push_back(edge);
        }
        ifp.close();
    }
}


/*
 * Args:
    -vertices_fpath:

    Returns:
    -vertex_map: map from vertex ID to VertexPGO objects.
 */
void PoseGraph2D::read_vertex_data(std::string vertices_fpath) {
    std::ifstream ifp(vertices_fpath);
    std::string line;
    if (ifp.is_open())
    {
        while ( std::getline (ifp,line) )
        {
            auto vertex_data = split_comma_delimited_string(line);
            size_t v_id = str2size_t(vertex_data[0]);
            size_t x_offset_idx = str2size_t(vertex_data[1]);
            size_t dim = str2size_t(vertex_data[2]);
            vertex_map_.insert(std::make_pair(v_id, new VertexPGO(v_id, x_offset_idx, dim) ));
        }
        ifp.close();
    }
}

/*
 * Each float is on its own line. This forms initial state vector "x".
 */
void PoseGraph2D::read_initial_state_vector(std::string initial_state_fpath)
{
    std::ifstream ifp(initial_state_fpath);
    std::string line;
    std::vector<double> state_data;
    if (ifp.is_open())
    {
        while ( std::getline (ifp,line) )
        {
            state_data.push_back(std::stod(line));
        }
        ifp.close();
    }

    size_t  state_vec_len = 0;
    x_.resize(state_vec_len, Eigen::NoChange);
    for (size_t i=0;i<state_vec_len;i++) {
        x_(i, 0) = state_data[i];
    }
}



/*
 * Extract the offset of the poses and the landmarks.
    Args:
    -   g

    Returns:
    -   poses
    -   landmarks
 */
void PoseGraph2D::get_poses_landmarks() {

    std::vector<size_t> pose_idx_vec;
    std::vector<size_t> landmark_idx_vec;

    // We don't know how many poses/landmarks we'll have, so we push back one at a time
    for (auto& kv : vertex_map_) {
        auto v_id = kv.first;
        auto v_pgo = kv.second;

        if (v_pgo->dim_ == 3){
            pose_idx_vec.push_back(v_pgo->x_offset_idx_);
        } else if (v_pgo->dim_==2){
            landmark_idx_vec.push_back(v_pgo->x_offset_idx_);
        }
    }
    // Now place them into an array
    size_t num_landmarks = landmark_idx_vec.size();
    size_t num_poses = pose_idx_vec.size();
    pose_start_idxs_.resize(num_poses);
    landmark_start_idxs_.resize(num_landmarks);
    for (size_t i=0; i < num_poses; i++)
    {
        pose_start_idxs_(i,0) = pose_idx_vec[i];
    }

    for (size_t i=0; i < num_landmarks; i++)
    {
        landmark_start_idxs_(i,0) = landmark_idx_vec[i];
    }
}



/*
 * Constructor
 */
VertexPGO::VertexPGO(size_t v_id, size_t x_offset_idx, size_t dim):
        v_id_(v_id),
        x_offset_idx_(x_offset_idx),
        dim_(dim)
{};


/*
 * Pose Graph Edge
 * "from_v_id" is "from" vertex ID
 * "to_v_id" is "to" vertex ID
 */
EdgePGO::EdgePGO(std::string edge_type, size_t from_v_id, size_t to_v_id,
            VectorXd measurement, MatrixXd information):
            edge_type_(edge_type), from_v_id_(from_v_id), to_v_id_(to_v_id),
                measurement_(measurement), information_(information)

    {};

