//
// Created by John Lambert-Admin on 4/19/19.
//

#pragma once

#include <cmath>
#include "pose_graph.h


/*
 * Move to [-pi,pi] range
 */
double clip_angle_to_arctan(double theta)
{
    double s = std::sin(theta);
    double c = std::cos(theta);
    return std::atan2(s,c);
}


/*
 * Normalize the angles between -PI and PI
 * Only use robot poses that do not correspond to landmark poses
 *      (those have a 3rd element theta)

    Args:
    -	g: graph, with poses OFF the manifold SO(2)

    Returns:
    -	g:graph, with poses ON the manifold SO(2)
 */
PoseGraph2D normalize_angles(PoseGraph2D pg)
{
    for (auto & kv: pg.vertex_map_) {
        size_t v_id = kv.first;
        VertexPGO v_pgo = kv.second;

        if (v_pgo.dim_ == 3) {
            // get the offset of the x element of this pose vector
            // in the whole state vector
            size_t offs = v_pgo.x_offset_idx_;

            // theta from the first robot pose -- each is (x,y,theta)
            size_t theta_idx = offs + 2;
            pg.x_(theta_idx,0) = clip_angle_to_arctan(pg.x_(theta_idx) );
        }
    }
    return pg;
}

