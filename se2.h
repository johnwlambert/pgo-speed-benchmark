//
// Created by John Lambert-Admin on 4/19/19.
//

#pragma once

#include <Eigen/Core>

using Vec2 = Eigen::Matrix<double, 2, 1>;
using Vec3 = Eigen::Matrix<double, 3, 1>;
using Mat2x2 = Eigen::Matrix<double, 2, 2>;
using Mat2x3 = Eigen::Matrix<double, 2, 3>;
using Mat3x3 = Eigen::Matrix<double, 3, 3>;

static constexpr size_t DOF = 3;

/*
 * Computes the homogeneous transform matrix A of the pose vector v.
 * Args:
    -	v: (3,1) vector

   Creates
    -	mat3x3: 3x3 Numpy array
    -	R: Numpy array of shape (2,2)
    -	t: Numpy array of shape (2,1)
 */
class SE2 {
    public:
        SE2(Vec3 v);
        Mat3x3 inverse();

        Mat2x2 R_;
        Mat2x2 dRT_dtheta_;
        Mat3x3 mat_3x3_;
        Vec3 v_;
        Vec2 t_;

    private:
        double theta_;
        double c_;
        double s_;
};

class SE2_mat{
    public:
        SE2_mat(Mat3x3 mat_3x3)
        {
            // access (2x2) block at (0,0)
            Mat2x2 R_ = mat_3x3.block(0,0,2,2);

            // access far right block (2,1) starting at (0,2)
            Vec2 t_ = mat_3x3.block(0,2,2,1);
        }
        Mat2x2 R_;
        Vec2 t_;

    /*
     *  Computes the pose vector v from a homogeneous transform A.
        Args:
        -	A: 3x3 Numpy array
        Returns:
        -	v
     */
    Vec3 as_pose_vector()
    {
        Vec3 v;
        v(0,0) = t_(0,0);
        v(1,0) = t_(1,0); // could use comma initializer instead
        double theta = std::atan2(R_(1,0),R_(0,0));
        v(2) = theta;
        return v;
    }
};