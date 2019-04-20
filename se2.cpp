//
// Created by John Lambert-Admin on 4/19/19.
//
#include "se2.h"

SE2::SE2(Vec3 v): v_(v)
{
    //assert(v.shape==(3,))
    t_(0,0) = v(0,0);
    t_(1,0) = v(1,0);
    theta_ = v(2);

    c_ = std::cos(theta_);
    s_ = std::sin(theta_);

    R_(0,0) = c_;
    R_(0,1) = -1.0*s_;
    R_(1,0) = s_;
    R_(1,1) = c_;

    dRT_dtheta_(0,0) =-1.0*s_;
    dRT_dtheta_(0,1) = c_;
    dRT_dtheta_(1,0) = -1.0*c_;
    dRT_dtheta_(1,1) = -1.0*s_;

    // assert(R.shape == (2,2))
    // assert(t.shape == (2,))

    mat_3x3_ = Mat3x3::Identity(3,3);
    // access (2x2) block at (0,0)
    mat_3x3_.block(0,0,2,2) = R_;
    // access far right block (2,1) starting at (0,2)
     mat_3x3_.block(0,2,2,1) = t_;
}


/*
 * """
    Inverts a homogenous transform.
    Args:
    -	None

            Returns:
    -	T_inv: 3x3 Numpy array
    """
 */
Mat3x3 SE2::inverse()
{
    Mat3x3 T_inv = Eigen::Matrix3d::Identity(3,3);

    // first 2x2 block
    T_inv.block(0,0,2,2) = R_.transpose() ;
    // set 3rd column
    T_inv.block(0,2,2,1) = -1.0 * R_.transpose() * t_;
    return T_inv;
}
