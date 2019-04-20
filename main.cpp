

#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <cmath>
#include <Eigen/Core>

#include "se2.h"


int main() {
    std::cout << "Hello, World!" << std::endl;
    Vec3 v;
    v(0,0) = 5.0;
    v(1,0) = -6.0;
    v(2,0) = M_PI / 4.0;
    SE2 w_T_i = SE2(v);
//    std::cout << w_T_i.v_ << std::endl;
    std::cout << w_T_i.mat_3x3_ << std::endl;
//    std::cout << w_T_i.R_ << std::endl;
//    std::cout << w_T_i.dRT_dtheta_ << std::endl;

    Mat3x3 inv = w_T_i.inverse();
    std::cout << w_T_i.inverse() << std::endl;
    return 0;
}