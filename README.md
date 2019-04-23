
## 2D Pose SLAM Benchmark
Comparison of sparse system solving in C (CSparse) vs. C++ (GTSAM)

### For CSparse
cmake  -DCMAKE_BUILD_TYPE=Release ..
make -j8

### For GTSAM
cmake -DGTSAM_TOOLBOX_INSTALL_PATH:PATH=$HOME/toolbox ..