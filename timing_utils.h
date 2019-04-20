//
// Created by John Lambert-Admin on 4/19/19.
//

#pragma once

#include <numeric>

/*
 * In order to obtain wall-clock performance benchmarks, we average the runtime for the algorithms specific phases
 * across each run.
 */
double inline compute_elapsed_wall_clock_time(std::chrono::steady_clock::time_point start_time,
                                              std::chrono::steady_clock::time_point current_time)
{
    std::chrono::duration<double> elapsed_secs = current_time - start_time;
    return static_cast<double>(elapsed_secs.count());
}
