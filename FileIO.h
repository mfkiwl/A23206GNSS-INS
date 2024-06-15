#pragma once
#include <Eigen/Dense>
#include <vector>

const double deg_to_rad = 0.01745329252;

void fileToMatrix(const std::string& input_profile_name, int no_columns,
    Eigen::MatrixXd& in_profile, int& no_epochs);

bool Read_profile(const std::string& input_profile_name, const std::string& input_profile_name1, const std::string& input_profile_name2,
    Eigen::MatrixXd& in_profile, int& no_epochs, Eigen::MatrixXd& in_profile1, int& no_epochs1, Eigen::MatrixXd& in_profile2, int& no_epochs2);