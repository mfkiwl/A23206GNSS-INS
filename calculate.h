#pragma once
#include <Eigen/Dense>
#include <vector>

void Calculate_errors_NED(
    double est_L_b, double est_lambda_b, double est_h_b,
    Eigen::Vector3d est_v_eb_n, Eigen::Matrix3d est_C_b_n,
    double true_L_b, double true_lambda_b, double true_h_b,
    Eigen::Vector3d true_v_eb_n, Eigen::Matrix3d true_C_b_n,
    Eigen::Vector3d& delta_r_eb_n, Eigen::Vector3d& delta_v_eb_n, Eigen::Vector3d& delta_eul_nb_n);

Eigen::Vector3d CTM_to_Euler(Eigen::Matrix3d C);

void ECEF_to_NED(Eigen::Vector3d r_eb_e, Eigen::Vector3d v_eb_e, Eigen::Matrix3d C_b_e, double& L_b, double& lambda_b, double& h_b, Eigen::Vector3d& v_eb_n, Eigen::Matrix3d& C_b_n);

Eigen::Matrix3d Euler_to_CTM(Eigen::Vector3d eul);

Eigen::Vector3d Gravity_ECEF(Eigen::Vector3d r_eb_e);

Eigen::Matrix3d Initialize_NED_attitude(Eigen::Matrix3d C_b_n, Eigen::Vector3d delta_eul_nb_n);

void Nav_equations_ECEF(const double tor_i,
    const Eigen::Vector3d& old_r_eb_e,
    const Eigen::Vector3d& old_v_eb_e,
    const Eigen::Matrix3d& old_C_b_e,
    const Eigen::Vector3d& f_ib_b,
    const Eigen::Vector3d& omega_ib_b,
    Eigen::Vector3d& r_eb_e,
    Eigen::Vector3d& v_eb_e,
    Eigen::Matrix3d& C_b_e);

void NED_to_ECEF(Eigen::Vector3d& r_eb_e, Eigen::Vector3d& v_eb_e, Eigen::Matrix3d& C_b_e, double L_b,
    double lambda_b, double h_b, Eigen::Vector3d v_eb_n, Eigen::Matrix3d C_b_n);

void Radii_of_curvature(double L, double& R_N, double& R_E);

Eigen::Matrix3d Skew_symmetric(Eigen::Vector3d a);