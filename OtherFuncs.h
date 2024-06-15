#pragma once
#include <Eigen/Dense>
#include <vector>

class LC_KF_config
{
public:
    double init_att_unc;
    double init_vel_unc;
    double init_pos_unc;
    double init_b_a_unc;
    double init_b_g_unc;
    double gyro_noise_PSD;
    double accel_noise_PSD;
    double accel_bias_PSD;
    double gyro_bias_PSD;
    double pos_meas_SD;
    double vel_meas_SD;
};

class GNSS_config
{
public:
    double epoch_interval;
    Eigen::Vector3d init_est_r_ea_e;
    double no_sat;
    double r_os;
    double inclination;
    double const_delta_lambda;
    double const_delta_t;
    double mask_angle;
    double SIS_err_SD;
    double zenith_iono_err_SD;
    double zenith_trop_err_SD;
    double code_track_err_SD;
    double rate_track_err_SD;
    double rx_clock_offset;
    double rx_clock_drift;
};

class IMU_errors
{
public:
    Eigen::Vector3d b_a;
    Eigen::Vector3d b_g;
    Eigen::Matrix3d M_a;
    Eigen::Matrix3d M_g;
    Eigen::Matrix3d G_g;
    double accel_noise_root_PSD;
    double gyro_noise_root_PSD;
    double accel_quant_level;
    double gyro_quant_level;
};

class InitializationErrors
{
public:
    Eigen::Vector3d delta_eul_nb_n;
};

Eigen::VectorXd generateRandomVector();

void IMU_model(double tor_i, Eigen::Vector3d true_f_ib_b, Eigen::Vector3d true_omega_ib_b, IMU_errors IMU_errors, Eigen::Matrix<double, 6, 1>& old_quant_residuals, Eigen::Vector3d& meas_f_ib_b, Eigen::Vector3d& meas_omega_ib_b, Eigen::Matrix<double, 6, 1>& quant_residuals);

Eigen::Matrix<double, 15, 15> Initialize_LC_P_matrix(const LC_KF_config LC_KF_config);

void LC_KF_Epoch(Eigen::Vector3d GNSS_r_eb_e, Eigen::Vector3d GNSS_v_eb_e, double tor_s,
    Eigen::Matrix3d& est_C_b_e_old, Eigen::Vector3d& est_v_eb_e_old,
    Eigen::Vector3d& est_r_eb_e_old, Eigen::VectorXd& est_IMU_bias_old,
    Eigen::Matrix<double, 15, 15>& P_matrix_old, Eigen::Vector3d meas_f_ib_b,
    double est_L_b_old, LC_KF_config LC_KF_config,
    Eigen::Matrix3d& est_C_b_e_new, Eigen::Vector3d& est_v_eb_e_new,
    Eigen::Vector3d& est_r_eb_e_new, Eigen::VectorXd& est_IMU_bias_new,
    Eigen::Matrix<double, 15, 15>& P_matrix_new);

void Loosely_coupled_INS_GNSS(const Eigen::MatrixXd& in_profile, int no_epochs, const Eigen::MatrixXd& in_profile1, int no_epochs1,
    const Eigen::MatrixXd& in_profile2, int no_epochs2, const InitializationErrors initialization_errors,
    const IMU_errors IMU_errors, const GNSS_config& gnss_config, const LC_KF_config LC_KF_config,
    Eigen::MatrixXd& out_profile1, Eigen::MatrixXd& out_errors, Eigen::MatrixXd& out_imu_bias_est,
    Eigen::MatrixXd& out_clock, Eigen::MatrixXd& out_kf_sd);
