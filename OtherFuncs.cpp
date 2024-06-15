#include <iostream>
#include <Eigen/Dense>
#include <iomanip> // 包含 setprecision
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include "calculate.h"
#include "OtherFuncs.h"

using namespace std;
using namespace Eigen;

VectorXd generateRandomVector()
{
	unsigned seed = chrono::system_clock::now().time_since_epoch().count();
	default_random_engine gen(seed);

	normal_distribution<> normal_dist(0, 1);

	VectorXd random_vector(6);
	for (int n = 0; n < 6; ++n)
	{
		random_vector(n) = normal_dist(gen);
	}

	return random_vector;
}

void IMU_model(double tor_i, Eigen::Vector3d true_f_ib_b, Eigen::Vector3d true_omega_ib_b, IMU_errors IMU_errors, Eigen::Matrix<double, 6, 1>& old_quant_residuals, Eigen::Vector3d& meas_f_ib_b, Eigen::Vector3d& meas_omega_ib_b, Eigen::Matrix<double, 6, 1>& quant_residuals)
{
	Vector3d accel_noise(0.0, 0.0, 0.0);
	Vector3d gyro_noise(0.0, 0.0, 0.0);

	VectorXd random_vector = generateRandomVector();

	Vector3d accel_distribution_samples(0.0, 0.0, 0.0);
	accel_distribution_samples = random_vector.segment(0, 3);
	//accel_distribution_samples << 0.154391271801951, 0.325876736354160, -0.355654503308562;

	Vector3d gyro_distribution_samples(0.0, 0.0, 0.0);
	gyro_distribution_samples = random_vector.segment(3, 3);
	//gyro_distribution_samples << -0.586295008994886, 1.224761441269211, -0.444573554832088;

	// cout << accel_distribution_samples << endl;
	// cout << gyro_distribution_samples << endl;

	// 输入的时间间隔 tor_i 大于0时，生成加速度计和陀螺仪的噪声项。
	if (tor_i > 0)
	{
		// 加速度计噪声存储在 accel_noise 中，陀螺仪噪声存储在 gyro_noise 中。
		accel_noise = accel_distribution_samples * IMU_errors.accel_noise_root_PSD / sqrt(tor_i);
		gyro_noise = gyro_distribution_samples * IMU_errors.gyro_noise_root_PSD / sqrt(tor_i);

		// cout << accel_noise << endl;
		// cout << gyro_noise << endl;
	}
	else
	{
		// 如果时间间隔 tor_i 小于等于0，则将噪声项设置为零向量。
		accel_noise = Vector3d::Zero();
		gyro_noise = Vector3d::Zero();
	}

	Vector3d uq_f_ib_b = IMU_errors.b_a + (Matrix3d::Identity() + IMU_errors.M_a) * true_f_ib_b + accel_noise;
	// cout << "uq_f_ib_b" << endl;
	// cout << uq_f_ib_b << endl;
	 //cout << "true_f_ib_b" << endl;
	 //cout << true_f_ib_b << endl;
	 //cout << "accel_noise" << endl;
	 //cout << accel_noise << endl;


	Vector3d uq_omega_ib_b = IMU_errors.b_g + (Matrix3d::Identity() + IMU_errors.M_g) * true_omega_ib_b +
		IMU_errors.G_g * true_f_ib_b + gyro_noise;
	// cout << "uq_omega_ib_b" << endl;
	// cout << uq_omega_ib_b << endl;

	// 对加速度计和陀螺仪的输出进行量化处理。
	if (IMU_errors.accel_quant_level > 0)
	{
		Vector3d temp_residuals(0.0, 0.0, 0.0);
		Vector3d rounded_residuals(0.0, 0.0, 0.0);
		// double round(double x);四舍五入浮点数到最接近的整数值

		/*cout << uq_f_ib_b << endl;*/
		temp_residuals = (uq_f_ib_b + old_quant_residuals.segment(0, 3)) / IMU_errors.accel_quant_level;
		/*cout << temp_residuals << endl;*/
		for (int i = 0; i < 3; i++)
		{
			rounded_residuals(i) = round(temp_residuals(i));
		}

		//cout<< rounded_residuals <<endl;

		meas_f_ib_b = IMU_errors.accel_quant_level * rounded_residuals;

		//cout << meas_f_ib_b << endl;

		// segment<Size>(Start) 其中，Size 是切片的大小，Start 是切片的起始位置。
		// quant_residuals.segment<3>(0)等价于quant_residuals.segment(0, 3)
		quant_residuals.segment(0, 3) = uq_f_ib_b + old_quant_residuals.segment(0, 3) - meas_f_ib_b;
	}
	else
	{
		meas_f_ib_b = uq_f_ib_b;
		// Vector3d实际上是 Matrix<double, 3, 1> 的别名。因此，将 Vector3d 的值赋给 Matrix<double, 6, 1> 对象时，实际上是将 Vector3d 的值复制到Matrix<double, 6, 1> 的前三个元素中，而其他元素保持不变。
		quant_residuals.segment(0, 3) = Vector3d::Zero();
	}

	if (IMU_errors.gyro_quant_level > 0)
	{
		Vector3d temp_residuals(0.0, 0.0, 0.0);
		Vector3d rounded_residuals(0.0, 0.0, 0.0);

		temp_residuals = (uq_omega_ib_b + old_quant_residuals.segment(3, 3)) / IMU_errors.gyro_quant_level;
		//cout << temp_residuals << endl;

		for (int i = 0; i < 3; i++)
		{
			rounded_residuals(i) = round(temp_residuals(i));
		}

		//cout << rounded_residuals << endl;

		meas_omega_ib_b = IMU_errors.gyro_quant_level * rounded_residuals;

		//cout << meas_omega_ib_b << endl;

		quant_residuals.segment(3, 3) = uq_omega_ib_b + old_quant_residuals.segment(3, 3) - meas_omega_ib_b;
	}
	else
	{
		meas_omega_ib_b = uq_omega_ib_b;
		quant_residuals.segment(3, 3) = Vector3d::Zero();
	}
}

Matrix<double, 15, 15> Initialize_LC_P_matrix(const LC_KF_config LC_KF_config)
{
	Matrix<double, 15, 15> P_matrix = Matrix<double, 15, 15>::Zero();

	// P.block<rows, cols>(i, j)
	// P.block(i, j, rows, cols)
	// 从i行j列开始的rows行cols列块
	//  Initialize error covariance matrix
	P_matrix.block<3, 3>(0, 0) = Matrix3d::Identity() * LC_KF_config.init_att_unc * LC_KF_config.init_att_unc;
	P_matrix.block<3, 3>(3, 3) = Matrix3d::Identity() * LC_KF_config.init_vel_unc * LC_KF_config.init_vel_unc;
	P_matrix.block<3, 3>(6, 6) = Matrix3d::Identity() * LC_KF_config.init_pos_unc * LC_KF_config.init_pos_unc;
	P_matrix.block<3, 3>(9, 9) = Matrix3d::Identity() * LC_KF_config.init_b_a_unc * LC_KF_config.init_b_a_unc;
	P_matrix.block<3, 3>(12, 12) = Matrix3d::Identity() * LC_KF_config.init_b_g_unc * LC_KF_config.init_b_g_unc;

	return P_matrix;
}

void LC_KF_Epoch(Vector3d GNSS_r_eb_e, Vector3d GNSS_v_eb_e, double tor_s,
	Matrix3d& est_C_b_e_old, Vector3d& est_v_eb_e_old,
	Vector3d& est_r_eb_e_old, VectorXd& est_IMU_bias_old,
	Matrix<double, 15, 15>& P_matrix_old, Vector3d meas_f_ib_b,
	double est_L_b_old, LC_KF_config LC_KF_config,
	Matrix3d& est_C_b_e_new, Vector3d& est_v_eb_e_new,
	Vector3d& est_r_eb_e_new, VectorXd& est_IMU_bias_new,
	Matrix<double, 15, 15>& P_matrix_new)
{

	const double c = 299792458;
	const double omega_ie = 7.292115E-5;
	const double R_0 = 6378137;
	const double e = 0.0818191908425;

	Vector3d a(0, 0, omega_ie);
	Matrix3d Omega_ie = Skew_symmetric(a);

	// 1. Determine transition matrix using (14.50) (first-order approx)
	Matrix<double, 15, 15> Phi_matrix = Matrix<double, 15, 15>::Identity();
	// matrix.block(i,j,p,q) 提取块大小为(p,q),起始于(i,j)
	Phi_matrix.block(0, 0, 3, 3) -= Omega_ie * tor_s;
	Phi_matrix.block(0, 12, 3, 3) = est_C_b_e_old * tor_s;
	Phi_matrix.block(3, 0, 3, 3) = -tor_s * Skew_symmetric(est_C_b_e_old * meas_f_ib_b);
	Phi_matrix.block(3, 3, 3, 3) -= 2 * Omega_ie * tor_s;

	double geocentric_radius = R_0 / sqrt(1.0 - e * e * sin(est_L_b_old) * sin(est_L_b_old)) * sqrt(cos(est_L_b_old) * cos(est_L_b_old) + (1.0 - e * e) * (1.0 - e * e) * sin(est_L_b_old) * sin(est_L_b_old));

	// transpose()转置
	Phi_matrix.block(3, 6, 3, 3) = -tor_s * 2 * Gravity_ECEF(est_r_eb_e_old) /
		geocentric_radius * est_r_eb_e_old.transpose() / sqrt(est_r_eb_e_old.transpose() * est_r_eb_e_old);

	Phi_matrix.block(3, 9, 3, 3) = est_C_b_e_old * tor_s;
	Phi_matrix.block(6, 3, 3, 3) = Matrix3d::Identity() * tor_s;
	// 相同
	//cout << "Phi_matrix" << endl;
	//cout << setprecision(15) << Phi_matrix << endl;

	// 2. Determine approximate system noise covariance matrix using (14.82)
	Matrix<double, 15, 15> Q_prime_matrix = Matrix<double, 15, 15>::Zero();
	Q_prime_matrix.block(0, 0, 3, 3) = Matrix3d::Identity() * LC_KF_config.gyro_noise_PSD * tor_s;
	Q_prime_matrix.block(3, 3, 3, 3) = Matrix3d::Identity() * LC_KF_config.accel_noise_PSD * tor_s;
	Q_prime_matrix.block(9, 9, 3, 3) = Matrix3d::Identity() * LC_KF_config.accel_bias_PSD * tor_s;
	Q_prime_matrix.block(12, 12, 3, 3) = Matrix3d::Identity() * LC_KF_config.gyro_bias_PSD * tor_s;
	// 相同
	//cout << "Q_prime_matrix" << endl;
	//cout << setprecision(15) << Q_prime_matrix << endl;
	// 3. Propagate state estimates using (3.14) noting that all states are zero
	Matrix<double, 15, 1> x_est_propagated = Matrix<double, 15, 1>::Zero();

	// 4. Propagate state estimation error covariance matrix using (3.46)
	Matrix<double, 15, 15> P_matrix_propagated = Matrix<double, 15, 15>::Zero();
	P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
	// 相同
	//cout << "P_matrix_propagated" << endl;
	//cout << setprecision(15) << P_matrix_propagated.block(0, 0, 15, 3) << endl;
	//cout << setprecision(15) << P_matrix_propagated.block(0, 3, 15, 3) << endl;
	//cout << setprecision(15) << P_matrix_propagated.block(0, 6, 15, 3) << endl;
	//cout << setprecision(15) << P_matrix_propagated.block(0, 9, 15, 3) << endl;
	//cout << setprecision(15) << P_matrix_propagated.block(0, 12, 15, 3) << endl;
	// MEASUREMENT UPDATE PHASE

	// 5. Set-up measurement matrix using (14.115)
	Matrix<double, 6, 15> H_matrix = MatrixXd::Zero(6, 15);
	H_matrix.block(0, 6, 3, 3) = -Matrix3d::Identity();
	H_matrix.block(3, 3, 3, 3) = -Matrix3d::Identity();
	//cout << "H_matrix" << endl;
	//cout << setprecision(15) << H_matrix << endl;
	// 6. Set-up measurement noise covariance matrix assuming all components of GNSS position and velocity are independent and have equal variance.
	Matrix<double, 6, 6> R_matrix = Matrix<double, 6, 6>::Zero();
	R_matrix.block(0, 0, 3, 3) = Matrix3d::Identity() * pow(LC_KF_config.pos_meas_SD, 2);
	R_matrix.block(0, 3, 3, 3) = Matrix3d::Zero();
	R_matrix.block(3, 0, 3, 3) = Matrix3d::Zero();
	R_matrix.block(3, 3, 3, 3) = Matrix3d::Identity() * pow(LC_KF_config.vel_meas_SD, 2);
	// 相同
	//cout << "R_matrix" << endl;
	//cout << setprecision(15) << R_matrix << endl;
	// 7. Calculate Kalman gain using (3.21)
	// inverse()计算矩阵的逆矩阵
	Matrix<double, 15, 6> K_matrix = Matrix<double, 15, 6>::Zero();
	K_matrix = P_matrix_propagated * H_matrix.transpose() * (H_matrix * P_matrix_propagated * H_matrix.transpose() + R_matrix).inverse();
	//cout << "K_matrix" << endl;
	//cout << setprecision(15) << K_matrix.block(0, 0, 15, 3) << endl;
	//cout << setprecision(15) << K_matrix.block(0, 3, 15, 3) << endl;

	// 8. Formulate measurement innovations using (14.102), noting that zero
	VectorXd delta_z(6);
	delta_z << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	delta_z.segment(0, 3) = GNSS_r_eb_e - est_r_eb_e_old;
	delta_z.segment(3, 3) = GNSS_v_eb_e - est_v_eb_e_old;
	//cout << "delta_z" << endl;
	//cout << setprecision(15) << delta_z << endl;
	// 9. Update state estimates using (3.24)
	VectorXd x_est_new(15);
	x_est_new = x_est_propagated + K_matrix * delta_z;
	//cout << "x_est_new" << endl;
	//cout << setprecision(15) << x_est_new << endl;

	// 10. Update state estimation error covariance matrix using (3.25)
	P_matrix_new = (Matrix<double, 15, 15>::Identity() - K_matrix * H_matrix) * P_matrix_propagated;
	//cout << "P_matrix_new" << endl;
	//cout << setprecision(15) << P_matrix_new << endl;
	// CLOSED-LOOP CORRECTION

	// Correct attitude, velocity, and position using (14.7-9)
	est_C_b_e_new = (Matrix3d::Identity() - Skew_symmetric(x_est_new.segment(0, 3))) * est_C_b_e_old;
	est_v_eb_e_new = est_v_eb_e_old - x_est_new.segment(3, 3);
	est_r_eb_e_new = est_r_eb_e_old - x_est_new.segment(6, 3);

	// Update IMU bias estimates
	est_IMU_bias_new = est_IMU_bias_old + x_est_new.segment(9, 6);
}

void Loosely_coupled_INS_GNSS(const MatrixXd& in_profile, int no_epochs, const MatrixXd& in_profile1, int no_epochs1,
	const MatrixXd& in_profile2, int no_epochs2, const InitializationErrors initialization_errors,
	const IMU_errors IMU_errors, const GNSS_config& gnss_config, const LC_KF_config LC_KF_config,
	MatrixXd& out_profile1, MatrixXd& out_errors, MatrixXd& out_imu_bias_est,
	MatrixXd& out_clock, MatrixXd& out_kf_sd)
{

	// Initialize true navigation solution
	double old_time = in_profile(0, 1);
	double true_L_b = in_profile(0, 2);
	double true_lambda_b = in_profile(0, 3);
	double true_h_b = in_profile(0, 4);
	// 在Eigen中，向量的默认情况下是列向量
	Vector3d true_v_eb_n = in_profile.row(0).segment(5, 3).transpose();
	Vector3d true_eul_nb = in_profile.row(0).segment(8, 3).transpose();
	Matrix3d true_C_b_n = Euler_to_CTM(true_eul_nb).transpose();

	Vector3d old_true_r_eb_e(0.0, 0.0, 0.0);
	Vector3d old_true_v_eb_e(0.0, 0.0, 0.0);
	Matrix3d old_true_C_b_e = Matrix3d::Zero();

	NED_to_ECEF(old_true_r_eb_e, old_true_v_eb_e, old_true_C_b_e, true_L_b, true_lambda_b, true_h_b, true_v_eb_n, true_C_b_n);

	Vector2d est_clock;
	est_clock = Vector2d::Zero();

	double old_est_L_b = in_profile2(0, 1);
	double old_est_lambda_b = in_profile2(0, 2);
	double old_est_h_b = in_profile2(0, 3);
	Vector3d old_est_v_eb_n = in_profile2.row(0).segment(4, 3).transpose();
	//cout << "old_est_v_eb_n" << endl;
	//cout << old_est_v_eb_n << endl;
	//cout << in_profile2.row(0) << endl;

	double est_L_b = old_est_L_b;

	Vector3d old_est_r_eb_e(0.0, 0.0, 0.0);
	Vector3d old_est_v_eb_e(0.0, 0.0, 0.0);
	Matrix3d C_b_e = Matrix3d::Zero();

	NED_to_ECEF(old_est_r_eb_e, old_est_v_eb_e, C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, Matrix3d::Zero());

	Matrix3d old_est_C_b_n = Initialize_NED_attitude(true_C_b_n, initialization_errors.delta_eul_nb_n);

	Vector3d temp1(0.0, 0.0, 0.0);
	Vector3d temp2(0.0, 0.0, 0.0);
	Matrix3d old_est_C_b_e = Matrix3d::Zero();
	NED_to_ECEF(temp1, temp2, old_est_C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, old_est_C_b_n);

	// Initialize output profile record and errors record
	out_profile1 = MatrixXd::Zero(no_epochs1, 10);
	out_errors = MatrixXd::Zero(no_epochs1, 10);

	// Generate output profile record for the first epoch
	out_profile1(0, 0) = old_time;
	out_profile1(0, 1) = old_est_L_b;
	out_profile1(0, 2) = old_est_lambda_b;
	out_profile1(0, 3) = old_est_h_b;
	out_profile1.row(0).segment(4, 3) = old_est_v_eb_n.transpose();
	out_profile1.row(0).segment(7, 3) = CTM_to_Euler(old_est_C_b_n.transpose()).transpose();

	// Determine errors for the first epoch
	Vector3d delta_r_eb_n(0.0, 0.0, 0.0);
	Vector3d delta_v_eb_n(0.0, 0.0, 0.0);
	Vector3d delta_eul_nb_n(0.0, 0.0, 0.0);

	Calculate_errors_NED(
		old_est_L_b, old_est_lambda_b, old_est_h_b,
		old_est_v_eb_n, old_est_C_b_n,
		true_L_b, true_lambda_b, true_h_b,
		true_v_eb_n, true_C_b_n,
		delta_r_eb_n, delta_v_eb_n, delta_eul_nb_n);

	out_errors(0, 0) = old_time;
	out_errors.row(0).segment(1, 3) = delta_r_eb_n.transpose();
	out_errors.row(0).segment(4, 3) = delta_v_eb_n.transpose();
	out_errors.row(0).segment(7, 3) = delta_eul_nb_n.transpose();

	// Initialize Kalman filter P matrix and IMU bias states
	Matrix<double, 15, 15> P_matrix = Initialize_LC_P_matrix(LC_KF_config);
	//cout << "LC_KF_config" << endl;
	//cout << LC_KF_config.accel_bias_PSD << endl;
	//cout << "P_matrix" << endl;
	//cout << P_matrix << endl;
	VectorXd est_IMU_bias = VectorXd::Zero(6);

	// Initialize IMU quantization residuals
	Matrix<double, 6, 1> quant_residuals;
	quant_residuals << 0, 0, 0, 0, 0, 0;

	out_imu_bias_est = Matrix<double, 1, 7>::Zero();
	out_imu_bias_est(0, 0) = old_time;
	out_imu_bias_est.row(0).segment(1, 6) = est_IMU_bias.transpose();
	out_clock = Matrix<double, 1, 3>::Zero();
	out_clock(0, 0) = old_time;
	out_clock.row(0).segment(1, 2) = est_clock;

	// Generate KF uncertainty record
	out_kf_sd = MatrixXd::Zero(1, 16);
	out_kf_sd(0, 0) = old_time;

	for (int i = 1; i <= 15; ++i)
	{
		out_kf_sd(0, i) = sqrt(P_matrix(i - 1, i - 1));
	}

	//cout << "out_kf_sd" << endl;
	//cout << out_kf_sd << endl;

	// Initialize GNSS model timing
	double time_last_GNSS = old_time;
	int GNSS_epoch = 1;

	// Progress bar
	string dots = "....................";
	string bars = "||||||||||||||||||||";
	string rewind = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
	cout << "Processing: " << dots;
	int progress_mark = 0;
	int progress_epoch = 0;

	// Main loop
	for (int epoch = 2; epoch <= no_epochs1; ++epoch)
	{
		// Update progress bar
		if ((epoch - progress_epoch) > (no_epochs1 / 20.0))
		{
			progress_mark++;
			progress_epoch = epoch;
			cout << rewind << bars.substr(0, progress_mark) << dots.substr(0, 20 - progress_mark);
		}

		// Input data from motion profile
		true_L_b = in_profile(epoch - 1, 1);
		true_lambda_b = in_profile(epoch - 1, 2);
		true_h_b = in_profile(epoch - 1, 3);
		true_v_eb_n = in_profile.row(epoch - 1).segment(4, 3).transpose();
		true_eul_nb = in_profile.row(epoch - 1).segment(7, 3).transpose();
		true_C_b_n = Euler_to_CTM(true_eul_nb).transpose();
		Vector3d true_r_eb_e = Vector3d::Zero();
		Vector3d true_v_eb_e = Vector3d::Zero();
		Matrix3d true_C_b_e = Matrix3d::Identity();
		NED_to_ECEF(true_r_eb_e, true_v_eb_e, true_C_b_e, true_L_b, true_lambda_b, true_h_b, true_v_eb_n, true_C_b_n);


		//cout << "true_r_eb_e" << endl;
		//cout << true_r_eb_e << endl;
		//cout << "true_v_eb_e" << endl;
		//cout << true_v_eb_e << endl;
		//cout << "true_C_b_e" << endl;
		//cout << true_C_b_e << endl;
		//cout << "true_L_b" << endl;
		//cout << true_L_b << endl;
		//cout << "true_lambda_b" << endl;
		//cout << true_lambda_b << endl;
		//cout << "true_h_b" << endl;
		//cout << true_h_b << endl;
		//cout << "true_v_eb_n" << endl;
		//cout << true_v_eb_n << endl;
		//cout << "true_C_b_n" << endl;
		//cout << true_C_b_n << endl;



		double time = in_profile1(epoch - 1, 0);
		Vector3d alpha_ib_b = in_profile1.row(epoch - 1).segment(1, 3).transpose();

		Vector3d v_ib_b = in_profile1.row(epoch - 1).segment(4, 3).transpose();

		double tor_i = time - old_time;
		Vector3d meas_f_ib_b_old = v_ib_b / tor_i;

		Vector3d meas_omega_ib_b_old = alpha_ib_b / tor_i;

		Vector3d meas_omega_ib_b(0.0, 0.0, 0.0);
		Vector3d meas_f_ib_b(0.0, 0.0, 0.0);

		Matrix<double, 6, 1> quant_residuals_new;
		quant_residuals_new << 0, 0, 0, 0, 0, 0;
		// Corrected the function call to use quant_residuals

		//meas_f_ib_b, meas_omega_ib_b, quant_residuals_new

		//cout << "IMU_errors.b_a" << endl;
		//cout << IMU_errors.b_a << endl;
		//cout << IMU_errors.b_g << endl;
		//cout << IMU_errors.M_a << endl;
		//cout << IMU_errors.M_g << endl;
		//cout << IMU_errors.G_g << endl;

		IMU_model(tor_i, meas_f_ib_b_old, meas_omega_ib_b_old, IMU_errors, quant_residuals, meas_f_ib_b, meas_omega_ib_b, quant_residuals_new);
		//cout << meas_f_ib_b << endl;
		//cout << meas_omega_ib_b << endl;
		//cout << quant_residuals_new << endl;
		quant_residuals = quant_residuals_new;
		meas_f_ib_b -= est_IMU_bias.segment(0, 3);
		meas_omega_ib_b -= est_IMU_bias.segment(3, 3);

		Vector3d est_r_eb_e(0.0, 0.0, 0.0);
		Vector3d est_v_eb_e(0.0, 0.0, 0.0);
		Matrix3d est_C_b_e = Matrix3d::Identity();

		double est_lambda_b = 0;
		double est_h_b = 0;

		Nav_equations_ECEF(tor_i, old_est_r_eb_e, old_est_v_eb_e, old_est_C_b_e, meas_f_ib_b, meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
		/*cout << est_r_eb_e << endl;
		cout << est_v_eb_e << endl;
		cout << est_C_b_e << endl;*/
		Vector3d est_v_eb_n(0.0, 0.0, 0.0);
		Matrix3d est_C_b_n = Matrix3d::Identity();
		// Determine whether to update GNSS simulation and run Kalman filter
		if ((time - time_last_GNSS) >= gnss_config.epoch_interval)
		{
			GNSS_epoch++;
			double tor_s = time - time_last_GNSS;
			time_last_GNSS = in_profile2(GNSS_epoch - 1, 0);

			if ((time_last_GNSS - time) > 0.5)
			{
				GNSS_epoch = GNSS_epoch - 1;
				time_last_GNSS = time_last_GNSS - 1;
				ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, est_L_b, est_lambda_b, est_h_b, est_v_eb_n, est_C_b_n);
				//cout << est_L_b << endl;
				//cout << est_lambda_b << endl;
				//cout << est_h_b << endl;
				//cout << est_v_eb_n << endl;
				//cout << est_C_b_n << endl;


				out_profile1(epoch - 1, 0) = time;
				out_profile1(epoch - 1, 1) = est_L_b;
				out_profile1(epoch - 1, 2) = est_lambda_b;
				out_profile1(epoch - 1, 3) = est_h_b;
				out_profile1.row(epoch - 1).segment(4, 3) = est_v_eb_n.transpose();
				out_profile1.row(epoch - 1).segment(7, 3) = CTM_to_Euler(est_C_b_n.transpose()).transpose();

				old_time = time;
				old_est_r_eb_e = est_r_eb_e;
				old_est_v_eb_e = est_v_eb_e;
				old_est_C_b_e = est_C_b_e;
			}

			old_est_L_b = in_profile2(GNSS_epoch - 1, 1);
			old_est_lambda_b = in_profile2(GNSS_epoch - 1, 2);
			old_est_h_b = in_profile2(GNSS_epoch - 1, 3);
			// 162
			double a, b, c;
			a = b = c = 0;
			Matrix3d A = Matrix3d::Zero();

			/*cout << "est_r_eb_e" << endl;
			cout << est_r_eb_e << endl;
			cout << "est_v_eb_e" << endl;
			cout << est_v_eb_e << endl;
			cout << "est_C_b_e" << endl;
			cout << est_C_b_e << endl;*/

			ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, a, b, c, est_v_eb_n, A);

			//cout << "est_v_eb_n" << endl;
			//cout << est_v_eb_n << endl;

			old_est_v_eb_n = est_v_eb_n;

			est_L_b = old_est_L_b;

			Vector3d GNSS_r_eb_e(0.0, 0.0, 0.0);
			Vector3d GNSS_v_eb_e(0.0, 0.0, 0.0);
			Matrix3d GNSS_C_b_e = Matrix3d::Zero();

			//cout << "old_est_L_b" << endl;
			//cout << old_est_L_b << endl;
			//cout << "old_est_lambda_b" << endl;
			//cout << old_est_lambda_b << endl;
			//cout << "old_est_h_b" << endl;
			//cout << old_est_h_b << endl;
			//cout << "old_est_v_eb_n" << endl;
			//cout << old_est_v_eb_n << endl;

			NED_to_ECEF(GNSS_r_eb_e, GNSS_v_eb_e, GNSS_C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, Matrix3d::Zero());

			Matrix<double, 15, 15> P_matrix_new = Matrix<double, 15, 15>::Zero();
			// Run Integration Kalman filter
			//cout << "GNSS_r_eb_e" << endl;
			//cout << GNSS_r_eb_e << endl;
			//cout << "GNSS_v_eb_e" << endl;
			//cout << GNSS_v_eb_e << endl;
			//cout << "tor_s" << endl;
			//cout << tor_s << endl;
			//cout << "est_C_b_e" << endl;
			//cout << est_C_b_e << endl;
			//cout << "est_v_eb_e" << endl;
			//cout << est_v_eb_e << endl;
			//cout << "est_r_eb_e" << endl;
			//cout << est_r_eb_e << endl;
			//cout << "est_IMU_bias" << endl;
			//cout << est_IMU_bias << endl;
			//cout << "P_matrix" << endl;
			//cout << P_matrix << endl;
			//cout << "meas_f_ib_b" << endl;
			//cout << meas_f_ib_b << endl;
			//cout << "est_L_b" << endl;
			//cout << est_L_b << endl;

			LC_KF_Epoch(GNSS_r_eb_e, GNSS_v_eb_e, tor_s, est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix,
				meas_f_ib_b, est_L_b, LC_KF_config, est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix_new);
			P_matrix = P_matrix_new;

			//cout << "est_C_b_e" << endl;
			//cout << est_C_b_e << endl;
			//cout << "est_v_eb_e" << endl;
			//cout << est_v_eb_e << endl;
			//cout << "est_r_eb_e" << endl;
			//cout << est_r_eb_e << endl;
			//cout << "est_IMU_bias" << endl;
			//cout << est_IMU_bias << endl;
			//cout << "P_matrix" << endl;
			//cout << P_matrix << endl;

			out_imu_bias_est.resize(GNSS_epoch, 7);
			out_imu_bias_est(GNSS_epoch - 1, 0) = time;
			out_imu_bias_est.row(GNSS_epoch - 1).segment(1, 6) = est_IMU_bias.transpose();

			out_clock.resize(GNSS_epoch, 3);
			out_clock(GNSS_epoch - 1, 0) = time;
			out_clock.row(GNSS_epoch - 1).segment(1, 2) = est_clock;

			out_kf_sd.resize(GNSS_epoch, 16);
			out_kf_sd(GNSS_epoch - 1, 0) = time;
			for (int i = 1; i <= 15; ++i)
			{
				out_kf_sd(GNSS_epoch - 1, i) = sqrt(P_matrix(i - 1, i - 1));
			}
		}

		// Convert navigation solution to NED
		ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, est_L_b, est_lambda_b, est_h_b, est_v_eb_n, est_C_b_n);
		//cout << est_L_b << endl;
		//cout << est_lambda_b << endl;
		//cout << est_h_b << endl;
		//cout << est_v_eb_n << endl;
		//cout << est_C_b_n << endl;
		// Generate output profile record
		out_profile1(epoch - 1, 0) = time;
		out_profile1(epoch - 1, 1) = est_L_b;
		out_profile1(epoch - 1, 2) = est_lambda_b;
		out_profile1(epoch - 1, 3) = est_h_b;
		out_profile1.row(epoch - 1).segment(4, 3) = est_v_eb_n.transpose();
		out_profile1.row(epoch - 1).segment(7, 3) = CTM_to_Euler(est_C_b_n.transpose()).transpose();

		// Determine errors and generate output record
		Calculate_errors_NED(
			est_L_b, est_lambda_b, est_h_b,
			est_v_eb_n, est_C_b_n,
			true_L_b, true_lambda_b, true_h_b,
			true_v_eb_n, true_C_b_n,
			delta_r_eb_n, delta_v_eb_n, delta_eul_nb_n);

		out_errors(epoch - 1, 0) = time;
		out_errors.row(epoch - 1).segment(1, 3) = delta_r_eb_n.transpose();
		out_errors.row(epoch - 1).segment(4, 3) = delta_v_eb_n.transpose();
		out_errors.row(epoch - 1).segment(7, 3) = delta_eul_nb_n.transpose();

		old_time = time;
		old_est_r_eb_e = est_r_eb_e;
		old_est_v_eb_e = est_v_eb_e;
		old_est_C_b_e = est_C_b_e;
		//cout << setprecision(15) << old_time << endl;
		//cout << setprecision(15) << old_est_r_eb_e << endl;
		//cout << setprecision(15) << old_est_v_eb_e << endl;
		//cout << setprecision(15) << old_est_C_b_e << endl;
		//cout << old_est_C_b_e << endl;
	}

	// Complete progress bar
	cout << rewind << bars << endl;
}