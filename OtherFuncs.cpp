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

VectorXd generateRandomVector()//生成一个包含6个元素的随机向量
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
//仿真IMU输出，在真值的基础上加上误差
void IMU_model(double tor_i, Eigen::Vector3d true_f_ib_b, Eigen::Vector3d true_omega_ib_b, IMU_errors IMU_errors, Eigen::Matrix<double, 6, 1>& old_quant_residuals, Eigen::Vector3d& meas_f_ib_b, Eigen::Vector3d& meas_omega_ib_b, Eigen::Matrix<double, 6, 1>& quant_residuals)
{
	Vector3d accel_noise(0.0, 0.0, 0.0);
	Vector3d gyro_noise(0.0, 0.0, 0.0);

	VectorXd random_vector = generateRandomVector();

	Vector3d accel_distribution_samples(0.0, 0.0, 0.0);//加速度计噪声
	accel_distribution_samples = random_vector.segment(0, 3);


	Vector3d gyro_distribution_samples(0.0, 0.0, 0.0);
	gyro_distribution_samples = random_vector.segment(3, 3);//陀螺仪噪声


	// 输入的时间间隔 tor_i 大于0时，生成加速度计和陀螺仪的噪声项。
	if (tor_i > 0)
	{
		// 加速度计噪声存储在 accel_noise 中，陀螺仪噪声存储在 gyro_noise 中。
		accel_noise = accel_distribution_samples * IMU_errors.accel_noise_root_PSD / sqrt(tor_i);
		gyro_noise = gyro_distribution_samples * IMU_errors.gyro_noise_root_PSD / sqrt(tor_i);
		//模拟随机噪声Wa=random*PSD/sqrt（T）
	
	}
	else
	{
		// 如果时间间隔 tor_i 小于等于0，则将噪声项设置为零向量。
		accel_noise = Vector3d::Zero();
		gyro_noise = Vector3d::Zero();
	}

	Vector3d uq_f_ib_b = IMU_errors.b_a + (Matrix3d::Identity() + IMU_errors.M_a) * true_f_ib_b + accel_noise;



	Vector3d uq_omega_ib_b = IMU_errors.b_g + (Matrix3d::Identity() + IMU_errors.M_g) * true_omega_ib_b +
		IMU_errors.G_g * true_f_ib_b + gyro_noise;


	// 对加速度计和陀螺仪的输出进行量化处理。
	if (IMU_errors.accel_quant_level > 0)
	{
		Vector3d temp_residuals(0.0, 0.0, 0.0);
		Vector3d rounded_residuals(0.0, 0.0, 0.0);
		

		temp_residuals = (uq_f_ib_b + old_quant_residuals.segment(0, 3)) / IMU_errors.accel_quant_level;
		
		for (int i = 0; i < 3; i++)
		{
			rounded_residuals(i) = round(temp_residuals(i));
		}

		;

		meas_f_ib_b = IMU_errors.accel_quant_level * rounded_residuals;

		

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
	

		for (int i = 0; i < 3; i++)
		{
			rounded_residuals(i) = round(temp_residuals(i));
		}

	

		meas_omega_ib_b = IMU_errors.gyro_quant_level * rounded_residuals;

		

		quant_residuals.segment(3, 3) = uq_omega_ib_b + old_quant_residuals.segment(3, 3) - meas_omega_ib_b;
	}
	else
	{
		meas_omega_ib_b = uq_omega_ib_b;
		quant_residuals.segment(3, 3) = Vector3d::Zero();
	}
}


Matrix<double, 15, 15> Initialize_LC_P_matrix(const LC_KF_config LC_KF_config)//初始化P矩阵
{
	Matrix<double, 15, 15> P_matrix = Matrix<double, 15, 15>::Zero();

	// P.block<rows, cols>(i, j)
	// P.block(i, j, rows, cols)
	// 从i行j列开始的rows行cols列块
	//  Initialize error covariance matrix 初始化P矩阵
	P_matrix.block<3, 3>(0, 0) = Matrix3d::Identity() * LC_KF_config.init_att_unc * LC_KF_config.init_att_unc;
	P_matrix.block<3, 3>(3, 3) = Matrix3d::Identity() * LC_KF_config.init_vel_unc * LC_KF_config.init_vel_unc;
	P_matrix.block<3, 3>(6, 6) = Matrix3d::Identity() * LC_KF_config.init_pos_unc * LC_KF_config.init_pos_unc;
	P_matrix.block<3, 3>(9, 9) = Matrix3d::Identity() * LC_KF_config.init_b_a_unc * LC_KF_config.init_b_a_unc;
	P_matrix.block<3, 3>(12, 12) = Matrix3d::Identity() * LC_KF_config.init_b_g_unc * LC_KF_config.init_b_g_unc;

	return P_matrix;
}
//卡尔曼滤波 ，进行组合导航
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
	//状态转移矩阵
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

	Matrix<double, 15, 15> Q_prime_matrix = Matrix<double, 15, 15>::Zero();
	Q_prime_matrix.block(0, 0, 3, 3) = Matrix3d::Identity() * LC_KF_config.gyro_noise_PSD * tor_s;
	Q_prime_matrix.block(3, 3, 3, 3) = Matrix3d::Identity() * LC_KF_config.accel_noise_PSD * tor_s;
	Q_prime_matrix.block(9, 9, 3, 3) = Matrix3d::Identity() * LC_KF_config.accel_bias_PSD * tor_s;
	Q_prime_matrix.block(12, 12, 3, 3) = Matrix3d::Identity() * LC_KF_config.gyro_bias_PSD * tor_s;
	//设置Q矩阵（系统状态方差阵）
	Matrix<double, 15, 1> x_est_propagated = Matrix<double, 15, 1>::Zero();


	Matrix<double, 15, 15> P_matrix_propagated = Matrix<double, 15, 15>::Zero();
	P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
	//Pk（-）当前时刻P矩阵的预测值  

	Matrix<double, 6, 15> H_matrix = MatrixXd::Zero(6, 15);
	H_matrix.block(0, 6, 3, 3) = -Matrix3d::Identity();
	H_matrix.block(3, 3, 3, 3) = -Matrix3d::Identity();
	//观测矩阵 Z为相反数需要取负号

	Matrix<double, 6, 6> R_matrix = Matrix<double, 6, 6>::Zero();
	R_matrix.block(0, 0, 3, 3) = Matrix3d::Identity() * pow(LC_KF_config.pos_meas_SD, 2);
	R_matrix.block(0, 3, 3, 3) = Matrix3d::Zero();
	R_matrix.block(3, 0, 3, 3) = Matrix3d::Zero();
	R_matrix.block(3, 3, 3, 3) = Matrix3d::Identity() * pow(LC_KF_config.vel_meas_SD, 2);
	//量测噪声方差矩阵

	Matrix<double, 15, 6> K_matrix = Matrix<double, 15, 6>::Zero();
	K_matrix = P_matrix_propagated * H_matrix.transpose() * (H_matrix * P_matrix_propagated * H_matrix.transpose() + R_matrix).inverse();
	//增益矩阵（对权重的计算）

	VectorXd delta_z(6);
	delta_z << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	delta_z.segment(0, 3) = GNSS_r_eb_e - est_r_eb_e_old;
	delta_z.segment(3, 3) = GNSS_v_eb_e - est_v_eb_e_old;
	//Z观测向量

	VectorXd x_est_new(15);
	x_est_new = x_est_propagated + K_matrix * delta_z;
	//状态向量更新

	P_matrix_new = (Matrix<double, 15, 15>::Identity() - K_matrix * H_matrix) * P_matrix_propagated;
	//当前时刻P矩阵最优估计Pk（+）

	est_C_b_e_new = (Matrix3d::Identity() - Skew_symmetric(x_est_new.segment(0, 3))) * est_C_b_e_old;
	est_v_eb_e_new = est_v_eb_e_old - x_est_new.segment(3, 3);
	est_r_eb_e_new = est_r_eb_e_old - x_est_new.segment(6, 3);
	//减去误差，更新结果

	// Update IMU bias estimates  更新IMU零偏误差
	est_IMU_bias_new = est_IMU_bias_old + x_est_new.segment(9, 6);
}

//组合导航，松耦合
void Loosely_coupled_INS_GNSS(const MatrixXd& in_profile, int no_epochs, const MatrixXd& in_profile1, int no_epochs1,
	const MatrixXd& in_profile2, int no_epochs2, const InitializationErrors initialization_errors,
	const IMU_errors IMU_errors, const GNSS_config& gnss_config, const LC_KF_config LC_KF_config,
	MatrixXd& out_profile1, MatrixXd& out_errors, MatrixXd& out_imu_bias_est,
	MatrixXd& out_clock, MatrixXd& out_kf_sd)
{

	// Initialize true navigation solution用真值做初始化
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
	//导航系转化到地心地固坐标系old_true_r_eb_e, old_true_v_eb_e, old_true_C_b_e

	Vector2d est_clock;//初始化时间
	est_clock = Vector2d::Zero();

	//GNSS数据est
	double old_est_L_b = in_profile2(0, 1);
	double old_est_lambda_b = in_profile2(0, 2);
	double old_est_h_b = in_profile2(0, 3);
	Vector3d old_est_v_eb_n = in_profile2.row(0).segment(4, 3).transpose();
	

	double est_L_b = old_est_L_b;
	
	//导航系转化到地心地固坐标系
	Vector3d old_est_r_eb_e(0.0, 0.0, 0.0);
	Vector3d old_est_v_eb_e(0.0, 0.0, 0.0);
	Matrix3d C_b_e = Matrix3d::Zero();

	NED_to_ECEF(old_est_r_eb_e, old_est_v_eb_e, C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, Matrix3d::Zero());
	

	Matrix3d old_est_C_b_n = Initialize_NED_attitude(true_C_b_n, initialization_errors.delta_eul_nb_n);
	//真值方向余弦矩阵b系转换为n系，加上初始化偏差等于仿真方向余弦矩阵（无法由GNSS提供）

	Vector3d temp1(0.0, 0.0, 0.0);
	Vector3d temp2(0.0, 0.0, 0.0);
	Matrix3d old_est_C_b_e = Matrix3d::Zero();
	NED_to_ECEF(temp1, temp2, old_est_C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, old_est_C_b_n);

	//old_est_r_eb_e, old_est_v_eb_e,old_est_C_b_e

	// Initialize output profile record and errors record   定义输出文件
	out_profile1 = MatrixXd::Zero(no_epochs1, 10);
	out_errors = MatrixXd::Zero(no_epochs1, 10);

	// Generate output profile record for the first epoch初始化输出文件
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

	//计算结果与真值的误差
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


	VectorXd est_IMU_bias = VectorXd::Zero(6);

	// Initialize IMU quantization residuals IMU量化残差
	Matrix<double, 6, 1> quant_residuals;
	quant_residuals << 0, 0, 0, 0, 0, 0;

	out_imu_bias_est = Matrix<double, 1, 7>::Zero(); //零偏
	out_imu_bias_est(0, 0) = old_time;
	out_imu_bias_est.row(0).segment(1, 6) = est_IMU_bias.transpose();
	out_clock = Matrix<double, 1, 3>::Zero();    //时刻
	out_clock(0, 0) = old_time;
	out_clock.row(0).segment(1, 2) = est_clock;

	// Generate KF uncertainty record卡尔曼滤波不确定度
	out_kf_sd = MatrixXd::Zero(1, 16);
	out_kf_sd(0, 0) = old_time;

	for (int i = 1; i <= 15; ++i)
	{
		out_kf_sd(0, i) = sqrt(P_matrix(i - 1, i - 1));
	}



	// Initialize GNSS model timing
	double time_last_GNSS = old_time;
	int GNSS_epoch = 1;

	// Progress bar 进度条
	string dots = "....................";
	string bars = "||||||||||||||||||||";
	string rewind = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
	cout << "Processing: " << dots;
	int progress_mark = 0;
	int progress_epoch = 0;

	// Main loop
	
	for (int epoch = 2; epoch <= no_epochs1; ++epoch)  //no_epochs1 IMU输出文件响应次数
	{
		// Update progress bar
		if ((epoch - progress_epoch) > (no_epochs1 / 20.0)) 
			//当第一个阶段IMU解算时间超过整个导航时间的1/20时，出现一个进度条，一共有20条
		{
			progress_mark++;
			progress_epoch = epoch;
			cout << rewind << bars.substr(0, progress_mark) << dots.substr(0, 20 - progress_mark);
		}
		

		// Input data from motion profile  读取该时刻的真值
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

		double time = in_profile1(epoch - 1, 0);
		Vector3d alpha_ib_b = in_profile1.row(epoch - 1).segment(1, 3).transpose();
		//角速度增量真值
		Vector3d v_ib_b = in_profile1.row(epoch - 1).segment(4, 3).transpose();
		//速度增量真值
		
		double tor_i = time - old_time;//计算时间间隔
		
		Vector3d meas_f_ib_b_old = v_ib_b / tor_i;
		//平均加速度

		Vector3d meas_omega_ib_b_old = alpha_ib_b / tor_i;
		//平均角加速度

		Vector3d meas_omega_ib_b(0.0, 0.0, 0.0);
		Vector3d meas_f_ib_b(0.0, 0.0, 0.0);

		Matrix<double, 6, 1> quant_residuals_new;
		quant_residuals_new << 0, 0, 0, 0, 0, 0;
		

		IMU_model(tor_i, meas_f_ib_b_old, meas_omega_ib_b_old, IMU_errors, quant_residuals, meas_f_ib_b, meas_omega_ib_b, quant_residuals_new);
	//仿真IMU输出，在真值的基础上加上随机误差，模拟现实场景下IMU输出

		quant_residuals = quant_residuals_new;
		meas_f_ib_b -= est_IMU_bias.segment(0, 3);//进行零偏修正
		meas_omega_ib_b -= est_IMU_bias.segment(3, 3);

		Vector3d est_r_eb_e(0.0, 0.0, 0.0);
		Vector3d est_v_eb_e(0.0, 0.0, 0.0);
		Matrix3d est_C_b_e = Matrix3d::Identity();

		double est_lambda_b = 0;
		double est_h_b = 0;

		//进行惯性导航，力学编排
		Nav_equations_ECEF(tor_i, old_est_r_eb_e, old_est_v_eb_e, old_est_C_b_e, meas_f_ib_b, meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
		

		Vector3d est_v_eb_n(0.0, 0.0, 0.0);
		Matrix3d est_C_b_n = Matrix3d::Identity();
		
		// Determine whether to update GNSS simulation and run Kalman filter
		//IMU响应时间0.005s，GNSS响应时间为1s
		//time为当前IMU历元，old time上一个历元time_last_GNSS=old time
		if ((time - time_last_GNSS) >= gnss_config.epoch_interval)//GNSS响应间隔1s
		{
			GNSS_epoch++;
			double tor_s = time - time_last_GNSS;//上一次GNSS数据进入与现在时间的间隔
			time_last_GNSS = in_profile2(GNSS_epoch - 1, 0);
			//更新GNSS数据加入的时间
			if ((time_last_GNSS - time) > 0.5)
				//如果加入的时间先于现在的时间超过0.5s，则不进入组合导航，很少满足
			{
				GNSS_epoch = GNSS_epoch - 1;//GNSS响应次数减去1
				time_last_GNSS = time_last_GNSS - 1;
				ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, est_L_b, est_lambda_b, est_h_b, est_v_eb_n, est_C_b_n);
		
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

			old_est_L_b = in_profile2(GNSS_epoch - 1, 1);//使用GNSS数据作为现在时刻导航结果
			old_est_lambda_b = in_profile2(GNSS_epoch - 1, 2);
			old_est_h_b = in_profile2(GNSS_epoch - 1, 3);
			
			double a, b, c;
			a = b = c = 0;
			Matrix3d A = Matrix3d::Zero();

			ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, a, b, c, est_v_eb_n, A);
			//通过力学编排求得的est_r_eb_e, est_v_eb_e, est_C_b_e求得est_v_eb_n
			old_est_v_eb_n = est_v_eb_n;

			est_L_b = old_est_L_b;

			Vector3d GNSS_r_eb_e(0.0, 0.0, 0.0);
			Vector3d GNSS_v_eb_e(0.0, 0.0, 0.0);
			Matrix3d GNSS_C_b_e = Matrix3d::Zero();

			NED_to_ECEF(GNSS_r_eb_e, GNSS_v_eb_e, GNSS_C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, Matrix3d::Zero());
			//从n系转换为e系
			Matrix<double, 15, 15> P_matrix_new = Matrix<double, 15, 15>::Zero();

			LC_KF_Epoch(GNSS_r_eb_e, GNSS_v_eb_e, tor_s, est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix,
				meas_f_ib_b, est_L_b, LC_KF_config, est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix_new);
			P_matrix = P_matrix_new;//更新P矩阵


			out_imu_bias_est.resize(GNSS_epoch, 7);
			out_imu_bias_est(GNSS_epoch - 1, 0) = time;//更新零偏
			out_imu_bias_est.row(GNSS_epoch - 1).segment(1, 6) = est_IMU_bias.transpose();

			out_clock.resize(GNSS_epoch, 3);
			out_clock(GNSS_epoch - 1, 0) = time;//更新时间
			out_clock.row(GNSS_epoch - 1).segment(1, 2) = est_clock;

			out_kf_sd.resize(GNSS_epoch, 16);
			out_kf_sd(GNSS_epoch - 1, 0) = time;
			for (int i = 1; i <= 15; ++i)
			{
				out_kf_sd(GNSS_epoch - 1, i) = sqrt(P_matrix(i - 1, i - 1));
			}
		}//更新卡尔曼滤波不确定度

		// Convert navigation solution to NED
		ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, est_L_b, est_lambda_b, est_h_b, est_v_eb_n, est_C_b_n);
	//将得到的e系中的数据转换为n系
		// Generate output profile record 输出结果
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
		old_est_C_b_e = est_C_b_e;    //new转变为old开始下一个循环

	}

	// Complete progress bar
	cout << rewind << bars << endl;
}