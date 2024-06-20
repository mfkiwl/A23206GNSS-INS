#include <iostream>
#include <Eigen/Dense>
#include <iomanip> // ���� setprecision
#include <cmath>
#include <fstream>
#include <random>
#include <chrono>
#include "calculate.h"
#include "OtherFuncs.h"

using namespace std;
using namespace Eigen;

VectorXd generateRandomVector()//����һ������6��Ԫ�ص��������
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
//����IMU���������ֵ�Ļ����ϼ������
void IMU_model(double tor_i, Eigen::Vector3d true_f_ib_b, Eigen::Vector3d true_omega_ib_b, IMU_errors IMU_errors, Eigen::Matrix<double, 6, 1>& old_quant_residuals, Eigen::Vector3d& meas_f_ib_b, Eigen::Vector3d& meas_omega_ib_b, Eigen::Matrix<double, 6, 1>& quant_residuals)
{
	Vector3d accel_noise(0.0, 0.0, 0.0);
	Vector3d gyro_noise(0.0, 0.0, 0.0);

	VectorXd random_vector = generateRandomVector();

	Vector3d accel_distribution_samples(0.0, 0.0, 0.0);//���ٶȼ�����
	accel_distribution_samples = random_vector.segment(0, 3);


	Vector3d gyro_distribution_samples(0.0, 0.0, 0.0);
	gyro_distribution_samples = random_vector.segment(3, 3);//����������


	// �����ʱ���� tor_i ����0ʱ�����ɼ��ٶȼƺ������ǵ������
	if (tor_i > 0)
	{
		// ���ٶȼ������洢�� accel_noise �У������������洢�� gyro_noise �С�
		accel_noise = accel_distribution_samples * IMU_errors.accel_noise_root_PSD / sqrt(tor_i);
		gyro_noise = gyro_distribution_samples * IMU_errors.gyro_noise_root_PSD / sqrt(tor_i);
		//ģ���������Wa=random*PSD/sqrt��T��
	
	}
	else
	{
		// ���ʱ���� tor_i С�ڵ���0��������������Ϊ��������
		accel_noise = Vector3d::Zero();
		gyro_noise = Vector3d::Zero();
	}

	Vector3d uq_f_ib_b = IMU_errors.b_a + (Matrix3d::Identity() + IMU_errors.M_a) * true_f_ib_b + accel_noise;



	Vector3d uq_omega_ib_b = IMU_errors.b_g + (Matrix3d::Identity() + IMU_errors.M_g) * true_omega_ib_b +
		IMU_errors.G_g * true_f_ib_b + gyro_noise;


	// �Լ��ٶȼƺ������ǵ����������������
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

		

		// segment<Size>(Start) ���У�Size ����Ƭ�Ĵ�С��Start ����Ƭ����ʼλ�á�
		// quant_residuals.segment<3>(0)�ȼ���quant_residuals.segment(0, 3)
		quant_residuals.segment(0, 3) = uq_f_ib_b + old_quant_residuals.segment(0, 3) - meas_f_ib_b;
	}
	else
	{
		meas_f_ib_b = uq_f_ib_b;
		// Vector3dʵ������ Matrix<double, 3, 1> �ı�������ˣ��� Vector3d ��ֵ���� Matrix<double, 6, 1> ����ʱ��ʵ�����ǽ� Vector3d ��ֵ���Ƶ�Matrix<double, 6, 1> ��ǰ����Ԫ���У�������Ԫ�ر��ֲ��䡣
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


Matrix<double, 15, 15> Initialize_LC_P_matrix(const LC_KF_config LC_KF_config)//��ʼ��P����
{
	Matrix<double, 15, 15> P_matrix = Matrix<double, 15, 15>::Zero();

	// P.block<rows, cols>(i, j)
	// P.block(i, j, rows, cols)
	// ��i��j�п�ʼ��rows��cols�п�
	//  Initialize error covariance matrix ��ʼ��P����
	P_matrix.block<3, 3>(0, 0) = Matrix3d::Identity() * LC_KF_config.init_att_unc * LC_KF_config.init_att_unc;
	P_matrix.block<3, 3>(3, 3) = Matrix3d::Identity() * LC_KF_config.init_vel_unc * LC_KF_config.init_vel_unc;
	P_matrix.block<3, 3>(6, 6) = Matrix3d::Identity() * LC_KF_config.init_pos_unc * LC_KF_config.init_pos_unc;
	P_matrix.block<3, 3>(9, 9) = Matrix3d::Identity() * LC_KF_config.init_b_a_unc * LC_KF_config.init_b_a_unc;
	P_matrix.block<3, 3>(12, 12) = Matrix3d::Identity() * LC_KF_config.init_b_g_unc * LC_KF_config.init_b_g_unc;

	return P_matrix;
}
//�������˲� ��������ϵ���
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
	//״̬ת�ƾ���
	// matrix.block(i,j,p,q) ��ȡ���СΪ(p,q),��ʼ��(i,j)
	Phi_matrix.block(0, 0, 3, 3) -= Omega_ie * tor_s;
	Phi_matrix.block(0, 12, 3, 3) = est_C_b_e_old * tor_s;
	Phi_matrix.block(3, 0, 3, 3) = -tor_s * Skew_symmetric(est_C_b_e_old * meas_f_ib_b);
	Phi_matrix.block(3, 3, 3, 3) -= 2 * Omega_ie * tor_s;

	double geocentric_radius = R_0 / sqrt(1.0 - e * e * sin(est_L_b_old) * sin(est_L_b_old)) * sqrt(cos(est_L_b_old) * cos(est_L_b_old) + (1.0 - e * e) * (1.0 - e * e) * sin(est_L_b_old) * sin(est_L_b_old));

	// transpose()ת��
	Phi_matrix.block(3, 6, 3, 3) = -tor_s * 2 * Gravity_ECEF(est_r_eb_e_old) /
		geocentric_radius * est_r_eb_e_old.transpose() / sqrt(est_r_eb_e_old.transpose() * est_r_eb_e_old);

	Phi_matrix.block(3, 9, 3, 3) = est_C_b_e_old * tor_s;
	Phi_matrix.block(6, 3, 3, 3) = Matrix3d::Identity() * tor_s;

	Matrix<double, 15, 15> Q_prime_matrix = Matrix<double, 15, 15>::Zero();
	Q_prime_matrix.block(0, 0, 3, 3) = Matrix3d::Identity() * LC_KF_config.gyro_noise_PSD * tor_s;
	Q_prime_matrix.block(3, 3, 3, 3) = Matrix3d::Identity() * LC_KF_config.accel_noise_PSD * tor_s;
	Q_prime_matrix.block(9, 9, 3, 3) = Matrix3d::Identity() * LC_KF_config.accel_bias_PSD * tor_s;
	Q_prime_matrix.block(12, 12, 3, 3) = Matrix3d::Identity() * LC_KF_config.gyro_bias_PSD * tor_s;
	//����Q����ϵͳ״̬������
	Matrix<double, 15, 1> x_est_propagated = Matrix<double, 15, 1>::Zero();


	Matrix<double, 15, 15> P_matrix_propagated = Matrix<double, 15, 15>::Zero();
	P_matrix_propagated = Phi_matrix * (P_matrix_old + 0.5 * Q_prime_matrix) * Phi_matrix.transpose() + 0.5 * Q_prime_matrix;
	//Pk��-����ǰʱ��P�����Ԥ��ֵ  

	Matrix<double, 6, 15> H_matrix = MatrixXd::Zero(6, 15);
	H_matrix.block(0, 6, 3, 3) = -Matrix3d::Identity();
	H_matrix.block(3, 3, 3, 3) = -Matrix3d::Identity();
	//�۲���� ZΪ�෴����Ҫȡ����

	Matrix<double, 6, 6> R_matrix = Matrix<double, 6, 6>::Zero();
	R_matrix.block(0, 0, 3, 3) = Matrix3d::Identity() * pow(LC_KF_config.pos_meas_SD, 2);
	R_matrix.block(0, 3, 3, 3) = Matrix3d::Zero();
	R_matrix.block(3, 0, 3, 3) = Matrix3d::Zero();
	R_matrix.block(3, 3, 3, 3) = Matrix3d::Identity() * pow(LC_KF_config.vel_meas_SD, 2);
	//���������������

	Matrix<double, 15, 6> K_matrix = Matrix<double, 15, 6>::Zero();
	K_matrix = P_matrix_propagated * H_matrix.transpose() * (H_matrix * P_matrix_propagated * H_matrix.transpose() + R_matrix).inverse();
	//������󣨶�Ȩ�صļ��㣩

	VectorXd delta_z(6);
	delta_z << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0;
	delta_z.segment(0, 3) = GNSS_r_eb_e - est_r_eb_e_old;
	delta_z.segment(3, 3) = GNSS_v_eb_e - est_v_eb_e_old;
	//Z�۲�����

	VectorXd x_est_new(15);
	x_est_new = x_est_propagated + K_matrix * delta_z;
	//״̬��������

	P_matrix_new = (Matrix<double, 15, 15>::Identity() - K_matrix * H_matrix) * P_matrix_propagated;
	//��ǰʱ��P�������Ź���Pk��+��

	est_C_b_e_new = (Matrix3d::Identity() - Skew_symmetric(x_est_new.segment(0, 3))) * est_C_b_e_old;
	est_v_eb_e_new = est_v_eb_e_old - x_est_new.segment(3, 3);
	est_r_eb_e_new = est_r_eb_e_old - x_est_new.segment(6, 3);
	//��ȥ�����½��

	// Update IMU bias estimates  ����IMU��ƫ���
	est_IMU_bias_new = est_IMU_bias_old + x_est_new.segment(9, 6);
}

//��ϵ����������
void Loosely_coupled_INS_GNSS(const MatrixXd& in_profile, int no_epochs, const MatrixXd& in_profile1, int no_epochs1,
	const MatrixXd& in_profile2, int no_epochs2, const InitializationErrors initialization_errors,
	const IMU_errors IMU_errors, const GNSS_config& gnss_config, const LC_KF_config LC_KF_config,
	MatrixXd& out_profile1, MatrixXd& out_errors, MatrixXd& out_imu_bias_est,
	MatrixXd& out_clock, MatrixXd& out_kf_sd)
{

	// Initialize true navigation solution����ֵ����ʼ��
	double old_time = in_profile(0, 1);
	double true_L_b = in_profile(0, 2);
	double true_lambda_b = in_profile(0, 3);
	double true_h_b = in_profile(0, 4);
	// ��Eigen�У�������Ĭ���������������
	Vector3d true_v_eb_n = in_profile.row(0).segment(5, 3).transpose();
	Vector3d true_eul_nb = in_profile.row(0).segment(8, 3).transpose();
	Matrix3d true_C_b_n = Euler_to_CTM(true_eul_nb).transpose();

	Vector3d old_true_r_eb_e(0.0, 0.0, 0.0);
	Vector3d old_true_v_eb_e(0.0, 0.0, 0.0);
	Matrix3d old_true_C_b_e = Matrix3d::Zero();

	NED_to_ECEF(old_true_r_eb_e, old_true_v_eb_e, old_true_C_b_e, true_L_b, true_lambda_b, true_h_b, true_v_eb_n, true_C_b_n);
	//����ϵת�������ĵع�����ϵold_true_r_eb_e, old_true_v_eb_e, old_true_C_b_e

	Vector2d est_clock;//��ʼ��ʱ��
	est_clock = Vector2d::Zero();

	//GNSS����est
	double old_est_L_b = in_profile2(0, 1);
	double old_est_lambda_b = in_profile2(0, 2);
	double old_est_h_b = in_profile2(0, 3);
	Vector3d old_est_v_eb_n = in_profile2.row(0).segment(4, 3).transpose();
	

	double est_L_b = old_est_L_b;
	
	//����ϵת�������ĵع�����ϵ
	Vector3d old_est_r_eb_e(0.0, 0.0, 0.0);
	Vector3d old_est_v_eb_e(0.0, 0.0, 0.0);
	Matrix3d C_b_e = Matrix3d::Zero();

	NED_to_ECEF(old_est_r_eb_e, old_est_v_eb_e, C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, Matrix3d::Zero());
	

	Matrix3d old_est_C_b_n = Initialize_NED_attitude(true_C_b_n, initialization_errors.delta_eul_nb_n);
	//��ֵ�������Ҿ���bϵת��Ϊnϵ�����ϳ�ʼ��ƫ����ڷ��淽�����Ҿ����޷���GNSS�ṩ��

	Vector3d temp1(0.0, 0.0, 0.0);
	Vector3d temp2(0.0, 0.0, 0.0);
	Matrix3d old_est_C_b_e = Matrix3d::Zero();
	NED_to_ECEF(temp1, temp2, old_est_C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, old_est_C_b_n);

	//old_est_r_eb_e, old_est_v_eb_e,old_est_C_b_e

	// Initialize output profile record and errors record   ��������ļ�
	out_profile1 = MatrixXd::Zero(no_epochs1, 10);
	out_errors = MatrixXd::Zero(no_epochs1, 10);

	// Generate output profile record for the first epoch��ʼ������ļ�
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

	//����������ֵ�����
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

	// Initialize IMU quantization residuals IMU�����в�
	Matrix<double, 6, 1> quant_residuals;
	quant_residuals << 0, 0, 0, 0, 0, 0;

	out_imu_bias_est = Matrix<double, 1, 7>::Zero(); //��ƫ
	out_imu_bias_est(0, 0) = old_time;
	out_imu_bias_est.row(0).segment(1, 6) = est_IMU_bias.transpose();
	out_clock = Matrix<double, 1, 3>::Zero();    //ʱ��
	out_clock(0, 0) = old_time;
	out_clock.row(0).segment(1, 2) = est_clock;

	// Generate KF uncertainty record�������˲���ȷ����
	out_kf_sd = MatrixXd::Zero(1, 16);
	out_kf_sd(0, 0) = old_time;

	for (int i = 1; i <= 15; ++i)
	{
		out_kf_sd(0, i) = sqrt(P_matrix(i - 1, i - 1));
	}



	// Initialize GNSS model timing
	double time_last_GNSS = old_time;
	int GNSS_epoch = 1;

	// Progress bar ������
	string dots = "....................";
	string bars = "||||||||||||||||||||";
	string rewind = "\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b";
	cout << "Processing: " << dots;
	int progress_mark = 0;
	int progress_epoch = 0;

	// Main loop
	
	for (int epoch = 2; epoch <= no_epochs1; ++epoch)  //no_epochs1 IMU����ļ���Ӧ����
	{
		// Update progress bar
		if ((epoch - progress_epoch) > (no_epochs1 / 20.0)) 
			//����һ���׶�IMU����ʱ�䳬����������ʱ���1/20ʱ������һ����������һ����20��
		{
			progress_mark++;
			progress_epoch = epoch;
			cout << rewind << bars.substr(0, progress_mark) << dots.substr(0, 20 - progress_mark);
		}
		

		// Input data from motion profile  ��ȡ��ʱ�̵���ֵ
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
		//���ٶ�������ֵ
		Vector3d v_ib_b = in_profile1.row(epoch - 1).segment(4, 3).transpose();
		//�ٶ�������ֵ
		
		double tor_i = time - old_time;//����ʱ����
		
		Vector3d meas_f_ib_b_old = v_ib_b / tor_i;
		//ƽ�����ٶ�

		Vector3d meas_omega_ib_b_old = alpha_ib_b / tor_i;
		//ƽ���Ǽ��ٶ�

		Vector3d meas_omega_ib_b(0.0, 0.0, 0.0);
		Vector3d meas_f_ib_b(0.0, 0.0, 0.0);

		Matrix<double, 6, 1> quant_residuals_new;
		quant_residuals_new << 0, 0, 0, 0, 0, 0;
		

		IMU_model(tor_i, meas_f_ib_b_old, meas_omega_ib_b_old, IMU_errors, quant_residuals, meas_f_ib_b, meas_omega_ib_b, quant_residuals_new);
	//����IMU���������ֵ�Ļ����ϼ��������ģ����ʵ������IMU���

		quant_residuals = quant_residuals_new;
		meas_f_ib_b -= est_IMU_bias.segment(0, 3);//������ƫ����
		meas_omega_ib_b -= est_IMU_bias.segment(3, 3);

		Vector3d est_r_eb_e(0.0, 0.0, 0.0);
		Vector3d est_v_eb_e(0.0, 0.0, 0.0);
		Matrix3d est_C_b_e = Matrix3d::Identity();

		double est_lambda_b = 0;
		double est_h_b = 0;

		//���й��Ե�������ѧ����
		Nav_equations_ECEF(tor_i, old_est_r_eb_e, old_est_v_eb_e, old_est_C_b_e, meas_f_ib_b, meas_omega_ib_b, est_r_eb_e, est_v_eb_e, est_C_b_e);
		

		Vector3d est_v_eb_n(0.0, 0.0, 0.0);
		Matrix3d est_C_b_n = Matrix3d::Identity();
		
		// Determine whether to update GNSS simulation and run Kalman filter
		//IMU��Ӧʱ��0.005s��GNSS��Ӧʱ��Ϊ1s
		//timeΪ��ǰIMU��Ԫ��old time��һ����Ԫtime_last_GNSS=old time
		if ((time - time_last_GNSS) >= gnss_config.epoch_interval)//GNSS��Ӧ���1s
		{
			GNSS_epoch++;
			double tor_s = time - time_last_GNSS;//��һ��GNSS���ݽ���������ʱ��ļ��
			time_last_GNSS = in_profile2(GNSS_epoch - 1, 0);
			//����GNSS���ݼ����ʱ��
			if ((time_last_GNSS - time) > 0.5)
				//��������ʱ���������ڵ�ʱ�䳬��0.5s���򲻽�����ϵ�������������
			{
				GNSS_epoch = GNSS_epoch - 1;//GNSS��Ӧ������ȥ1
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

			old_est_L_b = in_profile2(GNSS_epoch - 1, 1);//ʹ��GNSS������Ϊ����ʱ�̵������
			old_est_lambda_b = in_profile2(GNSS_epoch - 1, 2);
			old_est_h_b = in_profile2(GNSS_epoch - 1, 3);
			
			double a, b, c;
			a = b = c = 0;
			Matrix3d A = Matrix3d::Zero();

			ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, a, b, c, est_v_eb_n, A);
			//ͨ����ѧ������õ�est_r_eb_e, est_v_eb_e, est_C_b_e���est_v_eb_n
			old_est_v_eb_n = est_v_eb_n;

			est_L_b = old_est_L_b;

			Vector3d GNSS_r_eb_e(0.0, 0.0, 0.0);
			Vector3d GNSS_v_eb_e(0.0, 0.0, 0.0);
			Matrix3d GNSS_C_b_e = Matrix3d::Zero();

			NED_to_ECEF(GNSS_r_eb_e, GNSS_v_eb_e, GNSS_C_b_e, old_est_L_b, old_est_lambda_b, old_est_h_b, old_est_v_eb_n, Matrix3d::Zero());
			//��nϵת��Ϊeϵ
			Matrix<double, 15, 15> P_matrix_new = Matrix<double, 15, 15>::Zero();

			LC_KF_Epoch(GNSS_r_eb_e, GNSS_v_eb_e, tor_s, est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix,
				meas_f_ib_b, est_L_b, LC_KF_config, est_C_b_e, est_v_eb_e, est_r_eb_e, est_IMU_bias, P_matrix_new);
			P_matrix = P_matrix_new;//����P����


			out_imu_bias_est.resize(GNSS_epoch, 7);
			out_imu_bias_est(GNSS_epoch - 1, 0) = time;//������ƫ
			out_imu_bias_est.row(GNSS_epoch - 1).segment(1, 6) = est_IMU_bias.transpose();

			out_clock.resize(GNSS_epoch, 3);
			out_clock(GNSS_epoch - 1, 0) = time;//����ʱ��
			out_clock.row(GNSS_epoch - 1).segment(1, 2) = est_clock;

			out_kf_sd.resize(GNSS_epoch, 16);
			out_kf_sd(GNSS_epoch - 1, 0) = time;
			for (int i = 1; i <= 15; ++i)
			{
				out_kf_sd(GNSS_epoch - 1, i) = sqrt(P_matrix(i - 1, i - 1));
			}
		}//���¿������˲���ȷ����

		// Convert navigation solution to NED
		ECEF_to_NED(est_r_eb_e, est_v_eb_e, est_C_b_e, est_L_b, est_lambda_b, est_h_b, est_v_eb_n, est_C_b_n);
	//���õ���eϵ�е�����ת��Ϊnϵ
		// Generate output profile record ������
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
		old_est_C_b_e = est_C_b_e;    //newת��Ϊold��ʼ��һ��ѭ��

	}

	// Complete progress bar
	cout << rewind << bars << endl;
}