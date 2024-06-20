#include <iostream>
#include <vector>
#include <iomanip> // 包含头文件以使用 std::setprecision
#include <cmath>
#include <fstream>
#include "FileIO.h"
#include "OtherFuncs.h"

using namespace std;
using namespace Eigen;

const double rad_to_deg = 1.0 / deg_to_rad;
const double micro_g_to_meters_per_second_squared = 9.80665E-6;

int main()
{
	string input_profile_name = "truth.nav";//真值
	//string input_profile_name = "truth1.nav";

	string input_profile_name1 = "ADIS16465.txt";//IMU输出
	//string input_profile_name1 = "ADIS.txt";

	 string input_profile_name2 = "GNSS_RTK.pos";//高精度RTK数据
	//string input_profile_name2 = "GNSS_RTK1.pos";

	string output_profile_name = "INS_GNSS_Demo_3_Profile.csv";
	string output_errors_name = "INS_GNSS_Demo_3_Errors.csv";

	// Attitude initialization error
	InitializationErrors initialization_errors = {};
	initialization_errors.delta_eul_nb_n = Vector3d(-0.05 * deg_to_rad, 0.04 * deg_to_rad, 1.0 * deg_to_rad);

	// IMU errors
	IMU_errors imu_errors = {};
	imu_errors.b_a = Vector3d(600 * micro_g_to_meters_per_second_squared, 600 * micro_g_to_meters_per_second_squared, 600 * micro_g_to_meters_per_second_squared);
	imu_errors.b_g = Vector3d(75 * deg_to_rad / 3600, 75 * deg_to_rad / 3600, 75 * deg_to_rad / 3600);
	//零偏误差
	imu_errors.M_a << 500E-6, -300E-6, 200E-6,
		-150E-6, -600E-6, 250E-6,
		-250E-6, 100E-6, 450E-6;
	imu_errors.M_g << 400E-6, -300E-6, 250E-6,
		0, -300E-6, -150E-6,
		0, 0, -350E-6;
	//比例因子误差
	imu_errors.G_g << 0.9 * deg_to_rad / (3600 * 9.80665), -1.1 * deg_to_rad / (3600 * 9.80665), -0.6 * deg_to_rad / (3600 * 9.80665),
		-0.5 * deg_to_rad / (3600 * 9.80665), 1.9 * deg_to_rad / (3600 * 9.80665), -1.6 * deg_to_rad / (3600 * 9.80665),
		0.3 * deg_to_rad / (3600 * 9.80665), 1.1 * deg_to_rad / (3600 * 9.80665), -1.3 * deg_to_rad / (3600 * 9.80665);
	//交轴耦合误差误差
	imu_errors.accel_noise_root_PSD = 100 * micro_g_to_meters_per_second_squared;//加速度计随机噪声
	imu_errors.gyro_noise_root_PSD = 0.01 * deg_to_rad / 60;//陀螺仪随机噪声
	imu_errors.accel_quant_level = 1E-2;
	imu_errors.gyro_quant_level = 2E-4;

	// GNSS configuration GNSS相关参数
	GNSS_config gnss_config = {};
	gnss_config.epoch_interval = 1;
	gnss_config.init_est_r_ea_e = Vector3d(0, 0, 0);
	gnss_config.no_sat = 30;
	gnss_config.r_os = 2.656175E7;
	gnss_config.inclination = 55;
	gnss_config.const_delta_lambda = 0;
	gnss_config.const_delta_t = 0;
	gnss_config.mask_angle = 10;
	gnss_config.SIS_err_SD = 1;
	gnss_config.zenith_iono_err_SD = 2;
	gnss_config.zenith_trop_err_SD = 0.2;
	gnss_config.code_track_err_SD = 1;
	gnss_config.rate_track_err_SD = 0.02;
	gnss_config.rx_clock_offset = 10000;
	gnss_config.rx_clock_drift = 100;

	// Kalman filter configuration 卡尔曼滤波器配置
	LC_KF_config kf_config = {};
	kf_config.init_att_unc = deg_to_rad;
	kf_config.init_vel_unc = 0.1;
	kf_config.init_pos_unc = 0.1;
	kf_config.init_b_a_unc = 200 * micro_g_to_meters_per_second_squared;
	kf_config.init_b_g_unc = 25 * deg_to_rad / 3600;
	kf_config.gyro_noise_PSD = pow(0.1 * deg_to_rad / 60, 2);
	kf_config.accel_noise_PSD = pow(200 * micro_g_to_meters_per_second_squared, 2);
	kf_config.accel_bias_PSD = 1.0E-7;
	kf_config.gyro_bias_PSD = 2.0E-12;
	kf_config.pos_meas_SD = 0.01;
	kf_config.vel_meas_SD = 1;

	// Profile data
	int no_epochs, no_epochs1, no_epochs2;//Read_profile中进行初始化
	MatrixXd in_profile, in_profile1, in_profile2;
	in_profile.resize(1, 11);//2 时间 3-4 经纬高 6-8 速度 9-11 姿态
	in_profile1.resize(1, 7);//1时间 2—4角度增量 5—7加速度
	in_profile2.resize(1, 7);//1时间 2-4经纬高 5-7速度
	in_profile.setZero();//初始化
	in_profile1.setZero();
	in_profile2.setZero();

	bool read_profile_success = Read_profile(input_profile_name, input_profile_name1, input_profile_name2,
		in_profile, no_epochs, in_profile1, no_epochs1, in_profile2, no_epochs2);

	if (!read_profile_success)
	{
		cout << "Error reading the input profile." << endl;
		return 1;
	}

	// Output data
	MatrixXd out_profile1, out_errors, out_IMU_bias_est, out_clock, out_KF_SD;

	Loosely_coupled_INS_GNSS(in_profile, no_epochs,
		in_profile1, no_epochs1,
		in_profile2, no_epochs2, initialization_errors, imu_errors, gnss_config, kf_config,
		out_profile1, out_errors, out_IMU_bias_est, out_clock, out_KF_SD);

	//cout << setprecision(15) << out_profile1 << endl;

	ofstream fout("out_profile1.txt", ios::binary);
	fout << setprecision(15) << out_profile1 << endl;
	fout.flush(); //创建文件out_profile1.txt
	return 0;
}