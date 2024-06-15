#include <Eigen/Dense>
#include <iostream>
#include <iomanip> // 包含 std::setprecision
#include "calculate.h"
#include "FileIO.h"
#include "OtherFuncs.h"

using namespace Eigen;
using namespace std;

// 测试

// // 测试LC_KF_Epoch
// int main()
// {
//     // Test inputs
//     Vector3d GNSS_r_eb_e(333.0, 232.0, 357.0);
//     Vector3d GNSS_v_eb_e(3.0, 7.0, 13.0);
//     double tor_s = 1.0;
//     Matrix3d est_C_b_e_old = Matrix3d::Zero();
//     Vector3d est_v_eb_e_old(2.3, 6.0, 11.0);
//     Vector3d est_r_eb_e_old(311.0, 238.0, 301.0);
//     Eigen::VectorXd est_IMU_bias_old(6);
//     est_IMU_bias_old << 0.2, 0.2, 0.3, 0.4, 0.5, 0.6;
//     Matrix<double, 15, 15> P_matrix_old = Matrix<double, 15, 15>::Identity();
//     Vector3d meas_f_ib_b(0.01, 0.02, 0.02);
//     double est_L_b_old = 0.1;
//     LC_KF_config config;
//     config.gyro_noise_PSD = 2e-6;
//     config.accel_noise_PSD = 1e-3;
//     config.accel_bias_PSD = 5e-8;
//     config.gyro_bias_PSD = 1e-12;
//     config.pos_meas_SD = 0.5;
//     config.vel_meas_SD = 0.1;

//     // Test outputs
//     Matrix3d est_C_b_e_new;
//     Vector3d est_v_eb_e_new;
//     Vector3d est_r_eb_e_new;
//     VectorXd est_IMU_bias_new(6);
//     Matrix<double, 15, 15> P_matrix_new;

//     // Call the function
//     LC_KF_Epoch(GNSS_r_eb_e, GNSS_v_eb_e, tor_s, est_C_b_e_old, est_v_eb_e_old,
//                 est_r_eb_e_old, est_IMU_bias_old, P_matrix_old, meas_f_ib_b,
//                 est_L_b_old, config, est_C_b_e_new, est_v_eb_e_new,
//                 est_r_eb_e_new, est_IMU_bias_new, P_matrix_new);

//     // Print the results
//     cout << "est_C_b_e_new:\n"
//          << est_C_b_e_new << endl;
//     cout << "est_v_eb_e_new:\n"
//          << est_v_eb_e_new << endl;
//     cout << "est_r_eb_e_new:\n"
//          << est_r_eb_e_new << endl;
//     cout << "est_IMU_bias_new:\n"
//          << est_IMU_bias_new << endl;
//     cout << "P_matrix_new:\n"
//          << P_matrix_new << endl;

//     return 0;
// }

//// 测试IMU_model
//int main()
//{
//    // 定义输入参数
//    double tor_i = 0.1;
//    Vector3d true_f_ib_b(1, 2, 3);
//    Vector3d true_omega_ib_b(0.1, 0.2, 0.3);
//    IMU_errors IMU_errors;
//    IMU_errors.b_a = Vector3d(0.01, 0.02, 0.03);
//    IMU_errors.b_g = Vector3d(0.001, 0.002, 0.003);
//    IMU_errors.M_a << 1.1, 0.2, 0.3,
//        0.4, 1.5, 0.6,
//        0.7, 0.8, 1.9;
//
//    IMU_errors.M_g << 2.1, 0.2, 0.3,
//        0.4, 2.5, 0.6,
//        0.7, 0.8, 2.9;
//
//    IMU_errors.G_g << 0.001, 0.002, 0.003,
//        0.004, 0.005, 0.006,
//        0.007, 0.008, 0.009;
//    IMU_errors.accel_noise_root_PSD = 0.001;
//    IMU_errors.gyro_noise_root_PSD = 0.0001;
//    IMU_errors.accel_quant_level = 0.01;
//    IMU_errors.gyro_quant_level = 0.001;
//    Matrix<double, 6, 1> quant_residuals = Matrix<double, 6, 1>::Zero();
//    Matrix<double, 6, 1> quant_residuals_new = Matrix<double, 6, 1>::Zero();
//    Vector3d meas_f_ib_b, meas_omega_ib_b;
//
//    // 调用 IMU_model 函数
//    IMU_model(tor_i, true_f_ib_b, true_omega_ib_b, IMU_errors, quant_residuals, meas_f_ib_b, meas_omega_ib_b, quant_residuals_new);
//
//    quant_residuals = quant_residuals_new;
//
//    // 显示结果
//    cout << "meas_f_ib_b: " << endl
//         << meas_f_ib_b << endl;
//    cout << "meas_omega_ib_b: " << endl
//         << meas_omega_ib_b << endl;
//    cout << "quant_residuals: " << endl
//         << quant_residuals << endl;
//
//    return 0;
//}

// //测试generateRandomVector
// int main()
// {
//     Eigen::VectorXd random_vector = generateRandomVector();

//     std::cout << "Standard normal distribution (mean: 0, variance: 1):\n";
//     std::cout << random_vector << std::endl;

//     return 0;
// }

// //测试Calculate_errors_NED
// int main()
// {
//     // 输入参数
//     double est_L_b = 37.7749; // 估计的纬度（弧度）
//     double est_lambda_b = -122.4194; // 估计的经度（弧度）
//     double est_h_b = 0; // 估计的高度（米）
//     Vector3d est_v_eb_n(0, 0, 0); // 相对于ECEF框架的机体框架的速度，沿北、东和下方向（米/秒）
//     Matrix3d est_C_b_n = Matrix3d::Identity(); // 估计的机体到NED坐标变换矩阵

//     double true_L_b = 37.7749; // 真实纬度（弧度）
//     double true_lambda_b = -122.4194; // 真实经度（弧度）
//     double true_h_b = 0; // 真实高度（米）
//     Vector3d true_v_eb_n(0, 0, 0); // 真实相对于ECEF框架的机体框架的速度，沿北、东和下方向（米/秒）
//     Matrix3d true_C_b_n = Matrix3d::Identity(); // 真实的机体到NED坐标变换矩阵

//     // 输出参数
//     Vector3d delta_r_eb_n, delta_v_eb_n, delta_eul_nb_n;

//     // 调用函数
//     Calculate_errors_NED(est_L_b, est_lambda_b, est_h_b, est_v_eb_n, est_C_b_n,
//                          true_L_b, true_lambda_b, true_h_b, true_v_eb_n, true_C_b_n,
//                          delta_r_eb_n, delta_v_eb_n, delta_eul_nb_n);

//     // 显示结果
//     std::cout << "delta_r_eb_n: " << delta_r_eb_n.transpose() << std::endl;
//     std::cout << "delta_v_eb_n: " << delta_v_eb_n.transpose() << std::endl;
//     std::cout << "delta_eul_nb_n: " << delta_eul_nb_n.transpose() << std::endl;

//     return 0;
// }

// //输出
// delta_r_eb_n:  0  0 -0
// delta_v_eb_n: 0 0 0
// delta_eul_nb_n: -0  0 -0

// //测试NED_to_ECEF
// int main()
// {
//     // 输入参数
//     double L_b = 37.7749 * M_PI / 180.0;        // 纬度（弧度）
//     double lambda_b = -122.4194 * M_PI / 180.0; // 经度（弧度）
//     double h_b = 0;                             // 高度（米）
//     Vector3d v_eb_n(0, 0, 0);                   // 相对于ECEF框架的机体框架的速度，沿北、东和下方向（米/秒）
//     Matrix3d C_b_n = Matrix3d::Identity();      // 机体到NED坐标变换矩阵

//     // 输出变量
//     Vector3d r_eb_e; // ECEF框架中的机体框架的笛卡尔位置（米）
//     Vector3d v_eb_e; // ECEF框架中的机体框架的速度（米/秒）
//     Matrix3d C_b_e;  // 机体到ECEF框架的坐标变换矩阵

//     // 调用函数
//     NED_to_ECEF(r_eb_e, v_eb_e, C_b_e, L_b, lambda_b, h_b, v_eb_n, C_b_n);

//     // 显示结果
//     std::cout << "r_eb_e:" << std::endl
//               << r_eb_e << std::endl;
//     std::cout << "v_eb_e:" << std::endl
//               << v_eb_e << std::endl;
//     std::cout << "C_b_e:" << std::endl
//               << C_b_e << std::endl;

//     return 0;
// }

// //输出
// r_eb_e:
// -2.70617e+006
// -4.26106e+006
//  3.88573e+006
// v_eb_e:
// 0
// 0
// 0
// C_b_e:
//  0.328402  0.844146  0.423756
//  0.517091 -0.536113  0.667233
//  0.790423         0 -0.612561

// // 测试Nav_equations_ECEF
// int main()
// {
//     // 测试样例
//     double tor_i = 1;                          // 时间间隔
//     Vector3d old_r_eb_e(0, 0, 0);              // 之前的位置
//     Matrix3d old_C_b_e = Matrix3d::Identity(); // 之前的姿态矩阵
//     Vector3d old_v_eb_e(0, 0, 0);              // 之前的速度
//     Vector3d f_ib_b(0, 0, 0);                  // 体坐标系相对于ECEF坐标系的比力
//     Vector3d omega_ib_b(0, 0, 0);              // 体坐标系相对于ECEF坐标系的角速度

//     // 调用函数
//     Vector3d r_eb_e;
//     Vector3d v_eb_e;
//     Matrix3d C_b_e;
//     Nav_equations_ECEF(tor_i, old_r_eb_e, old_v_eb_e, old_C_b_e, f_ib_b, omega_ib_b, r_eb_e, v_eb_e, C_b_e);

//     // 显示结果
//     std::cout << "r_eb_e:\n"
//               << r_eb_e << std::endl;
//     std::cout << "v_eb_e:\n"
//               << v_eb_e << std::endl;
//     std::cout << "C_b_e:\n"
//               << C_b_e << std::endl;

//     return 0;
// }
// 输出结果
// r_eb_e:
// 0
// 0
// 0
// v_eb_e:
// 0
// 0
// 0
// C_b_e:
//             1  7.29211e-005             0
// -7.29211e-005             1             0
//             0             0             1

// // 测试Read_profile
// int main()
// {
//     // 测试样例
//     string input_profile_name = "truth.nav";
//     string input_profile_name1 = "ADIS16465.txt";
//     string input_profile_name2 = "GNSS_RTK.pos";
//     // string input_profile_name = "test.txt";
//     // string input_profile_name1 = "file1.txt";
//     // string input_profile_name2 = "file2.txt";
//     MatrixXd in_profile, in_profile1, in_profile2;
//     int no_epochs, no_epochs1, no_epochs2;

//     bool success = Read_profile(input_profile_name, input_profile_name1, input_profile_name2,
//                                 in_profile, no_epochs, in_profile1, no_epochs1, in_profile2, no_epochs2);

//     if (success)
//     {
//         cout << "读取文件成功" << endl;
//         // 在这里可以对读取的数据进行进一步处理或使用
//     }
//     else
//     {
//         cout << "读取文件失败" << endl;
//     }

//     // cout << in_profile << endl;
//     // cout << in_profile1 << endl;
//     // cout << in_profile2 << endl;
//     cout << no_epochs << endl;
//     cout << no_epochs1 << endl;
//     cout << no_epochs2 << endl;

//     return 0;
// }

// // 测试样例Radii_of_curvature
// int main()
// {
//     // Test case 1
//     double latitude1 = 45; // Example latitude value
//     double R_N1, R_E1;
//     Radii_of_curvature(latitude1, R_N1, R_E1);
//     std::cout << "Test case 1:" << std::endl;
//     std::cout << "Latitude: " << latitude1 << std::endl;
//     std::cout << std::setprecision(15) << "R_N: " << R_N1 << std::endl;
//     std::cout << std::setprecision(15) << "R_E: " << R_E1 << std::endl;

//     // Test case 2
//     double latitude2 = 30; // Example latitude value
//     double R_N2, R_E2;
//     Radii_of_curvature(latitude2, R_N2, R_E2);
//     std::cout << "Test case 2:" << std::endl;
//     std::cout << "Latitude: " << latitude2 << std::endl;
//     std::cout << "R_N: " << R_N2 << std::endl;
//     std::cout << "R_E: " << R_E2 << std::endl;

//     return 0;
// }
// // 输出结果
// // 6378137
// // 0.0818191908425
// // Test case 1:
// // Latitude: 45
// // R_N: 6381781.58646704
// // R_E: 6393650.76230677
// // 6378137
// // 0.0818191908425
// // Test case 2:
// // Latitude: 30
// // R_N: 6398054.61169313
// // R_E: 6399080.57994098

// //测试Initialize_LC_P_matrix
//  int main()
//  {
//      // 定义初始化配置结构体
//      LC_KF_config config;
//      config.init_att_unc = 0.01;
//      config.init_vel_unc = 0.1;
//      config.init_pos_unc = 1;
//      config.init_b_a_unc = 0.01;
//      config.init_b_g_unc = 0.001;

//     // 调用函数初始化状态估计误差协方差矩阵
//     Matrix<double, 15, 15> P_matrix = Initialize_LC_P_matrix(config);

//     // 打印初始化得到的状态估计误差协方差矩阵
//     std::cout << "P_matrix:" << std::endl;
//     std::cout << P_matrix << std::endl;

//     return 0;
// }

// // P_matrix:
// // 0.0001      0      0      0      0      0      0      0      0      0      0      0      0      0      0
// //      0 0.0001      0      0      0      0      0      0      0      0      0      0      0      0      0
// //      0      0 0.0001      0      0      0      0      0      0      0      0      0      0      0      0
// //      0      0      0   0.01      0      0      0      0      0      0      0      0      0      0      0
// //      0      0      0      0   0.01      0      0      0      0      0      0      0      0      0      0
// //      0      0      0      0      0   0.01      0      0      0      0      0      0      0      0      0
// //      0      0      0      0      0      0      1      0      0      0      0      0      0      0      0
// //      0      0      0      0      0      0      0      1      0      0      0      0      0      0      0
// //      0      0      0      0      0      0      0      0      1      0      0      0      0      0      0
// //      0      0      0      0      0      0      0      0      0 0.0001      0      0      0      0      0
// //      0      0      0      0      0      0      0      0      0      0 0.0001      0      0      0      0
// //      0      0      0      0      0      0      0      0      0      0      0 0.0001      0      0      0
// //      0      0      0      0      0      0      0      0      0      0      0      0 1e-006      0      0
// //      0      0      0      0      0      0      0      0      0      0      0      0      0 1e-006      0
// //      0      0      0      0      0      0      0      0      0      0      0      0      0      0 1e-006

// Gravity_ECEF
//  int main()
//  {
//      Vector3d r_eb_e(5000, 1000, 200);
//      Vector3d g = Gravity_ECEF(r_eb_e);

//     cout << "Gravity Acceleration (ECEF):" << endl;
//     cout << g << endl;

//     return 0;
// }

// // Gravity Acceleration (ECEF):
// // -3.77735e+010
// //  -7.5547e+009
// //   -4.555e+009

// 测试fileToMatrix
//  测试样例
//  int main()
//  {
//      string input_profile_name = "test.txt";
//      Eigen::MatrixXd in_profile;
//      int no_epochs;

//     fileToMatrix(input_profile_name, 10, in_profile, no_epochs);
//     cout << input_profile_name << endl;
//     cout << no_epochs << endl;
//     cout << in_profile << endl;
//     return 0;
// }

// 测试Euler_to_CTM
//  int main()
//  {
//      // 输入欧拉角
//      Vector3d eul(0.5, 0.3, 0.2);

//     // 调用函数计算坐标变换矩阵
//     Matrix3d C = Euler_to_CTM(eul);

//     // 打印计算得到的坐标变换矩阵
//     cout << "Coordinate Transformation Matrix (CTM):" << endl;
//     cout << C << endl;

//     return 0;
// }
// //输出
// // Coordinate Transformation Matrix (CTM):
// //  0.936293  0.189796  -0.29552
// // -0.035493  0.888237  0.458013
// //  0.349421 -0.418345  0.838387

// ECEF_to_NED
//  //测试样例
//  int main()
//  {
//      Vector3d r_eb_e1(3967283, 337505, 4966824);
//      Vector3d v_eb_e1(0, 0, 0);
//      Matrix3d C_b_e1 = Matrix3d::Identity();
//      double L_b1, lambda_b1, h_b1;
//      Vector3d v_eb_n1;
//      Matrix3d C_b_n1;

//     ECEF_to_NED(r_eb_e1, v_eb_e1, C_b_e1, L_b1, lambda_b1, h_b1, v_eb_n1, C_b_n1);

//     cout << "测试样例 1:" << endl;
//     cout << "r_eb_e:" << endl
//          << r_eb_e1 << endl;
//     cout << "v_eb_e:" << endl
//          << v_eb_e1 << endl;
//     cout << "C_b_e:" << endl
//          << C_b_e1 << endl;
//     cout << "L_b: " << L_b1 << endl;
//     cout << "lambda_b: " << lambda_b1 << endl;
//     cout << "h_b: " << h_b1 << endl;
//     cout << "v_eb_n:" << endl
//          << v_eb_n1 << endl;
//     cout << "C_b_n:" << endl
//          << C_b_n1 << endl;

//     Vector3d r_eb_e2(4567890, -1234567, 7890123);
//     Vector3d v_eb_e2(100, -200, 300);
//     Matrix3d C_b_e2;
//     C_b_e2 << 0.866, -0.500, 0.000,
//         0.500, 0.866, 0.000,
//         0.000, 0.000, 1.000;
//     double L_b2, lambda_b2, h_b2;
//     Vector3d v_eb_n2;
//     Matrix3d C_b_n2;

//     ECEF_to_NED(r_eb_e2, v_eb_e2, C_b_e2, L_b2, lambda_b2, h_b2, v_eb_n2, C_b_n2);

//     cout << "测试样例 2:" << endl;
//     cout << "r_eb_e:" << endl
//          << r_eb_e2 << endl;
//     cout << "v_eb_e:" << endl
//          << v_eb_e2 << endl;
//     cout << "C_b_e:" << endl
//          << C_b_e2 << endl;
//     cout << std::setprecision(15)<< "L_b: " << L_b2 << endl;
//     cout << std::setprecision(15)<< "lambda_b: " << lambda_b2 << endl;
//     cout << std::setprecision(15)<< "h_b: " << h_b2 << endl;
//     cout << std::setprecision(15) << "v_eb_n:" << endl
//          << v_eb_n2 << endl;
//     cout << std::setprecision(15) << "C_b_n:" << endl
//          << C_b_n2 << endl;

//     return 0;
// }

// 测试样例 1:
// r_eb_e:
// 3.96728e+006
//       337505
// 4.96682e+006
// v_eb_e:
// 0
// 0
// 0
// C_b_e:
// 1 0 0
// 0 1 0
// 0 0 1
// L_b: 0.89833
// lambda_b: 0.0848677
// h_b: 642.411
// v_eb_n:
//  0
//  0
// -0
// C_b_n:
//  -0.779472 -0.0663113   0.622917
// -0.0847659   0.996401          0
//  -0.620675 -0.0528021  -0.782288
// 测试样例 2:
// r_eb_e:
//  4.56789e+006
// -1.23457e+006
//  7.89012e+006
// v_eb_e:
//  100
// -200
//  300
// C_b_e:
// 0.866  -0.5     0
//   0.5 0.866     0
//     0     0     1
// L_b: 1.03263787271023
// lambda_b: -0.263964159229175
// h_b: 2837810.76986713
// v_eb_n:
//  26.0691866073423
//  -166.98171130922
// -333.822566040518
// C_b_n:
// -0.605823246400733  0.608467336268173  0.512555610834937
//  0.708629219627819  0.705549877109811                  0
// -0.361633548236534  0.363211882521821 -0.858654031494423

// 测试CTM_to_Euler
//  测试样例
//  控制台输出
//  0.866   0.5     0
//   -0.5 0.866     0
//      0     0     1
//         0
//        -0
//  0.523611

//     1     0     0
//     0 0.866  -0.5
//     0   0.5 0.866
// -0.523611
//        -0
//         0

// int main()
// {
//     // Matrix3d C;
//     // C << 0.866, 0.500, 0.000,
//     //     -0.500, 0.866, 0.000,
//     //     0.000, 0.000, 1.000;
//     // Vector3d eul = CTM_to_Euler(C);
//     // std::cout << C << std::endl;
//     // std::cout << eul << std::endl;
//     // return 0;

//     Matrix3d C;
//     C << 1.000, 0.000, 0.000,
//         0.000, 0.866, -0.500,
//         0.000, 0.500, 0.866;
//     Vector3d eul = CTM_to_Euler(C);
//     std::cout << C << std::endl;
//     std::cout << eul << std::endl;
//     return 0;
// }

// // 测试Skew_symmetric
// int main()
// {
//     // 测试样例
//     Vector3d a(1, 2, 3); // 输入向量 a

//     Matrix3d A_actual = Skew_symmetric(a); // 调用函数计算输出矩阵 A

//     cout << A_actual << endl;
//     return 0;
// }
// // 输出
// //  0 -3  2
// //  3  0 -1
// //  -2  1  0

// int main()
// {
// 	printf("Hello World!");
// 	return 0;
// }
// 测试样例，在Eigen中，向量的默认情况下是列向量，不需要进行转置
// int main()
// {
//     Matrix3d a;
//     a << 0, 1, 2, 3, 4, 5, 6, 7, 8;
//     // cout << a(0, 0) << endl;
//     // cout << a(0, 2) << endl;
//     // cout << a(2, 2) << endl;
//     // Assertion failed!
//     // cout << a(0, 3) << endl;
//     Vector3d b = a.row(0).segment(0, 3);
//     cout << b << endl;
// }

// 测试Initialize_NED_attitude
//  int main()
//  {
//      Eigen::Matrix3d Euler_to_CTM(Eigen::Vector3d eul);

//     Matrix3d C_b_n = Matrix3d::Identity();

//     Vector3d delta_eul_nb_n(0.1, -0.2, 0.3);

//     Matrix3d est_C_b_n = Initialize_NED_attitude(C_b_n, delta_eul_nb_n);

//     std::cout << "est_C_b_n: " << std::endl;
//     std::cout << est_C_b_n << std::endl;

//     return 0;
// }