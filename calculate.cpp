#define _USE_MATH_DEFINES
#include <iostream>
#include <Eigen/Dense>
#include <iomanip> // 包含 std::setprecision
#include <cmath>
#include <fstream>
#include <random>
#include "calculate.h"

using namespace std;
using namespace Eigen;

void Calculate_errors_NED(
    double est_L_b, double est_lambda_b, double est_h_b,
    Vector3d est_v_eb_n, Matrix3d est_C_b_n,
    double true_L_b, double true_lambda_b, double true_h_b,
    Vector3d true_v_eb_n, Matrix3d true_C_b_n,
    Vector3d& delta_r_eb_n, Vector3d& delta_v_eb_n, Vector3d& delta_eul_nb_n)
{

    // Calculate position error, using (2.119)
    double R_N, R_E;
    // 计算给定纬度下的曲率半径
    Radii_of_curvature(true_L_b, R_N, R_E);

    // 位置误差
    delta_r_eb_n(0) = (est_L_b - true_L_b) * (R_N + true_h_b);
    delta_r_eb_n(1) = (est_lambda_b - true_lambda_b) * (R_E + true_h_b) * cos(true_L_b);
    delta_r_eb_n(2) = -(est_h_b - true_h_b);

    // Calculate velocity error
    // 速度误差
    delta_v_eb_n = est_v_eb_n - true_v_eb_n;

    // Calculate attitude error, using (5.109) and (5.111)
    // 姿态误差
    Matrix3d delta_C_b_n = est_C_b_n * true_C_b_n.transpose();
    delta_eul_nb_n = -CTM_to_Euler(delta_C_b_n);
}

Vector3d CTM_to_Euler(Matrix3d C)
{
    Vector3d eul;

    // Calculate Euler angles using (2.23)
    eul(0) = atan2(C(1, 2), C(2, 2));   // roll滚转
    eul(1) = -asin(C(0, 2));            // pitch俯仰
    double a = atan2(C(0, 1), C(0, 0)); // yaw偏航
    if (a < 0)
    {
        eul(2) = a + 2 * M_PI;
    }
    else
    {
        eul(2) = a;
    }

    return eul;
}

void ECEF_to_NED(Vector3d r_eb_e, Vector3d v_eb_e, Matrix3d C_b_e, double& L_b, double& lambda_b, double& h_b, Vector3d& v_eb_n, Matrix3d& C_b_n)
{

    double R_0 = 6378137;
    double e = 0.0818191908425;

    lambda_b = atan2(r_eb_e(1), r_eb_e(0));

    // double sqrt(double x) 返回 x 的平方根
    // int abs(int x) 返回整数 x 的绝对值
    double k1 = sqrt(1.0 - e * e) * abs(r_eb_e(2));
    double k2 = e * e * R_0;
    double beta = sqrt(r_eb_e(0) * r_eb_e(0) + r_eb_e(1) * r_eb_e(1));
    double E = (k1 - k2) / beta;
    double F = (k1 + k2) / beta;

    double P = 4.0 / 3.0 * (E * F + 1.0);

    double Q = 2.0 * (E * E - F * F);

    double D = P * P * P + Q * Q;

    // double pow(double x, double y) 返回 x 的 y 次幂
    double V = pow((sqrt(D) - Q), 1.0 / 3.0) - pow((sqrt(D) + Q), 1.0 / 3.0);

    double G = 0.5 * (sqrt(E * E + V) + E);

    double T = sqrt(G * G + (F - V * G) / (2.0 * G - E)) - G;

    // copysign()函数接受两个参数，并返回一个值，该值具有第一个参数的大小和第二个参数的符号
    // 能否确定值的范围？
    // double fabs(double x)，求绝对值
    L_b = copysign(atan((1.0 - T * T) / (2.0 * T * sqrt(1.0 - e * e))), r_eb_e(2));

    // R_0 * sqrt(1 - e * e)大于0
    h_b = (beta - R_0 * T) * cos(L_b) + (r_eb_e(2) - copysign(R_0 * sqrt(1.0 - e * e), r_eb_e(2))) * sin(L_b);

    double cos_lat = cos(L_b);
    double sin_lat = sin(L_b);
    double cos_long = cos(lambda_b);
    double sin_long = sin(lambda_b);
    Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
        -sin_long, cos_long, 0,
        -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;

    v_eb_n = C_e_n * v_eb_e;

    C_b_n = C_e_n * C_b_e;
}

Matrix3d Euler_to_CTM(Vector3d eul)
{

    double sin_phi = sin(eul(0));
    double cos_phi = cos(eul(0));

    double sin_theta = sin(eul(1));
    double cos_theta = cos(eul(1));

    double sin_psi = sin(eul(2));
    double cos_psi = cos(eul(2));

    Matrix3d C;

    C(0, 0) = cos_theta * cos_psi;
    C(0, 1) = cos_theta * sin_psi;
    C(0, 2) = -sin_theta;

    C(1, 0) = -cos_phi * sin_psi + sin_phi * sin_theta * cos_psi;
    C(1, 1) = cos_phi * cos_psi + sin_phi * sin_theta * sin_psi;
    C(1, 2) = sin_phi * cos_theta;

    C(2, 0) = sin_phi * sin_psi + cos_phi * sin_theta * cos_psi;
    C(2, 1) = -sin_phi * cos_psi + cos_phi * sin_theta * sin_psi;
    C(2, 2) = cos_phi * cos_theta;

    return C;
}

Vector3d Gravity_ECEF(Vector3d r_eb_e)
{

    double R_0 = 6378137;
    double mu = 3.986004418E14;
    double J_2 = 1.082627E-3;
    double omega_ie = 7.292115E-5;

    // norm()计算向量的模
    double mag_r = r_eb_e.norm();

    Vector3d g;

    // If the input position is 0,0,0, produce a dummy output
    if (mag_r == 0)
    {
        g << 0, 0, 0;
    }
    else
    {
        // Calculate gravitational acceleration using (2.142)
        double z_scale = 5 * pow(r_eb_e(2) / mag_r, 2);
        Vector3d gamma = -mu / pow(mag_r, 3) * (r_eb_e + 1.5 * J_2 * pow(R_0 / mag_r, 2) * Vector3d((1 - z_scale) * r_eb_e(0), (1 - z_scale) * r_eb_e(1), (3 - z_scale) * r_eb_e(2)));

        // Add centripetal acceleration using (2.133)
        g.head(2) = gamma.head(2) + pow(omega_ie, 2) * r_eb_e.head(2);
        g(2) = gamma(2);
    }

    return g;
}

void Nav_equations_ECEF(const double tor_i,
    const Vector3d& old_r_eb_e,
    const Vector3d& old_v_eb_e,
    const Matrix3d& old_C_b_e,
    const Vector3d& f_ib_b,
    const Vector3d& omega_ib_b,
    Vector3d& r_eb_e,
    Vector3d& v_eb_e,
    Matrix3d& C_b_e)
{

    double omega_ie = 7.292115E-5; // Earth rotation rate (rad/s)

    // ATTITUDE UPDATE
    double alpha_ie = omega_ie * tor_i;
    Matrix3d C_Earth;

    C_Earth << cos(alpha_ie), sin(alpha_ie), 0,
        -sin(alpha_ie), cos(alpha_ie), 0,
        0, 0, 1;
    //cout << "C_Earth " << endl;
    //cout << setprecision(15) << C_Earth << endl;

    Vector3d alpha_ib_b = omega_ib_b * tor_i;
    double mag_alpha = sqrt(alpha_ib_b.transpose() * alpha_ib_b);
    Matrix3d Alpha_ib_b = Skew_symmetric(alpha_ib_b);

    //cout << "alpha_ib_b " << endl;
    //cout << setprecision(15) << alpha_ib_b << endl;
    //cout << "mag_alpha " << endl;
    //cout << setprecision(15) << mag_alpha << endl;
    //cout << "Alpha_ib_b " << endl;
    //cout << setprecision(15) << Alpha_ib_b << endl;

    Matrix3d C_new_old;
    if (mag_alpha > 1.E-8)
    {
        C_new_old = Matrix3d::Identity() + sin(mag_alpha) / mag_alpha * Alpha_ib_b + (1.0 - cos(mag_alpha)) / (mag_alpha * mag_alpha) * Alpha_ib_b * Alpha_ib_b;
    }
    else
    {
        C_new_old = Matrix3d::Identity() + Alpha_ib_b;
    }

    // Update attitude
    C_b_e = C_Earth * old_C_b_e * C_new_old;
    //cout << "C_new_old "<< endl;
    //cout << setprecision(15) << C_new_old <<endl;

    // SPECIFIC FORCE FRAME TRANSFORMATION
    Matrix3d ave_C_b_e;
    if (mag_alpha > 1.E-8)
    {
        ave_C_b_e = old_C_b_e * (Matrix3d::Identity() + (1 - cos(mag_alpha)) / (mag_alpha * mag_alpha) * Alpha_ib_b + (1 - sin(mag_alpha) / mag_alpha) / (mag_alpha * mag_alpha) * Alpha_ib_b * Alpha_ib_b) - 0.5 * Skew_symmetric(Vector3d(0, 0, alpha_ie)) * old_C_b_e;
    }
    else
    {
        ave_C_b_e = old_C_b_e - 0.5 * Skew_symmetric(Vector3d(0, 0, alpha_ie)) * old_C_b_e;
    }

    // Transform specific force to ECEF-frame resolving axes
    Vector3d f_ib_e = ave_C_b_e * f_ib_b;
    //cout << "ave_C_b_e "<< endl;
    //cout<< ave_C_b_e <<endl;
    //cout << "f_ib_e "<< endl;
    //cout<< f_ib_e <<endl;

    // UPDATE VELOCITY
    v_eb_e = old_v_eb_e + tor_i * (f_ib_e + Gravity_ECEF(old_r_eb_e) - 2.0 * Skew_symmetric(Vector3d(0, 0, omega_ie)) * old_v_eb_e);

    // UPDATE CARTESIAN POSITION
    r_eb_e = old_r_eb_e + (v_eb_e + old_v_eb_e) * 0.5 * tor_i;
}

void NED_to_ECEF(Vector3d& r_eb_e, Vector3d& v_eb_e, Matrix3d& C_b_e, double L_b,
    double lambda_b, double h_b, Vector3d v_eb_n, Matrix3d C_b_n)
{
    // Parameters
    double R_0 = 6378137;       // WGS84 Equatorial radius in meters
    double e = 0.0818191908425; // WGS84 eccentricity

    // Calculate transverse radius of curvature using (2.105)
    double R_E = R_0 / sqrt(1.0 - (e * sin(L_b)) * (e * sin(L_b)));

    // Convert position using (2.112)
    double cos_lat = cos(L_b);
    double sin_lat = sin(L_b);
    double cos_long = cos(lambda_b);
    double sin_long = sin(lambda_b);

    r_eb_e << (R_E + h_b) * cos_lat * cos_long,
        (R_E + h_b)* cos_lat* sin_long,
        ((1.0 - e * e) * R_E + h_b)* sin_lat;

    // Calculate ECEF to NED coordinate transformation matrix using (2.150)
    Matrix3d C_e_n;
    C_e_n << -sin_lat * cos_long, -sin_lat * sin_long, cos_lat,
        -sin_long, cos_long, 0,
        -cos_lat * cos_long, -cos_lat * sin_long, -sin_lat;

    // Transform velocity using (2.73)
    v_eb_e = C_e_n.transpose() * v_eb_n;

    // Transform attitude using (2.15)
    C_b_e = C_e_n.transpose() * C_b_n;
}

void Radii_of_curvature(double L, double& R_N, double& R_E)
{

    const double R_0 = 6378137;       // WGS84 Equatorial radius in meters
    const double e = 0.0818191908425; // WGS84 eccentricity

    // std::cout << std::setprecision(15) << R_0 << std::endl;
    // std::cout << std::setprecision(15) << e << std::endl;

    // Calculate meridian radius of curvature using (2.105)
    double temp = 1 - (e * sin(L)) * (e * sin(L));
    R_N = R_0 * (1 - e * e) / pow(temp, 1.5);

    // Calculate transverse radius of curvature using (2.105)
    R_E = R_0 / sqrt(temp);
}

Matrix3d Skew_symmetric(Vector3d a)
{
    Matrix3d A;

    A << 0, -a(2), a(1),
        a(2), 0, -a(0),
        -a(1), a(0), 0;

    return A;
}

Matrix3d Initialize_NED_attitude(Matrix3d C_b_n, Vector3d delta_eul_nb_n)
{
    Matrix3d delta_C_b_n = Euler_to_CTM(-delta_eul_nb_n);
    Matrix3d est_C_b_n = delta_C_b_n * C_b_n;

    return est_C_b_n;
}