#include <iostream>
#include <Eigen/Dense>
#include <iomanip> // 包含 std::setprecision
#include <cmath>
#include <fstream>
#include <random>
#include "FileIO.h"
#include <sstream>
#include <vector>


using namespace std;
using namespace Eigen;


/*void  fileToMatrix(const string& input_profile_name, int no_columns,
    MatrixXd& in_profile, int& no_epochs)
{

    ifstream file(input_profile_name);

    if (!file.is_open())
    {
        cout << "Error opening input files." << endl;
        return;
    }

    string line;
    no_epochs = 0;

    while (getline(file, line))
    {
        istringstream iss(line);
        VectorXd row(no_columns);
        //string value_str;

        for (int i = 0; i < no_columns; i++)
        {
            string value_str;
            iss >> value_str;
            //double value=stod(value_str);
            //row(i) = value;
            row(i) = stod(value_str);
        }

        if (no_epochs == 0)
        {
            // 调用 resize() 函数将矩阵的大小调整为1行no_columns列
            in_profile.resize(1, no_columns);
        }
        else
        {
            in_profile.conservativeResize(no_epochs + 1, no_columns);
        }

        in_profile.row(no_epochs) = row;
        no_epochs++;
    }

    file.close();
}*/
void  fileToMatrix(const string& input_profile_name, int no_columns,
    MatrixXd& in_profile, int& no_epochs)
{

    ifstream file(input_profile_name);

    if (!file.is_open())
    {
        cout << "Error opening input files." << endl;
        return;
    }

    string line;
    //no_epochs = 325017;
    int i = 0;
    no_epochs = 0;
    while (getline(file, line)) 
    {
        ++no_epochs;
    }

    file.clear();
    file.seekg(0, ios::beg);

    in_profile.resize(no_epochs, no_columns);

        while (getline(file, line))
        {
            
            istringstream iss(line);

            for (int j = 0; j < no_columns; ++j)
            {
                double value_str;
                iss >> value_str;
                in_profile(i, j) = value_str;

            }
            //cout << setprecision(15) << in_profile.row(i) << endl;
            ++i;
        }
    
    
    file.close();
}

bool Read_profile(const string& input_profile_name, const string& input_profile_name1, const string& input_profile_name2,
    MatrixXd& in_profile, int& no_epochs, MatrixXd& in_profile1, int& no_epochs1, MatrixXd& in_profile2, int& no_epochs2)
{

    no_epochs = 0;
    no_epochs1 = 0;
    no_epochs2 = 0;

    fileToMatrix(input_profile_name, 11, in_profile, no_epochs);
    fileToMatrix(input_profile_name1, 7, in_profile1, no_epochs1);
    fileToMatrix(input_profile_name2, 7, in_profile2, no_epochs2);

    // Check the number of columns in each file
    if (in_profile.cols() > 0 && in_profile.cols() != 11)
    {
        cout << input_profile_name << " has the wrong number of columns." << endl;
        return false;
    }

    if (in_profile1.cols() > 0 && in_profile1.cols() != 7)
    {
        cout << input_profile_name1 << " has the wrong number of columns." << endl;
        return false;
    }

    if (in_profile2.cols() > 0 && in_profile2.cols() != 7)
    {
        cout << input_profile_name2 << " has the wrong number of columns." << endl;
        return false;
    }

    // Convert degrees to radians
    in_profile.block(0, 2, no_epochs, 2) *= deg_to_rad;
    in_profile.block(0, 8, no_epochs, 3) *= deg_to_rad;

    in_profile2.block(0, 1, no_epochs2, 2) *= deg_to_rad;

    return true;
}