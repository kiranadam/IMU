// run the file using  these instructions 
// $ g++ -std=c++11 $(pkg-config --cflags eigen3) imu.cpp -o imu 
// $ ./imu

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>

using namespace std;
using Eigen::MatrixXf;
using Eigen::Matrix3f;
using Eigen::VectorXf;
using Eigen::Vector3f;
using Eigen::Quaternionf;
using std::vector;
		
//function defination to get the matrix of raw data from CSV file
vector<Quaternionf> raw_data(string filename);

//function for the quaternion multiplication
Quaternionf quatMult(Quaternionf& q1, Quaternionf& q2);

//function with Kalman prediction
void KF_prediction(VectorXf& xkk_1, MatrixXf& Pkk_1, VectorXf u, VectorXf& xkk, MatrixXf& Pkk);

//function with Kalman correction
void KF_correction(VectorXf& xkk_1, MatrixXf& Pkk_1, Vector3f z, VectorXf& xkk, MatrixXf& Pkk);

//function to get sensor model;
Vector3f measurement(Quaternionf accel);

//function for quaternion to Euler angles
Vector3f quat2euler(Quaternionf q);

// check zero update system
bool ZVI(Quaternionf accel, Quaternionf gyro);

// initializer 
void initialize();


// some global variables 

float deltaT = 1.0/400.0;  // timestamp
float gyro_sigma = 3.141593 * (0.2/180.0); // gyro measurement error 2 degree / sec assumption
float accel_sigma = 0.02;   // accelrometer measurement error  0.02 m/sec^2 assumption
float g = 9.81; // accelration due to gravity

MatrixXf F; // dynamic matrix
MatrixXf G; // control matrix
MatrixXf Q; // plant covariance
MatrixXf H; // observation matrix
MatrixXf R; // measurement covariance

Vector3f z; // sensor measuremnt;
Vector3f e; // quaternion to euler angles
VectorXf u; // control gravity

VectorXf xkk; // corrected state
VectorXf xkk_1; // predicted state

MatrixXf Pkk; // corrected covariance
MatrixXf Pkk_1; // predicted covariance

Matrix3f gravity; // gravity matrix

Quaternionf q;  // initial Quaternion for earth frame 


int main(int , char**)
{
	// Read the csv files 
	vector<Quaternionf> gyro_list = raw_data("../../imu/set1/Gyro_Raw_2017-08-02-18-12-46-921.csv"); //gyro measurements
	vector<Quaternionf> acce_list = raw_data("../../imu/set1/Accelero_Raw_2017-08-02-18-12-46-922.csv"); // accelro measurements

	// initialize eath frame as below
	q.w() = 1.0;
	q.x() = 0.0;
	q.y() = 0.0;
	q.z() = 0.0;

        // calculate the estimate of angular velocity
	q = quatMult(q, gyro_list[2]);
	q.w() = 0.5*q.w();
	q.x() = 0.5*q.x();
	q.y() = 0.5*q.y();
	q.z() = 0.5*q.z();	

	// state and covariance initialization
	xkk = VectorXf(9);    // corrected state
	xkk_1 = VectorXf(9);  // predicted state
	
	xkk << 0, 0, 0, 0, 0, 0, quat2euler(q); // [x, y, z, vx, vy, vz, wx, wy, wz]
	xkk_1 = xkk;  
	
	// control vector to be given as an input of the system
	u = VectorXf(6); 
	
	// initialization with the values of the system to be checked with
	initialize();

	// writing the outputs in the file systems
	ofstream state;
	ofstream measure;
	ofstream pose;

	state.open("state.csv");
	measure.open("measure.csv");
	pose.open("pose.csv");

	for(unsigned i=0; i < acce_list.size(); i++)
	{
		// the gyroscope measuremt in the eath framw
		q = quatMult(q, gyro_list[i+2]);
		q.w() = 0.5*q.w();
		q.x() = 0.5*q.x();
		q.y() = 0.5*q.y();
		q.z() = 0.5*q.z();

		// the control values for the system as an input with accelrometer values and gyrovalues converted to euler angles
		u << acce_list[i].vec(), quat2euler(q);  // u = [ ax, ay, az, ex, ey ez]

		// Estimate the state with 
		KF_prediction(xkk_1, Pkk_1, u, xkk, Pkk);
		
		float dx,dy,dz;
		
		// get measurements 
		z = measurement(acce_list[i]);

		// writing measurements from sensor in the file measure.csv
		dx = z(0);
		dy = z(1);
		dz = z(2);
		measure << dx <<", "<<dy<<", "<<dx<<endl;

		//check zero velocity and angular velocity update
		if(ZVI(acce_list[i], gyro_list[i]) == true)
		{
			KF_correction(xkk_1, Pkk_1, z, xkk, Pkk);
			
			// writing state x,y,z in the file state.csv
			dx = xkk(0);
			dy = xkk(1);
			dz = xkk(2);
			state << dx <<", "<<dy<<", "<<dz<<endl;
			
			// writing state x,y,z in the file pose.csv
			dx = xkk(6)*deltaT;
			dy = xkk(7)*deltaT;
			dz = xkk(8)*deltaT;
			pose << dx <<", "<<dy<<", "<<dz<<endl;

		}
		else
		{
			Pkk << 8 * Matrix3f::Identity(), Matrix3f::Zero(), Matrix3f::Zero(),
	       		       Matrix3f::Zero(), 4 * Matrix3f::Identity(), Matrix3f::Zero(),
	                       Matrix3f::Zero(), Matrix3f::Zero(), 0.68 * Matrix3f::Identity();

			// writing state x,y,z in the file state.csv
			dx = xkk(0);
			dy = xkk(1);
			dz = xkk(2);
			state << dx <<", "<<dy<<", "<<dz<<endl;
			
			// writing state x,y,z in the file pose.csv
			dx = xkk(6)*deltaT;
			dy = xkk(7)*deltaT;
			dz = xkk(8)*deltaT;
			pose << dx <<", "<<dy<<", "<<dz<<endl;	
		}
	}
	state.close();
	measure.close();
	pose.close();
	
	return 0;
}


// function to read CSV file to convert into a raw data matrix
vector<Quaternionf> raw_data(string filename)
{
	// input file as stream
	ifstream infile(filename);
	
	if(!infile.is_open())
	{
		cout<<"file not found"<<endl;
		return {};
	}

	// Quaternion vector list
	vector<Quaternionf> vec;
	 
	// run till eof
	while(infile)
	{
		// token
    		string line;
    		// get line as a single token
		while(getline(infile, line))
		{
			stringstream ss(line);

			vector<double> t;
			// break the lines in multiple tokens with delimiter as ","
        		for(line; getline(ss, line, ','); )
			{
             			t.push_back(stof(line));	
        		}
			
			Quaternionf temp;
			temp.w() = 0;
			temp.x() = t[1];
			temp.y() = t[2];
			temp.z() = t[3];
			vec.push_back(temp);
		}	
	}

	return vec;
}

// Quaternion multiplication of two quaternions
Quaternionf quatMult(Quaternionf& q1, Quaternionf& q2) 
{
    Quaternionf result;
    result.setIdentity();

    result.w() = q1.w() * q2.w() - q1.vec().dot(q2.vec());
    result.vec() = q1.w() * q2.vec() + q2.w() * q1.vec() + q1.vec().cross(q2.vec());

    return result;
}

// Kalman Prediction 
void KF_prediction(VectorXf& xkk_1, MatrixXf& Pkk_1, VectorXf u, VectorXf& xkk, MatrixXf& Pkk)
{
	// gravity matrix for cross product
	gravity << 0, g, 0,
          	  -g, 0, 0,
                   0, 0, 0;

	// Dynamic matrix
	F = MatrixXf(9,9);
	
	F << Matrix3f::Identity(), Matrix3f::Identity()*deltaT, Matrix3f::Zero(),
             Matrix3f::Zero(), Matrix3f::Identity(), gravity*deltaT,
	     Matrix3f::Zero(), Matrix3f::Zero(), Matrix3f::Identity();

	// normalize q
	q.normalize();

	// Control matrix
	G = MatrixXf(9,6);
	
	G << Matrix3f::Zero(), Matrix3f::Zero(),
             q.toRotationMatrix()*deltaT, Matrix3f::Zero(),
	     Matrix3f::Zero(), q.toRotationMatrix();

	// kalman state and covariance prediction
	xkk_1 = F*xkk + G*u;
	Pkk_1 = F*Pkk*F.transpose() + G*Q*G.transpose();
}

// Kalman Correction
void KF_correction(VectorXf& xkk_1, MatrixXf& Pkk_1, Vector3f z, VectorXf& xkk, MatrixXf& Pkk)
{
	VectorXf nu = z - H * xkk_1;
	MatrixXf S = H * Pkk_1 * H.transpose() + R;
	MatrixXf W = Pkk_1 * H.transpose() * S.inverse();

	xkk = xkk_1 + W * nu;
	Pkk = (MatrixXf::Identity(Pkk_1.rows(),Pkk_1.cols()) - W * H ) * Pkk;
}

// measurement models
Vector3f measurement(Quaternionf accel)
{
	z(0) = z(0) + 0.5 * accel.x()*deltaT*deltaT;
	z(1) = z(1) + 0.5 * accel.y()*deltaT*deltaT;
	z(2) = z(2) + 0.5 * (accel.z()-g) * deltaT*deltaT;
	return z;
}

// quternions to euler angles
Vector3f quat2euler(Quaternionf q)
{
	Vector3f angle;
	
	angle(0) = atan(2.0f * (q.w()*q.x() + q.y()*q.z())/(1.0f - 2.0f*(q.x()*q.x() + q.y()*q.y())));
	angle(1) = asin(2.0f * (q.w()*q.y() - q.z()*q.x()));
	angle(2) = atan(2.0f * (q.w()*q.z() + q.x()*q.y())/(1.0f - 2.0f*(q.y()*q.y() + q.z()*q.z()))); 

	return angle;
}

// general initializer
void initialize()
{
	Pkk = MatrixXf(9,9);
	Pkk_1 = MatrixXf(9,9);

	Pkk << 8 * Matrix3f::Identity(), Matrix3f::Zero(), Matrix3f::Zero(),
	       Matrix3f::Zero(), 4 * Matrix3f::Identity(), Matrix3f::Zero(),
	       Matrix3f::Zero(), Matrix3f::Zero(), 0.68 * Matrix3f::Identity();

	Pkk_1 = Pkk;  
 
	// plant noise covaraince
	Q = MatrixXf(6,6);

	Q << 4 * Matrix3f::Identity(), Matrix3f::Zero(),
	     Matrix3f::Zero(), 0.68 * Matrix3f::Identity();

	// observational noise covariance

	R = MatrixXf(3,3);

	R << 8 * Matrix3f::Identity();

	H = MatrixXf(3,9);
	
	H << Matrix3f::Identity(), Matrix3f::Zero(), Matrix3f::Zero();

	// measurement initializer
	z << 0, 0, 0;
}

// function to check the zero update interval	
bool ZVI(Quaternionf accel, Quaternionf gyro)
{
	double a_norm = sqrt(accel.x()*accel.x() + accel.y()*accel.y() + accel.z()*accel.z());
	//cout<<a_norm<<endl;
	
	double w_norm = sqrt(gyro.x()*gyro.x() + gyro.y()*gyro.y() + gyro.z()*gyro.z());
	//cout<<w_norm<<endl;

	if((a_norm-g) <= sqrt(3)*accel_sigma && w_norm <= sqrt(3)*gyro_sigma)
		return false;
	else
		return true;
}
