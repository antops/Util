/*
*	created by liyun 2017/10/23
*   function 图形的变换（平移和旋转）
*   version 1.0
*/
#ifndef GRAPHTRANS_H
#define GRAPHTRANS_H

#include <fstream>
#include <sstream>
#include "Vector3.h"
#include "Matrix4D.h"
#include "Vector3D.h"
#include "Constant_Var.h"
#include "Definition.h"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/eigenvalues"
using namespace std;
using namespace Eigen;
class GraphTrans
{
public:
	GraphTrans()
	{
		trans_x = 0;
		trans_y = 0;
		trans_z = 0;
		rotate_x = 1;
		rotate_y = 0;
		rotate_z = 0;
		rotate_theta = 0;
		U = Vector3(1, 0, 0);
		V = Vector3(0, 1, 0);
		N = Vector3(0, 0, 1);

	}
	void setGraphTransPar(const double& _trans_x, const double& _trans_y, const double& _trans_z,
		const double& _rotate_x, const double& _rotate_y, const double& _rotate_z, const double& _rotate_theta)
	{
		trans_x = _trans_x;
		trans_y = _trans_y;
		trans_z = _trans_z;
		rotate_x = _rotate_x;
		rotate_y = _rotate_y;
		rotate_z = _rotate_z;
		rotate_theta = _rotate_theta;
		Vector3D RotateAsix(rotate_x, rotate_y, rotate_z);
		Matrix4D matrix = Matrix4D::getRotateMatrix(rotate_theta, RotateAsix) *
			Matrix4D::getTranslateMatrix(trans_x, trans_y, trans_z);
		U = Matrix4D::getRotateMatrix(rotate_theta, RotateAsix) * Vector3(1, 0, 0);
		V = Matrix4D::getRotateMatrix(rotate_theta, RotateAsix) * Vector3(0, 1, 0);
		N = Matrix4D::getRotateMatrix(rotate_theta, RotateAsix) * Vector3(0, 0, 1);
	}
	void getGraphTransPar(double &_trans_x, double &_trans_y, double &_trans_z,
		double &_rotate_x, double &_rotate_y, double &_rotate_z, double &_rotate_theta) const
	{
		_trans_x = trans_x;
		_trans_y = trans_y;
		_trans_z = trans_z;
		_rotate_x = rotate_x;
		_rotate_y = rotate_y;
		_rotate_z = rotate_z;
		_rotate_theta = rotate_theta;
	}

	bool setUV(const Vector3& _U, const Vector3& _V)
	{ 
		
		//double temp = _U.Cross(_V)
		Vector3 tmpU = _U;
		Vector3 tmpV = _V;
		tmpU.Normalization();
		tmpV.Normalization();
		U = tmpU;
		V = tmpV;
		N = U.Cross(V);
		N.Normalization();
		Matrix3d A;
		A << U.x, V.x, N.x, U.y, V.y, N.y, U.z, V.z, N.z;
		EigenSolver<Matrix3d> es(A);
		Matrix3d DM = es.pseudoEigenvalueMatrix();
		Matrix3d VM = es.pseudoEigenvectors();
		Vector3 rotate_axis;
		for (int i = 0; i < 3; i++) {
			double eigenvalue = DM(i, i);
			if (eigenvalue > 0.999 && eigenvalue < 1.001) {
				rotate_axis.set(VM(0, i), VM(1, i), VM(2, i));
			}
		}
		rotate_axis.Normalization();
		//update Rotation Angle:
		double theta;
		Vector3 tempx(1, 0, 0);
		Vector3 tempy(0, 1, 0);
		Vector3 tempz(0, 0, 1);
		double rotate_theta;

		if (abs(tempz.Dot(rotate_axis)) < 0.9) {
			Vector3 v = tempz.Cross(rotate_axis);
			v.Normalization();
			Vector3 rv = Vector3(A(0,0)*v.x + A(0,1)*v.y + A(0,2)*v.z,
								 A(1,0)*v.x + A(1,1)*v.y + A(1,2)*v.z,
								 A(2,0)*v.x + A(2,1)*v.y + A(2,2)*v.z);
			rotate_theta = acos(v.Dot(rv) / v.Length() / rv.Length());
			rotate_theta = rotate_theta*180/Pi;
			//if (rotate_theta > 180) rotate_theta = rotate_theta - 360;
			//if (rotate_theta < -180) rotate_theta = rotate_theta + 360;
			updateRotate(rotate_axis, rotate_theta);
			//updateRotate(Vector3(1, 0, 0), 180);
			Matrix4D CheckM;
			CheckM = Matrix4D::getRotateMatrix(rotate_theta, rotate_axis.x,rotate_axis.y,rotate_axis.z);
			Vector3 Uc, Vc, Nc;
			Uc = CheckM*Vector3(1, 0, 0);	Vc = CheckM*Vector3(0, 1, 0);	Nc = CheckM*Vector3(0, 0, 1);
			if((Uc.Dot(U) > 0.999) && (Vc.Dot(V) > 0.999) && (Nc.Dot(N) > 0.999) )	return true;
			else {
				rotate_axis = rotate_axis*(-1);
				v = tempz.Cross(rotate_axis);
				v.Normalization();
				rv = Vector3(A(0, 0)*v.x + A(0, 1)*v.y + A(0, 2)*v.z,
					A(1, 0)*v.x + A(1, 1)*v.y + A(1, 2)*v.z,
					A(2, 0)*v.x + A(2, 1)*v.y + A(2, 2)*v.z);
				rotate_theta = acos(v.Dot(rv) / v.Length() / rv.Length());
				rotate_theta = rotate_theta * 180 / Pi;
				//if (rotate_theta > 180) rotate_theta = rotate_theta - 360;
				//if (rotate_theta < -180) rotate_theta = rotate_theta + 360;
				updateRotate(rotate_axis, rotate_theta);
				return true;
			}
		}
		else if(abs(tempx.Dot(rotate_axis)) < 0.9) {
			Vector3 v = tempx.Cross(rotate_axis);
			v.Normalization();
			Vector3 rv = Vector3(A(0, 0)*v.x + A(0, 1)*v.y + A(0, 2)*v.z,
								A(1, 0)*v.x + A(1, 1)*v.y + A(1, 2)*v.z,
				                A(2, 0)*v.x + A(2, 1)*v.y + A(2, 2)*v.z);
			rotate_theta = acos(v.Dot(rv) / v.Length() / rv.Length());
			rotate_theta = rotate_theta * 180 / Pi;
			//if (rotate_theta > 180) rotate_theta = rotate_theta - 360;
			//if (rotate_theta < -180) rotate_theta = rotate_theta + 360;
			updateRotate(rotate_axis, rotate_theta);
			//updateRotate(Vector3(1, 0, 0), 180);
			Matrix4D CheckM;
			CheckM = Matrix4D::getRotateMatrix(rotate_theta, rotate_axis.x, rotate_axis.y, rotate_axis.z);
			Vector3 Uc, Vc, Nc;
			Uc = CheckM*Vector3(1, 0, 0);	Vc = CheckM*Vector3(0, 1, 0);	Nc = CheckM*Vector3(0, 0, 1);
			if ((Uc.Dot(U) > 0.999) && (Vc.Dot(V) > 0.999) && (Nc.Dot(N) > 0.999))	return true;
			else {
				rotate_axis = rotate_axis*(-1);
				v = tempx.Cross(rotate_axis);
				v.Normalization();
				rv = Vector3(A(0, 0)*v.x + A(0, 1)*v.y + A(0, 2)*v.z,
					A(1, 0)*v.x + A(1, 1)*v.y + A(1, 2)*v.z,
					A(2, 0)*v.x + A(2, 1)*v.y + A(2, 2)*v.z);
				rotate_theta = acos(v.Dot(rv) / v.Length() / rv.Length());
				rotate_theta = rotate_theta * 180 / Pi;
				//if (rotate_theta > 180) rotate_theta = rotate_theta - 360;
				//if (rotate_theta < -180) rotate_theta = rotate_theta + 360;
				updateRotate(rotate_axis, rotate_theta);
				return true;
			}
		} 
		else if(abs(tempy.Dot(rotate_axis)) < 0.9) {
			Vector3 v = tempy.Cross(rotate_axis);
			v.Normalization();
			Vector3 rv = Vector3(A(0, 0)*v.x + A(0, 1)*v.y + A(0, 2)*v.z,
								A(1, 0)*v.x + A(1, 1)*v.y + A(1, 2)*v.z,
				                A(2, 0)*v.x + A(2, 1)*v.y + A(2, 2)*v.z);
			rotate_theta = acos(v.Dot(rv) / v.Length() / rv.Length());
			rotate_theta = rotate_theta * 180 / Pi;
			//if (rotate_theta > 180) rotate_theta = rotate_theta - 360;
			//if (rotate_theta < -180) rotate_theta = rotate_theta + 360;
			updateRotate(rotate_axis, rotate_theta);
			//updateRotate(Vector3(1, 0, 0), 180);
			Matrix4D CheckM;
			CheckM = Matrix4D::getRotateMatrix(rotate_theta, rotate_axis.x, rotate_axis.y, rotate_axis.z);
			Vector3 Uc, Vc, Nc;
			Uc = CheckM*Vector3(1, 0, 0);	Vc = CheckM*Vector3(0, 1, 0);	Nc = CheckM*Vector3(0, 0, 1);
			if ((Uc.Dot(U) > 0.999) && (Vc.Dot(V) > 0.999) && (Nc.Dot(N) > 0.999))	return true;
			else {
				rotate_axis = rotate_axis*(-1);
				v = tempy.Cross(rotate_axis);
				v.Normalization();
				rv = Vector3(A(0, 0)*v.x + A(0, 1)*v.y + A(0, 2)*v.z,
					A(1, 0)*v.x + A(1, 1)*v.y + A(1, 2)*v.z,
					A(2, 0)*v.x + A(2, 1)*v.y + A(2, 2)*v.z);
				rotate_theta = acos(v.Dot(rv) / v.Length() / rv.Length());
				rotate_theta = rotate_theta * 180 / Pi;
				//if (rotate_theta > 180) rotate_theta = rotate_theta - 360;
				//if (rotate_theta < -180) rotate_theta = rotate_theta + 360;
				updateRotate(rotate_axis, rotate_theta);
				return true;
			}
		}
		

		/*
		if (tmpU.Length() < THRESHOLD)
			return false;
		if (tmpV.Length() < THRESHOLD)
			return false;
		double tmp = tmpU.Dot(tmpV);
		if (abs(abs(tmp) - 1) < THRESHOLD)
			return false;

		U = _U;
		V = _V;
		V = V - U*(U.Dot(V));
		U.Normalization();
		V.Normalization();
		N = U.Cross(V);
		tmp = N.Dot(Vector3(0, 0, 1));

		if (abs(abs(tmp) -1)  > THRESHOLD)
		{
			Vector3 rotate_axis = Vector3(0, 0, 1).Cross(N); // 旋转轴
			double rotate_theta = acos(Vector3(0, 0, 1).Dot(N));
			rotate_theta = rotate_theta / Pi * 180;
			updateRotate(rotate_axis, rotate_theta);
		}
		else if (tmp > 0)
		{
			updateRotate(Vector3(0, 0, 1), 0);
		}
		else
		{
			updateRotate(Vector3(1, 0, 0), 180);
		}
		return true;
		*/
			
     	}

	void getUVN(Vector3& _U, Vector3& _V, Vector3& _N) const
	{
		_U = U;
		_V = V;
		_N = N;
	}

	void normalization(double fator) //变为标准单位 m
	{
		trans_x *= fator;
		trans_y *= fator;
		trans_z *= fator;
	}

	double getTrans_x() const { return trans_x; }
	double getTrans_y() const { return trans_y; }
	double getTrans_z() const { return trans_z; }

	double getRotate_x() const { return rotate_x; }
	double getRotate_y() const { return rotate_y; }
	double getRotate_z() const { return rotate_z; }

	double getRotate_theta() const { return rotate_theta; }

	void updateRotate(Vector3 rotate, double theta)
	{
		rotate_x = rotate.x;
		rotate_y = rotate.y;
		rotate_z = rotate.z;
		rotate_theta = theta;
		U = Matrix4D::getRotateMatrix(rotate_theta, rotate_x, rotate_y, rotate_z)*Vector3(1.0, 0.0, 0.0);
		V = Matrix4D::getRotateMatrix(rotate_theta, rotate_x, rotate_y, rotate_z)*Vector3(0.0, 1.0, 0.0);
		N = Matrix4D::getRotateMatrix(rotate_theta, rotate_x, rotate_y, rotate_z)*Vector3(0.0, 0.0, 1.0);
	}

	void updateTranslate(Vector3 tran)
	{
		trans_x = tran.x;
		trans_y = tran.y;
		trans_z = tran.z;
	}

	string getTransString() const
	{
		string ss;
		stringstream stream;
		stream << "(" << trans_x
			<< "," << trans_y
			<< "," << trans_z
			<< ")";
		stream >> ss;
		return ss;
	}

	string getRotateString() const
	{
		string ss;
		stringstream stream;
		if (rotate_theta == 0.0)
			return "0";

		stream << "(" << rotate_x
			<< "," << rotate_y
			<< "," << rotate_z
			<< ") " << rotate_theta;
		stream >> ss;
		return ss;
	}

	void save(ofstream & savefile) const
	{
		savefile << trans_x << " "
			<< trans_y << " "
			<< trans_z << " "
			<< rotate_x << " "
			<< rotate_y << " "
			<< rotate_z << " "
			<< rotate_theta << endl;
	}

	void open(ifstream & readfile)
	{
		readfile >> trans_x >> trans_y >> trans_z >> rotate_x
			>> rotate_y >> rotate_z >> rotate_theta;

	}

	static void updateSource_n(const Vector3& new_n, GraphTrans & _gt)
	{
		//new_n.Normalization();
		if (new_n.x != 0 || new_n.y != 0 || new_n.z != 1)
		{
			Vector3 rotate_axis = Vector3(0, 0, 1).Cross(new_n); // 旋转轴
			rotate_axis.Normalization();
			double rotate_theta = acos(Vector3(0, 0, 1).Dot(new_n));
			rotate_theta = rotate_theta / Pi * 180;
			_gt.updateRotate(rotate_axis, rotate_theta);
			_gt.U = Matrix4D::getRotateMatrix(rotate_theta, rotate_axis.x, rotate_axis.y, rotate_axis.z)*Vector3(1.0, 0.0, 0.0);
			_gt.V = Matrix4D::getRotateMatrix(rotate_theta, rotate_axis.x, rotate_axis.y, rotate_axis.z)*Vector3(0.0, 1.0, 0.0);
			_gt.N = Matrix4D::getRotateMatrix(rotate_theta, rotate_axis.x, rotate_axis.y, rotate_axis.z)*Vector3(0.0, 0.0, 1.0);
		}
		else
		{
			_gt.updateRotate(Vector3(0, 0, 1), 0);
		}
	}
private:
	double trans_x, trans_y, trans_z; // 平移的点
	double rotate_x, rotate_y, rotate_z; // 旋转轴
	double rotate_theta; // 旋转角度 右手定则

	Vector3 U, V, N;
};

#endif // !GRAPHTRANS_H

