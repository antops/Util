#pragma once
/*
	Format Unified IO for FieldBase;
	Ming Jin 
*/
#include <vector>
#include <complex>
#include "../Util/GraphTrans.h"
#include "FieldBase.h"
#include <fstream>
#include <string>
#include <cmath>

using namespace std;

double Pii = atan(1.)*4.;

static int toFile(std::string _filename, FieldBase &_fieldbase, double frequency = 0) {
	//写文件
	double pi = atan(1.)*4.;
	double rad2deg = 180.0 / pi;
	ofstream outfile(_filename);
	stringstream temp;
	int flag= 0;
	if (_fieldbase.Hx.size() > 0) flag = 1;

	temp << _fieldbase.graphTransField.getTrans_x() << " "
		<< _fieldbase.graphTransField.getTrans_y() << " "
		<< _fieldbase.graphTransField.getTrans_z() << " "
		<< _fieldbase.graphTransField.getRotate_x() << " "
		<< _fieldbase.graphTransField.getRotate_y() << " "
		<< _fieldbase.graphTransField.getRotate_z() << " "
		<< _fieldbase.graphTransField.getRotate_theta() << " "
		<< _fieldbase.N_width << " " << _fieldbase.M_depth << " " << _fieldbase.ds_x << " "<<_fieldbase.ds_y <<" "
		<<flag <<" "<< frequency << endl;
	if (flag == 0) {
		for (int i = 0; i < _fieldbase.N_width; i++)
			for (int j = 0; j < _fieldbase.M_depth; j++)
			{
				temp
					<< abs(_fieldbase.Ex[i][j]) << " " << arg(_fieldbase.Ex[i][j]) * rad2deg << " "
					<< abs(_fieldbase.Ey[i][j]) << " " << arg(_fieldbase.Ey[i][j]) * rad2deg << " " << endl;
			}

	}
	else if (flag == 1) {
		for (int i = 0; i < _fieldbase.N_width; i++)
			for (int j = 0; j < _fieldbase.M_depth; j++)
			{
				temp
					<< abs(_fieldbase.Ex[i][j]) << " " << arg(_fieldbase.Ex[i][j]) * rad2deg << " "
					<< abs(_fieldbase.Ey[i][j]) << " " << arg(_fieldbase.Ey[i][j]) * rad2deg << " " 
					<< abs(_fieldbase.Ez[i][j]) << " " << arg(_fieldbase.Ez[i][j]) * rad2deg << " "
					<< abs(_fieldbase.Hx[i][j]) << " " << arg(_fieldbase.Hx[i][j]) * rad2deg << " " 
					<< abs(_fieldbase.Hy[i][j]) << " " << arg(_fieldbase.Hy[i][j]) * rad2deg << " "
					<< abs(_fieldbase.Hz[i][j]) << " " << arg(_fieldbase.Hz[i][j]) * rad2deg << " " << endl;
			}
	}
	outfile << temp.str();
	outfile.close();
	return 0;
}

static double fromFile(std::string _filename, FieldBase &_fieldbase) {
	double pi = atan(1.)*4.;
	double deg2rad = pi / 180;
	double frequency;
	int flag;
	//读文件
	fstream file;
	string filename = _filename;	file.open(_filename, ios::in);
	stringstream tempfile;	tempfile << file.rdbuf();
	string temp;		getline(tempfile, temp);
	istringstream tempLine(temp);
	double tx, ty, tz, rx, ry, rz, rth;	double du,dv;	
	int Nu, Nv;
	tempLine >> tx >> ty >> tz >> rx >> ry >> rz >> rth;
	tempLine >> Nu >> Nv >> du >> dv;
	tempLine >> flag >> frequency;

	_fieldbase.graphTransField.setGraphTransPar(tx, ty, tz, rx, ry, rz, rth);
	_fieldbase.ds_x = du;
	_fieldbase.ds_y = dv;
	_fieldbase.N_width = Nu;
	_fieldbase.M_depth = Nv;

	if (flag == 0) {
		_fieldbase.Ex.resize(_fieldbase.N_width);		_fieldbase.Ey.resize(_fieldbase.N_width);
		for (int i = 0; i < _fieldbase.N_width; i++) {
			_fieldbase.Ex[i].resize(_fieldbase.M_depth);	_fieldbase.Ey[i].resize(_fieldbase.M_depth);
		}
		//挨行读场分布
		double exa, exd, eya, eyd;
		double tempf;
		for (int i = 0; i < _fieldbase.N_width; i++)
			for (int j = 0; j < _fieldbase.M_depth; j++) {

				getline(tempfile, temp);
				istringstream perline(temp);
				//只读切向场
				perline >> exa >> exd >> eya >> eyd;
				_fieldbase.Ex[i][j] = exa*exp(complex<double>(0.0, exd*deg2rad));
				_fieldbase.Ey[i][j] = eya*exp(complex<double>(0.0, eyd*deg2rad));
			}

	}
	else if (flag == 1) {
		_fieldbase.Ex.resize(_fieldbase.N_width);		_fieldbase.Ey.resize(_fieldbase.N_width);	_fieldbase.Ez.resize(_fieldbase.N_width);
		_fieldbase.Hx.resize(_fieldbase.N_width);		_fieldbase.Hy.resize(_fieldbase.N_width);	_fieldbase.Hz.resize(_fieldbase.N_width);
		for (int i = 0; i < _fieldbase.N_width; i++) {
			_fieldbase.Ex[i].resize(_fieldbase.M_depth);	_fieldbase.Ey[i].resize(_fieldbase.M_depth);	_fieldbase.Ez[i].resize(_fieldbase.M_depth);
			_fieldbase.Hx[i].resize(_fieldbase.M_depth);	_fieldbase.Hy[i].resize(_fieldbase.M_depth);	_fieldbase.Hz[i].resize(_fieldbase.M_depth);
		}
		double a1, p1, a2, p2, a3, p3;		double a4, p4, a5, p5, a6, p6;
		double tempf;
		for (int i = 0; i < _fieldbase.N_width; i++)
			for (int j = 0; j < _fieldbase.M_depth; j++) {


				tempLine >> a1 >> p1 >> a2 >> p2 >> a3 >> p3 >> a4 >> p4 >> a5 >> p5 >> a6 >> p6;

				_fieldbase.Ex[i][j] = complex<double>(a1*cos(p1 * deg2rad), a1*sin(p1 * deg2rad));
				_fieldbase.Ey[i][j] = complex<double>(a2*cos(p2 * deg2rad), a2*sin(p2 * deg2rad));
				_fieldbase.Ez[i][j] = complex<double>(a3*cos(p3 * deg2rad), a3*sin(p3 * deg2rad));
				_fieldbase.Hx[i][j] = complex<double>(a4*cos(p4 * deg2rad), a4*sin(p4 * deg2rad));
				_fieldbase.Hy[i][j] = complex<double>(a5*cos(p5 * deg2rad), a5*sin(p5 * deg2rad));
				_fieldbase.Hz[i][j] = complex<double>(a6*cos(p6 * deg2rad), a6*sin(p6 * deg2rad));
			}
	}
	else {
		return -100;
	}


	file.close();
	return frequency;
}

static void fromFileExEy(std::string _filename, FieldBase &_fieldbase, int _cut) {
	//只取其中一小段
	double pi = atan(1.)*4.;
	FieldBase tempbase;
	//读文件
	fstream file;
	string filename = _filename;	file.open(_filename, ios::in);
	stringstream tempfile;		tempfile << file.rdbuf();
	string temp;
	getline(tempfile, temp);
	istringstream tempLine(temp);
	double tx, ty, tz, rx, ry, rz, rth;	double du,dv;
	int Nu, Nv;
	tempLine >> tx >> ty >> tz >> rx >> ry >> rz >> rth;
	tempLine >> Nu >> Nv >> du >> dv;

	tempbase.graphTransField.setGraphTransPar(tx, ty, tz, rx, ry, rz, rth);
	_fieldbase.graphTransField.setGraphTransPar(tx, ty, tz, rx, ry, rz, rth);
	tempbase.ds_x = du;	_fieldbase.ds_y = dv;
	tempbase.N_width = Nu;	_fieldbase.N_width = Nu - _cut * 2;
	tempbase.M_depth = Nv;	_fieldbase.M_depth = Nv - _cut * 2;

	//注意，大软件的面文件顺序是先Nx 再Ny, 而FDTD计算中是先Ny再Nx，在此进行转换
	tempbase.Ex.resize(tempbase.N_width);		tempbase.Ey.resize(tempbase.N_width);
	for (int i = 0; i < tempbase.N_width; i++) {
		tempbase.Ex[i].resize(tempbase.M_depth);	tempbase.Ey[i].resize(tempbase.M_depth);
	}
	_fieldbase.Ex.resize(_fieldbase.N_width);		_fieldbase.Ey.resize(_fieldbase.N_width);
	for (int i = 0; i < _fieldbase.N_width; i++) {
		_fieldbase.Ex[i].resize(_fieldbase.M_depth);	_fieldbase.Ey[i].resize(_fieldbase.M_depth);
	}
	//挨行读场分布
	float exa, exd, eya, eyd, hxa, hxd, hya, hyd;
	float tempf;
	for (int i = 0; i < tempbase.N_width; i++)
		for (int j = 0; j <tempbase.M_depth; j++) {

			getline(tempfile, temp);
			istringstream perline(temp);
			//只读切向场
			perline >> exa >> exd >> eya >> eyd >> tempf >> tempf >> hxa >> hxd >> hya >> hyd >> tempf >> tempf;

			tempbase.Ex[i][j] = exa*exp(complex<float>(0.0, exd*pi / (180.0)));
			tempbase.Ey[i][j] = eya*exp(complex<float>(0.0, eyd*pi / (180.0)));
		}
	file.close();

	for (int i = 0; i < _fieldbase.N_width; i++)
		for (int j = 0; j < _fieldbase.M_depth; j++)
		{
			_fieldbase.Ex[i][j] = tempbase.Ex[i + _cut][j + _cut];
			_fieldbase.Ey[i][j] = tempbase.Ey[i + _cut][j + _cut];
		}
}

static int toFileBinary(std::string _filename, FieldBase &_fieldbase, double frequency = 0) {
	//写文件
	FILE* writefile;
	errno_t err;
	err = fopen_s(&writefile, _filename.c_str(), "wb");
	if (err != 0) return err;
	double tx, ty, tz, rx, ry, rz, rth;
	double dx, dy;
	int flag;
	tx = _fieldbase.graphTransField.getTrans_x();
	ty = _fieldbase.graphTransField.getTrans_y();
	tz = _fieldbase.graphTransField.getTrans_z();
	rx = _fieldbase.graphTransField.getRotate_x();
	ry = _fieldbase.graphTransField.getRotate_y();
	rz = _fieldbase.graphTransField.getRotate_z();
	rth = _fieldbase.graphTransField.getRotate_theta();

	int Nx = _fieldbase.N_width;
	int Ny = _fieldbase.M_depth;
	if (Nx < 1 && _fieldbase.Ex.size() < 1) return -1;//nothing to writh, or empty;
	else if (Nx < 1 && _fieldbase.Ex.size() > 0) {	//有数据但没有更新数组大小
		_fieldbase.N_width = _fieldbase.Ex.size();
		_fieldbase.M_depth = _fieldbase.Ex[0].size();
	}
	else if (Nx > 0 && _fieldbase.Ex.size() < 1) {	//有数据大小但没有数据
		_fieldbase.Ex.resize(Nx);	_fieldbase.Ey.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			_fieldbase.Ex[i].resize(Ny);	_fieldbase.Ey[i].resize(Ny);
		}
	}

	dx = _fieldbase.ds_x;	dy = _fieldbase.ds_y;

	fwrite(&frequency, sizeof(double), 1, writefile);

	fwrite(&tx, sizeof(double), 1, writefile);
	fwrite(&ty, sizeof(double), 1, writefile);
	fwrite(&tz, sizeof(double), 1, writefile);

	fwrite(&rx, sizeof(double), 1, writefile);
	fwrite(&ry, sizeof(double), 1, writefile);
	fwrite(&rz, sizeof(double), 1, writefile);

	fwrite(&rth, sizeof(double), 1, writefile);

	fwrite(&Nx, sizeof(int), 1, writefile);
	fwrite(&Ny, sizeof(int), 1, writefile);

	fwrite(&dx, sizeof(double), 1, writefile);
	fwrite(&dy, sizeof(double), 1, writefile);

	if (_fieldbase.Hx.empty()) flag = 0;
	else flag = 1;

	fwrite(&flag, sizeof(int), 1, writefile);

	if (flag == 0) {
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Ex[i][j], sizeof(complex<double>), 1, writefile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Ey[i][j], sizeof(complex<double>), 1, writefile);
			}
	}
	else {
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Ex[i][j], sizeof(complex<double>), 1, writefile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Ey[i][j], sizeof(complex<double>), 1, writefile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Ez[i][j], sizeof(complex<double>), 1, writefile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Hx[i][j], sizeof(complex<double>), 1, writefile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Hy[i][j], sizeof(complex<double>), 1, writefile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Ny; j++) {
				fwrite(&_fieldbase.Hz[i][j], sizeof(complex<double>), 1, writefile);
			}
	}
	fclose(writefile);
	return err;
}

static double fromFileBinary(std::string _filename, FieldBase &_fieldbase)
{

	FILE* readfile = NULL;
	errno_t err;
	err = fopen_s(&readfile, _filename.c_str(), "rb");
	if (err != 0) return -100;
	double tx, ty, tz, rx, ry, rz, rth;
	double freq;
	int Nx, Ny;
	double dx = 0.0;	double dy = 0.0;
	fread(&freq, sizeof(double), 1, readfile);

	fread(&tx, sizeof(double), 1, readfile);
	fread(&ty, sizeof(double), 1, readfile);
	fread(&tz, sizeof(double), 1, readfile);

	fread(&rx, sizeof(double), 1, readfile);
	fread(&ry, sizeof(double), 1, readfile);
	fread(&rz, sizeof(double), 1, readfile);

	fread(&rth, sizeof(double), 1, readfile);

	fread(&Nx, sizeof(int), 1, readfile);
	fread(&Ny, sizeof(int), 1, readfile);

	fread(&dx, sizeof(double), 1, readfile);
	fread(&dy, sizeof(double), 1, readfile);

	_fieldbase.graphTransField.setGraphTransPar(tx, ty, tz, rx, ry, rz, rth);
	_fieldbase.ds_x = dx;	_fieldbase.ds_y = dy;
	_fieldbase.N_width = Nx;	_fieldbase.M_depth = Ny;

	int flag;
	fread(&flag, sizeof(int), 1, readfile);
	if (flag == 0) {
		_fieldbase.Ex.resize(Nx);		_fieldbase.Ey.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			_fieldbase.Ex[i].resize(Ny);		_fieldbase.Ey[i].resize(Ny);
		}
	}
	else {
		_fieldbase.Ex.resize(Nx);		_fieldbase.Ey.resize(Nx);		_fieldbase.Ez.resize(Nx);
		_fieldbase.Hx.resize(Nx);		_fieldbase.Hy.resize(Nx);		_fieldbase.Hz.resize(Nx);
		for (int i = 0; i < Nx; i++) {
			_fieldbase.Ex[i].resize(Ny);		_fieldbase.Ey[i].resize(Ny);		_fieldbase.Ez[i].resize(Ny);
			_fieldbase.Hx[i].resize(Ny);		_fieldbase.Hy[i].resize(Ny);		_fieldbase.Hz[i].resize(Ny);
		}
	}

	//下面是读取环节
	if (flag == 0) {
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j<Ny; j++) {
				fread(&_fieldbase.Ex[i][j], sizeof(complex<double>), 1, readfile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j<Ny; j++) {
				fread(&_fieldbase.Ey[i][j], sizeof(complex<double>), 1, readfile);
			}
	}
	else {
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++) {
				fread(&_fieldbase.Ex[i][j], sizeof(complex<double>), 1, readfile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++) {
				fread(&_fieldbase.Ey[i][j], sizeof(complex<double>), 1, readfile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++) {
				fread(&_fieldbase.Ez[i][j], sizeof(complex<double>), 1, readfile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++) {
				fread(&_fieldbase.Hx[i][j], sizeof(complex<double>), 1, readfile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++) {
				fread(&_fieldbase.Hy[i][j], sizeof(complex<double>), 1, readfile);
			}
		for (int i = 0; i < Nx; i++)
			for (int j = 0; j < Nx; j++) {
				fread(&_fieldbase.Hz[i][j], sizeof(complex<double>), 1, readfile);
			}
	}
	fclose(readfile);
	return  freq;
}

static FieldBase cutField(FieldBase _fieldbase, int _cut) {
	FieldBase result;
	result.graphTransField = _fieldbase.graphTransField;
	result.ds_x = _fieldbase.ds_x;
	result.ds_y = _fieldbase.ds_y;
	result.M_depth = _fieldbase.M_depth - _cut * 2;
	result.N_width = _fieldbase.N_width - _cut * 2;
	if (_fieldbase.Hx.size() > 0) {
		result.Ex.resize(result.N_width);		result.Ey.resize(result.N_width);	result.Ez.resize(result.N_width);
		result.Hx.resize(result.N_width);		result.Hy.resize(result.N_width);	result.Hz.resize(result.N_width);
		for (int i = 0; i < result.N_width; i++) {
			result.Ex[i].resize(result.M_depth);	result.Ey[i].resize(result.M_depth);	result.Ez[i].resize(result.M_depth);
			result.Hx[i].resize(result.M_depth);	result.Hy[i].resize(result.M_depth);	result.Hz[i].resize(result.M_depth);
		}
		for (int i = 0; i < result.N_width; i++)
			for (int j = 0; j < result.M_depth; j++)
			{
				result.Ex[i][j] = _fieldbase.Ex[i + _cut][j + _cut];
				result.Ey[i][j] = _fieldbase.Ey[i + _cut][j + _cut];
				result.Ez[i][j] = _fieldbase.Ez[i + _cut][j + _cut];
				result.Hx[i][j] = _fieldbase.Hx[i + _cut][j + _cut];
				result.Hy[i][j] = _fieldbase.Hy[i + _cut][j + _cut];
				result.Hz[i][j] = _fieldbase.Hz[i + _cut][j + _cut];
			}
	}
	else {
		result.Ex.resize(result.N_width);		result.Ey.resize(result.N_width);
		for (int i = 0; i < result.N_width; i++) {
			result.Ex[i].resize(result.M_depth);	result.Ey[i].resize(result.M_depth);
		}
		for (int i = 0; i < result.N_width; i++)
			for (int j = 0; j < result.M_depth; j++)
			{
				result.Ex[i][j] = _fieldbase.Ex[i + _cut][j + _cut];
				result.Ey[i][j] = _fieldbase.Ey[i + _cut][j + _cut];
			}
	}

	return result;
}

static FieldBase addField(FieldBase _fieldbase, int _add) {
	//场的两边补0.。。
	FieldBase result;
	result.graphTransField = _fieldbase.graphTransField;
	result.ds_x = _fieldbase.ds_x;
	result.ds_y = _fieldbase.ds_y;
	result.M_depth = _fieldbase.M_depth + _add * 2;
	result.N_width = _fieldbase.N_width + _add * 2;
	complex<double> zero(0.0,0.0);
	if (_fieldbase.Hx.size() > 0) {
		result.Ex.resize(result.N_width);		result.Ey.resize(result.N_width);	result.Ez.resize(result.N_width);
		result.Hx.resize(result.N_width);		result.Hy.resize(result.N_width);	result.Hz.resize(result.N_width);
		for (int i = 0; i < result.N_width; i++) {
			result.Ex[i].resize(result.M_depth);	result.Ey[i].resize(result.M_depth);	result.Ez[i].resize(result.M_depth);
			result.Hx[i].resize(result.M_depth);	result.Hy[i].resize(result.M_depth);	result.Hz[i].resize(result.M_depth);
		}
		for (int i = 0; i < result.N_width; i++)
			for (int j = 0; j < result.M_depth; j++)
			{
				result.Ex[i][j] = zero;
				result.Ey[i][j] = zero;
				result.Ez[i][j] = zero;
				result.Hx[i][j] = zero;
				result.Hy[i][j] = zero;
				result.Hz[i][j] = zero;
			}
		for (int i = 0; i < _fieldbase.N_width; i++)
			for (int j = 0; j <_fieldbase.M_depth; j++)
			{
				result.Ex[i + _add][j + _add] = _fieldbase.Ex[i][j];
				result.Ey[i + _add][j + _add] = _fieldbase.Ey[i][j];
				result.Ez[i + _add][j + _add] = _fieldbase.Ez[i][j];
				result.Hx[i + _add][j + _add] = _fieldbase.Hx[i][j];
				result.Hy[i + _add][j + _add] = _fieldbase.Hy[i][j];
				result.Hz[i + _add][j + _add] = _fieldbase.Hz[i][j];
			}
	}
	else {
		result.Ex.resize(result.N_width);		result.Ey.resize(result.N_width);	
		for (int i = 0; i < result.N_width; i++) {
			result.Ex[i].resize(result.M_depth);	result.Ey[i].resize(result.M_depth);	
		}
		for (int i = 0; i < result.N_width; i++)
			for (int j = 0; j < result.M_depth; j++)
			{
				result.Ex[i][j] = zero;
				result.Ey[i][j] = zero;
			}
		for (int i = 0; i < _fieldbase.N_width; i++)
			for (int j = 0; j <_fieldbase.M_depth; j++)
			{
				result.Ex[i + _add][j + _add] = _fieldbase.Ex[i][j];
				result.Ey[i + _add][j + _add] = _fieldbase.Ey[i][j];

			}
	}

	return result;
}