/*
*	created by liyun 2018/7/8
*   function 
*   version 
*/
#ifndef FIELDBASE_H
#define FIELDBASE_H

#include <vector>
#include <fstream>
#include <complex>

#include "GraphTrans.h"


struct FieldBase
{
	GraphTrans graphTransField; // 旋转平移参数

	// 场分量
	vector<vector<complex<double>>> Ex, Ey, Ez;
	vector<vector<complex<double>>> Hx, Hy, Hz;

	int N_width, M_depth;  //N_width = para[0] /ds
	double ds_x, ds_y;
	
};

#endif // !DRAW_H


