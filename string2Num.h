#pragma once
#include <iostream>
#include <sstream>
using namespace std;

template <class Type>
Type string2Num(const string& str) {
	istringstream iss(str);
	Type num;
	iss >> num;
	return num;
}
