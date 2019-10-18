﻿#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include "R3.h"


using namespace std;

void READ(const char*);

//=================================================================================================================================================
//=================================================================================================================================================
//===============================БЛОК КОНСТАНТ=====================================================================================================

const double LYAMBDA = 0.001;      //Длина вволны рентгеновского источника, в Ангстремах

//Матрица ориентации { a, b, c }, все координаты в Ангстремах
matrix M = {    
	{0, 0, 0},
	{0, 0, 0},
	{0, 0, 0}
};

//'координаты' площади 'hkl'
double h = (double) 0, k = (double) 0, l =(double) 0;

//=================================================================================================================================================
//=================================================================================================================================================
//=================================================================================================================================================

//Находим матрицу ориентации в ОБРАТНОМ ПРОСТРАНСТВЕ
matrix _M_ = {
	{0, 0, 0},
	{0, 0, 0},
	{0, 0, 0}
};

//MAIN
int main() {
	READ("C:\\Users\\yater\\source\\repos\\Reflection\\configuration.file");

	_M_ = {
	   (M.b * M.c) / M.V(),
	   (M.c * M.a) / M.V(),
	   (M.a * M.b) / M.V()
	};

	R3 d = _M_.a * h + _M_.b * k + _M_.c * l;
	double sin_theta = LYAMBDA * d.length();
	double alpha = (d.y >= 0) ?
		phi({ d.x, d.y,0 }, { 1, 0, 0 }) :
		-phi({ d.x, d.y,0 }, { 1, 0, 0 });
	double beta = asin((d.operator^({ 0,0,1 }) / d.length()));
	cout << d << endl;
	cout << rad_to_degrees(alpha) << " "
		<< rad_to_degrees(beta) << endl;

	cout << rad_to_degrees(asin(sin_theta)) << " degrees" << endl;
 	return 0;
}

void READ(const char* file_name)
{
	ifstream in(file_name);
	{ string s; std::getline(in >> s,s); }
	in >> M;

	{ string s; std::getline(in >> s,s); }
	in >> h >> k >> l;
	in.close();
}
