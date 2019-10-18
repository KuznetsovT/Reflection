﻿#include <iostream>
#include <cmath>
#include "R3.h"

#define PI 3.1415926
using namespace std;



//=================================================================================================================================================
//=================================================================================================================================================
//===============================БЛОК КОНСТАНТ=====================================================================================================

const double LYAMBDA = 1000;      //Длина вволны рентгеновского источника, в Ангстремах

//Матрица ориентации { a, b, c }, все координаты в Ангстремах
const matrix M = {    
	{1, 0, 0},
	{0, 1, 0},
	{0, 0, 1}
};

//'координаты' площади 'hkl'
const double h = (double) 1/2000, k = (double) 1/3000, l =(double) 1/4000;

//=================================================================================================================================================
//=================================================================================================================================================
//=================================================================================================================================================

//Находим матрицу ориентации в ОБРАТНОМ ПРОСТРАНСТВЕ
const matrix _M_ = {
	(M.b * M.c) / M.V(),
	(M.c * M.a) / M.V(),
	(M.a * M.b) / M.V()
};


//MAIN
int main() {
	R3 d = _M_.a * h + _M_.b * k + _M_.c * l;
	double s = LYAMBDA*d.length() ;
	cout << d << endl;
	cout << s << endl;
	cout << asin(s)*180/PI << " degrees" << endl;
 	return 0;
}
