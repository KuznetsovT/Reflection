//utf-8

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <clocale>


/*
	класс векторов в трехмероном пространстве и базисов в пространстве
	рекомендую предварительно посмотреть файл "R3.h".
*/

#include "Diffraction.h"
#include "R3.h"  
#include "Gonio.h"

using namespace std;

//функция, считывающая данные из 'configuration.file' из файла по его полному или относительному пути.
void READ(const char* path_configuration_file, matrix& _M_, double& h, double& k, double& l);






//MAIN
int main(int argn, char* argv[]) { 

	matrix M;      //Матрица ориентации { a, b, c }, все координаты в Ангстремах
	matrix _M_;    //Определяем матрицу ориентации{  a* b* c* } в ОБРАТНОМ ПРОСТРАНСТВЕ 
	double h = 0, k = 0, l = 0;    //hkl 


	if (argn == 1) READ("configuration.file", _M_, h, k, l);   //нет дополнительных аргументов - считываем матрицу ориентации и hkl из configuration.file в текущей папкe
	else READ(argv[1], _M_, h, k, l);                          //считываем матрицу ориентации и hkl из файл, указанного вторым параметром 


	M = _M_._M_();    //находим матрицу ориентации в "прямом" пространстве (пространстве объекта)
	cout << endl << M << endl;
	cout << M.a.length() << " " << M.b.length() << " " << M.c.length() << endl;

	
	//sin th не может быть по модулю больше единицы!
	if (abs(Diff::sin_th(_M_, h, k, l)) > 1) {
		cout << "|sin th| >1 diffraction doesn't valid!\n";
		return 0;
	}
	cout << "sin th = " << Diff::sin_th(_M_, h, k, l) << endl;

	std::pair<double, double> delta = Gonio(_M_).omega_rotation_angle(h, k, l);
	cout << " [ " << rad_to_degrees(delta.first) << " || " << rad_to_degrees(delta.second) << " ] \n";
	return 0;
}


//функция считывает матрицу ориентации и hkl с из файла "PATH\\configuration.file"
void READ(const char* path_configuration_file, matrix& _M_, double& h, double& k, double& l)
{
	ifstream in(path_configuration_file);
	if (!in.is_open()) {
		cerr << "CONFIG FILE NOT FOUND!\n";
	}
	{ string s; std::getline(in >> s, s); }   // "//матрица ориентации"
	// матрица ориентации (UB-matrix) описывает базис обратной решетки
	in >> _M_;

	{ string s; std::getline(in >> s, s); }   // "s(hkl)"
	in >> h >> k >> l;
	in.close();
}
