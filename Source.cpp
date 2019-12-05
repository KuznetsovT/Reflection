//utf-8

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>



#include "Diffraction.h"       //Собраны воедино все связанные с дифракцией константы и функции
#include "R3.h"                //класс векторов в пространстве и 3х3 матриц.
#include "Gonio.h"             //класс, в котором реализовано вращение кристалла в гониометре.

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

	/*
	M = _M_._M_();    //находим матрицу ориентации в "прямом" пространстве (пространстве объекта)
	cout << endl << M << endl;
	cout << M.a.length() << " " << M.b.length() << " " << M.c.length() << endl;
	*/
	
	//sin th не может быть по модулю больше единицы!
	if (fabs(Diff::sin_th(_M_, h, k, l)) > 1) {
		cout << "|sin th| >1 diffraction doesn't valid!\n";
		return 0;
	}
	cout << "sin th = " << Diff::sin_th(_M_, h, k, l) << endl;

	Gonio g(_M_);
	cout << "{ ";
	for (auto d : g.diff_rotation(1, 1, h, k, l)) {
		cout << "[ " << rad_to_degrees(d.omega) << " " << rad_to_degrees(d.phi) << " " << rad_to_degrees(d.chi) << "]\n";
	}
	cout << "}\n";


	/*cout << "[ ";
	for (auto d : g.set_omega(g.diff_rotation(0,0,h, k, l)[0].omega).delta_omega_rotation(h, k, l)) {
		cout << rad_to_degrees(d.omega) << " ";
	}
	cout << "]\n";
	*/
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
