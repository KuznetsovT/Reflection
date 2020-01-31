//utf-8

#include <iostream>
#include <fstream>
#include <string>
#include <cmath>



#include "Diffraction.h"       //Собраны воедино все связанные с дифракцией константы и функции
#include "R3.h"                //класс векторов в пространстве и 3х3 матриц.
#include "Gonio.h"             //класс, в котором реализовано вращение кристалла в гониометре.
#include "Readers.h"           //хранит два класса, UBreader матрицу ориентации, HKLreader считывает массивы hkl.

using namespace std;




void type0(const matrix & _M_, const HKL hkl, const double ksi, const double psi);
void type1(const matrix& _M_, const HKL hkl, const double ksi, const double psi);
void type2(const matrix& _M_, const HKL hkl, const double ksi, const double psi);
void type_def(const matrix& _M_, const HKL hkl, const double ksi, const double psi);


//MAIN
int main(int argn, char* argv[]) { 


	matrix M;      //Матрица ориентации { a, b, c }, все координаты в Ангстремах
	matrix _M_;    //Определяем матрицу ориентации{  a* b* c* } в ОБРАТНОМ ПРОСТРАНСТВЕ 
	HKL hkl = { 0, 0, 0 };    //hkl 
	double ksi = 0, psi = 0;

	HKLreader reader;

	switch (argn) {
	case 3:
		_M_ = UBreader(argv[1]).readUB();
		reader.open(argv[2]);
		break;
	case 2:
		_M_ = UBreader(argv[1]).readUB();
		reader.open("hkl.file");
		break;
	default:
		_M_ = UBreader("UB.file").readUB();
		reader.open("hkl.file");
	}
	//нет дополнительных аргументов - считываем матрицу ориентации и hkl в текущей папкe
	//считываем матрицу ориентации и hkl из файл, указанного вторым параметром 

	
	//sin th не может быть по модулю больше единицы!
	//if (fabs(Diff::sin_th(_M_, hkl)) > 1) {
	//	cout << "|sin th| >1 diffraction doesn't valid!\n";
	//	return 0;
	//}
	//cout << "sin th = " << Diff::sin_th(_M_, hkl) << endl;


	/*

	Считываем hkl далее выводим для всех psi или ksi с шагом 1 градус все получаемые opc.

	*/

	cout << "-h- -k- -l- -ksi- psi-  -2theta- [   -omega-   -phi-    -chi-    ]\n\n";
	Gonio g(_M_);
	//ofstream out("output 1,0; 2,0; 3,0; psi = 180; .txt");        //файл, куда выводим траекторию
	while (!reader.eof()) {
		switch (reader.read(hkl, ksi, psi)) {
		case 0:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			type0(_M_, hkl, ksi, psi);
			break;
		case 1:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			type1(_M_, hkl, ksi, psi);
			break;

		case 2:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			type2(_M_, hkl, ksi, psi);
			break;

		default:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			type_def(_M_, hkl, ksi, psi);
			
			break;

		}
		

	}

	return 0;
}






void type0(const matrix& _M_, const HKL hkl, const double ksi, const double psi)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;

	for (double k = 0; k <= 360; k += ksi) {

		for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(k), hkl)) {
			cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << k << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
				<< "		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
		}
		cout << endl;
	}
}

void type1(const matrix& _M_, const HKL hkl, const double ksi, const double psi)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << "-psi-" << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;

	for (double p = 0; p <= 360; p += psi) {

		for (auto d : g.diff_rotation(degrees_to_rades(p), degrees_to_rades(ksi), hkl)) {
			cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << p << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
				<< "		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
		}
		cout << endl;
	}
}

void type2(const matrix& _M_, const HKL hkl, const double ksi, const double psi)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << "-psi-" << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;

	for (double k = 0; k <= 360; k += ksi) {
		for (double p = 0; p <= 360; p += psi) {

			for (auto d : g.diff_rotation(degrees_to_rades(p), degrees_to_rades(k), hkl)) {
				cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << k << " " << p << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
					<< "		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
			}
			cout << endl;
		}
	}
}

void type_def(const matrix& _M_, const HKL hkl, const double ksi, const double psi)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;

	for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(ksi), hkl)) {
		cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
			<< "		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
	}
	cout << endl;
}
