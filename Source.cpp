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




//--------------------@@@DEPRECATED@@@-------------------------------
//функция, считывающая данные из 'configuration.file' из файла по его полному или относительному пути.
void READ(const char* path_configuration_file, matrix& _M_, HKL & hkl, double & ksi, double & psi);



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
	ofstream out("output 1,0; 2,0; 3,0; psi = 180; .txt");        //файл, куда выводим траекторию
	out << "-h- -k- -l- -ksi- psi-  -2theta- [   -omega-   -phi-    -chi-    ]\n\n";
		reader.read( hkl, ksi, psi);
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) return 0;


		out << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;
		for (double ksi = 0; ksi <= 360; ksi += 1) {

			for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(ksi), hkl)) {
				cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
					<< "		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
				out << rad_to_degrees(d.omega) << "; " << rad_to_degrees(d.phi) << "; " << rad_to_degrees(d.chi) << " \n";
			}
		}
	
	out.close();

	//for (double i = 0; i < 2 * PI; i += 0.1) {
		//cout << "ksi = "<< i << " { ";
		//for (auto d : g.diff_rotation(i, 1, hkl)) {
		//	cout << "[ " << rad_to_degrees(d.omega) << " " << rad_to_degrees(d.phi) << " " << rad_to_degrees(d.chi) << "]\n";
		//}
		//cout << "}\n";
	//}

	//auto d = g.diff_rotation(0, 0, hkl)[0];
	//cout << "[ " << rad_to_degrees(d.omega) << " " << rad_to_degrees(d.phi) << " " << rad_to_degrees(d.chi) << "]\n";
	return 0;
}









//---------------@@@DEPRECATED@@@@-----------------------------

//функция считывает матрицу ориентации и hkl с из файла "PATH\\configuration.file"
void READ(const char* path_configuration_file, matrix& _M_, HKL & hkl, double & ksi, double & psi)
{
	ifstream in(path_configuration_file);
	if (!in.is_open()) {
		cerr << "CONFIG FILE NOT FOUND!\n";
	}
	{ string s; std::getline(in >> s, s); }   // "//матрица ориентации" считается ПОСТРОЧНО
	// матрица ориентации (UB-matrix) описывает базис обратной решетки
	in >> _M_;

	{ string s; std::getline(in >> s, s); }   // "s(hkl)"

	cout << "-h- -k- -l- -ksi- psi-  -2theta- [   -omega-   -phi-    -chi-    ]\n\n";
	Gonio g(_M_);
	while (!in.eof()) {

		in >> hkl.h >> hkl.k >> hkl.l >> ksi >> psi;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) return in.close();
		for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(ksi), hkl)) {
			cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << psi << "  " << 2*rad_to_degrees(Diff::th(_M_, hkl))
			<<	"		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
		}
	}
	in.close();
}

//--------------------@@@@@@@---------------------------
