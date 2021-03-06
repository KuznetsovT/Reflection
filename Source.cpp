﻿//utf-8

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cmath>



#include "Diffraction.h"       //Собраны воедино все связанные с дифракцией константы и функции
#include "R3.h"                //класс векторов в пространстве и 3х3 матриц.
#include "Gonio.h"             //класс, в котором реализовано вращение кристалла в гониометре.
#include "Readers.h"           //хранит два класса, UBreader матрицу ориентации, HKLreader считывает массивы hkl.

using namespace std;




//несколько вариантов обхода позиций

//выводит позиции с разными psi с известным шагом
void type0(const matrix & _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour);
//выводит позиции с разными ksi с некоторым шагом
void type1(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour);
//выводит все позиции с шагом psi и шагом ksi
void type2(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour);
//выводит позицию для данного рефлекса
void type_def(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour);

int main_program(int argn, char* argv[]);


int main(int argn, char*argv[]) {

	//не обращать на это внимания.
	//это всё debug и связанные вещи

	//смотреть main program!

/*
	Gonio::BOOM_predictor predictor(0, 35);

	cout << predictor.is_available({ 330,0,0 });
*/

	/*
	ifstream in("input.txt");
	ofstream out("output.txt");

	while (!in.eof()) {
		string line;
		getline(in, line);
		stringstream ss(line);
		double chi;
		ss >> chi;
		out << chi;
		while (!ss.eof()) {
			double num = 0;
			ss >> num;
			if (num != 0) out << "	" << rad_to_degrees(arc(cos(degrees_to_rades(num)), sin(degrees_to_rades(num))));
		}
		out << endl;
	}*/



	//for (;;) {
	//	Gonio::OPC opc, rad_opc;
	//	cin >> opc.chi >> opc.omega; opc.phi = 0;
	//	rad_opc = { degrees_to_rades(opc.omega), degrees_to_rades(opc.phi), degrees_to_rades(opc.chi) };

	//	for (auto opk : Gonio::EulerToKappa(rad_opc)) {
	//		cout << "kappa: " << rad_to_degrees(opk.kappa) << "  omega: " << rad_to_degrees(opk.omega) << "  phi: " << rad_to_degrees(opk.phi) << endl;
	//	}
	//	cout << endl;
	//}
	//return 0;

	return main_program(argn, argv);
}

//MAIN
int main_program(int argn, char* argv[]) { 

	
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

	


	cout << "-h- -k- -l- -ksi- psi-  -2theta- [   -omega-   -phi-    -chi-    ]\n\n";
	Gonio g(_M_);
	while (!reader.eof()) {
		ofstream endeavour;        //файл, куда выводим траекторию

		switch (reader.read(hkl, ksi, psi)) {
		case 0:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;
			
			{
				stringstream strstr;
				strstr << hkl.h << "_" << hkl.k << "_" << hkl.l << "_" << ksi << "_" << psi << "_0.xyz";
				endeavour.open(strstr.str());
			}

			type0(_M_, hkl, ksi, psi, endeavour);
			endeavour.close();
			break;
		case 1:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			{
				stringstream strstr;
				strstr << hkl.h << "_" << hkl.k << "_" << hkl.l << "_" << ksi << "_" << psi << "_1.xyz";
				endeavour.open(strstr.str());
			}

			type1(_M_, hkl, ksi, psi, endeavour);
			endeavour.close();
			break;

		case 2:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			{
				stringstream strstr;
				strstr << hkl.h << "_" << hkl.k << "_" << hkl.l << "_" << ksi << "_" << psi << "_2.xyz";
				endeavour.open(strstr.str());
			}

			type2(_M_, hkl, ksi, psi, endeavour);
			break;

		default:

			if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) break;

			{
				stringstream strstr;
				strstr << hkl.h << "_" << hkl.k << "_" << hkl.l << "_" << ksi << "_" << psi << "_def.xyz";
				endeavour.open(strstr.str());
			}

			type_def(_M_, hkl, ksi, psi, endeavour);
			
			break;

		}
		

	}

	return 0;
}






void type0(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;
	//сначала записываем в файл число точек
	endeavour << 2*(int)(360 / ksi)+8 << endl;
	//комментарий для файла
	endeavour << "ksi trajectory with " << ksi << " deg step" << " !! manualy cropped, in 10th degrees !!" << endl;
	endeavour << "U -18 -18 -10\nW 18 -18 -10\nF -18 18 -10\nH -18 -18 10\nB -18 18 10\nB 18 -18 10\nB 18 18 -10\nB 18 18 10\n";
	for (double k = 0; k <= 360; k += ksi) {
		char flag = 'C';
		for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(k), hkl, Gonio::Euler)) {
			if (flag == 'C') {
				endeavour << "C ";
				flag = 'O';
			}
			else {
				endeavour << "O ";
				flag = 'C';
			}
			cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << k << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
				<< "		[   " << rad_to_degrees(d.co.omega) << "  " << rad_to_degrees(d.co.phi) << "  " << rad_to_degrees(d.co.chi) << "	]	\n";
			endeavour << rad_to_degrees(d.co.omega)/10 << "  " << rad_to_degrees(d.co.phi)/10 << "  " << rad_to_degrees(d.co.chi)/10 << endl;

		}
		cout << endl;
	}
}

void type1(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << "-psi-" << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;
	//сначала записываем в файл число точек
	endeavour << 2 * (int)(360 / psi)+8 << endl;
	//комментарий для файла
	endeavour << "psi trajectory with " << psi << " deg step" << " !! manualy cropped, in 10th degrees !!" << endl;
	endeavour << "U -18 -18 -10\nW 18 -18 -10\nF -18 18 -10\nH -18 -18 10\nB -18 18 10\nB 18 -18 10\nB 18 18 -10\nB 18 18 10\n";
	for (double p = 0; p <= 360; p += psi) {
		char flag = 'C';
		for (auto d : g.diff_rotation(degrees_to_rades(p), degrees_to_rades(ksi), hkl, Gonio::Euler)) {
			if (flag == 'C') {
				endeavour << "C ";
				flag = 'O';
			}
			else {
				endeavour << "O ";
				flag = 'C';
			}
			cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << p << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
				<< "		[   " << rad_to_degrees(d.co.omega) << "  " << rad_to_degrees(d.co.phi) << "  " << rad_to_degrees(d.co.chi) << "	]	\n";

			endeavour << rad_to_degrees(d.co.omega)/10 << "  " << rad_to_degrees(d.co.phi)/10 << "  " << rad_to_degrees(d.co.chi)/10 << endl;

		}
		cout << endl;
	}
}

void type2(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << "-psi-" << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;

	//сначала записываем в файл число точек
	endeavour << 2 * (int)(360 / psi) * 2 * (int)(360/ksi) +8 << endl;
	//комментарий для файла
	endeavour << "ksi&psi trajectory with " << ksi<<" & "<< psi << " deg step" << " !! manualy cropped, in 10th degrees !!" << endl;
	endeavour << "U -18 -18 -10\nW 18 -18 -10\nF -18 18 -10\nH -18 -18 10\nB -18 18 10\nB 18 -18 10\nB 18 18 -10\nB 18 18 10\n";
	for (double k = 0; k <= 360; k += ksi) {
		for (double p = 0; p <= 360; p += psi) {
			char flag = 'C';
			for (auto d : g.diff_rotation(degrees_to_rades(p), degrees_to_rades(k), hkl, Gonio::Euler)) {

				if (flag == 'C') {
					endeavour << "C ";
					flag = 'O';
				}
				else {
					endeavour << "O ";
					flag = 'C';
				}

				cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << k << " " << p << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
					<< "		[   " << rad_to_degrees(d.co.omega) << "  " << rad_to_degrees(d.co.phi) << "  " << rad_to_degrees(d.co.chi) << "	]	\n";

				endeavour << rad_to_degrees(d.co.omega)/10 << "  " << rad_to_degrees(d.co.phi)/10 << "  " << rad_to_degrees(d.co.chi)/10 << endl;
			
			}
			cout << endl;
		}
	}
}

void type_def(const matrix& _M_, const HKL hkl, const double ksi, const double psi, ostream & endeavour)
{
	Gonio g(_M_);
	cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << "-ksi-" << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl)) << endl << endl;

	//сначала записываем в файл число точек
	endeavour << 2+8 << endl;
	//комментарий для файла
	endeavour << "position with " << ksi << ' & ' << psi << " !! manualy cropped, in 10th degrees !!" << endl;
	endeavour << "U -18 -18 -10\nW 18 -18 -10\nF -18 18 -10\nH -18 -18 10\nB -18 18 10\nB 18 -18 10\nB 18 18 -10\nB 18 18 10\n";
	char flag = 'C';
	for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(ksi), hkl, Gonio::Euler)) {

		if (flag == 'C') {
			endeavour << "C ";
			flag = 'O';
		}
		else {
			endeavour << "O ";
			flag = 'C';
		}

		cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
			<< "		[   " << rad_to_degrees(d.co.omega) << "  " << rad_to_degrees(d.co.phi) << "  " << rad_to_degrees(d.co.chi) << "	]	\n";

		endeavour << rad_to_degrees(d.co.omega)/10 << "  " << rad_to_degrees(d.co.phi)/10 << "  " << rad_to_degrees(d.co.chi)/10 << endl;
	}
	cout << endl;
}
