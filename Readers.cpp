#include "Readers.h"

#include <string>
using std::string;
using std::cerr;


UBreader::UBreader(const char* ub_file)
{
	in.open(ub_file);
	if (!in.is_open()) {
		cerr << "UB FILE NOT FOUND!\n";
	}
}

matrix UBreader::readUB()
{
	{ string s; std::getline(in >> s, s); }   // "//матрица ориентации" считается ПОСТРОЧНО
	// матрица ориентации (UB-matrix) описывает базис обратной решетки
	matrix _M_;
	in >> _M_;
	return _M_.T();
}

UBreader::~UBreader()
{
	in.close();
}


int HKLreader::read(HKL& hkl)
{
	if (in.is_open() && !in.eof()) {
		in >> hkl.h >> hkl.k >> hkl.l;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) in.close();

		int type = -1;
		in >> type;
		return type;
	}
	return -1;
}

bool HKLreader::open(const char* hkl_file)
{
	in.open(hkl_file);
	{ string s; std::getline(in >> s, s); }   // "s(hkl)" 
	return in.is_open();
}

int HKLreader::read(HKL& hkl, double& ksi)
{
	if (in.is_open() && !in.eof()) {
		in >> hkl.h >> hkl.k >> hkl.l >> ksi;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) in.close();

		int type = -1;
		in >> type;
		return type;
	}
	return -1;
}

int HKLreader::read(HKL& hkl, double& ksi, double& psi)
{
	if (in.is_open() && !in.eof()) {
		in >> hkl.h >> hkl.k >> hkl.l >> ksi >> psi;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) in.close();

		int type = -1;
		in >> type;
		return type;
	}
	return -1;
}

HKLreader::~HKLreader()
{
	if (in.is_open())
		in.close();
}

bool HKLreader::eof() const
{
	return in.is_open() && in.eof();
}

#include "Gonio.h"

//---------------@@@DEPRECATED@@@@-----------------------------

//функция считывает матрицу ориентации и hkl с из файла "PATH\\configuration.file"
void READ(const char* path_configuration_file, matrix& _M_, HKL& hkl, double& ksi, double& psi)
{
	using namespace std;
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
		for (auto d : g.diff_rotation(degrees_to_rades(psi), degrees_to_rades(ksi), hkl, Gonio::Euler)) {
			cout << hkl.h << " " << hkl.k << " " << hkl.l << " " << ksi << " " << psi << "  " << 2 * rad_to_degrees(Diff::th(_M_, hkl))
				<< "		[   " << rad_to_degrees(d.omega) << "  " << rad_to_degrees(d.phi) << "  " << rad_to_degrees(d.chi) << "	]	\n";
		}
	}
	in.close();
}

//--------------------@@@@@@@---------------------------

