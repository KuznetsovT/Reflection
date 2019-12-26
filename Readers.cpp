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
	return _M_;
}

UBreader::~UBreader()
{
	in.close();
}


void HKLreader::read(HKL& hkl)
{
	if (in.is_open() && !in.eof()) {
		in >> hkl.h >> hkl.k >> hkl.l;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) in.close();
	}
}

bool HKLreader::open(const char* hkl_file)
{
	in.open(hkl_file);
	{ string s; std::getline(in >> s, s); }   // "s(hkl)" 
	return in.is_open();
}

void HKLreader::read(HKL& hkl, double& ksi)
{
	if (in.is_open() && !in.eof()) {
		in >> hkl.h >> hkl.k >> hkl.l >> ksi;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) in.close();
	}
}

void HKLreader::read(HKL& hkl, double& ksi, double& psi)
{
	if (in.is_open() && !in.eof()) {
		in >> hkl.h >> hkl.k >> hkl.l >> ksi >> psi;
		if (hkl.h == 0 && hkl.k == 0 && hkl.l == 0) in.close();
	}
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
