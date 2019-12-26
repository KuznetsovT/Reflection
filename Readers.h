#pragma once
#include "R3.h"
#include "Diffraction.h"
#include <fstream>
using std::ifstream;
using std::ofstream;

//считывает матрицу ориентации из файла, по умолчанию "UB.file"
class UBreader {
	ifstream in;
public:
	UBreader(const char* ub_file = "UB.file");
	matrix readUB();
	~UBreader();
};



//cчитывает массивы hkl из файла, по умолчанию "HKL.file"
class HKLreader {
	ifstream in;
public:
	bool open(const char* hkl_file = "HKL.file");
	void read(HKL& hkl);
	void read(HKL& hkl, double& ksi);
	void read(HKL& hkl, double& ksi, double& psi);
	~HKLreader();

	bool eof() const;
};
