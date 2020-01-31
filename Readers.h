#pragma once
#include "R3.h"
#include "Diffraction.h"
#include <fstream>
using std::ifstream;
using std::ofstream;



//--------------------@@@DEPRECATED@@@-------------------------------
//функция, считывающая данные из 'configuration.file' из файла по его полному или относительному пути.
void READ(const char* path_configuration_file, matrix& _M_, HKL& hkl, double& ksi, double& psi);



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
	int read(HKL& hkl);
	int read(HKL& hkl, double& ksi);
	int read(HKL& hkl, double& ksi, double& psi);
	~HKLreader();

	bool eof() const;
};
