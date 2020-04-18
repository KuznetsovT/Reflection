#include "Gonio.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include <vector>
#include <utility>
#include <string>

using namespace std;




Gonio::BOOM_predictor::Predictions Gonio::BOOM_predictor::is_available(const OPC & opc)
{
	Gonio::OPC deg_opc = rad_to_degrees(opc);

	//в файле данные находятся в виде строки. Мы находим два сhi, которые максимально хорошо подходят для данного chi
	//далее мы считываем соответствующие строки
	string up_line, down_line;

	//проверяем существование выданной нам директории
	if (!std::ifstream(file_constructor("chi"))) return Error;


	//мы находим две строки максимально похожие на исходное положение chi
	find_up_down(opc, file_constructor("chi"), up_line, down_line);
	

	//создаем векторы,  характеризующие возможные позиции omega 
	vector<pair<double, double>> up = create_position_data(up_line), down = create_position_data(down_line);

	

	//безопасное прохождение - это когда координаты поворота попадают и в верхний и в нижний строки
	//невозможное положение - координаты не попадают ни в одну из строк

	if (safety_prediction(deg_opc, up, down)) return Safe;
	if (impossible_prediction(deg_opc, up, down)) return Impossible;


	//далее будем производить линейную интерполяцию
	//продолжение следует
	return Unsafe;
}
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
std::string Gonio::BOOM_predictor::file_constructor(const char * geom)
{
	stringstream ss;
	ss << PATH << geom << "\\d_theta = " << d_theta << " d = " << d << " mm.txt";
	return ss.str();
}

void Gonio::BOOM_predictor::find_up_down(const Gonio::OPC & opc, string filename, string & up_line, string & down_line)
{
	Gonio::OPC deg_opc = rad_to_degrees(opc);
	double chi_up = 0xdafeca; double chi_down = -0xdafeca;
	up_line = ""; down_line = "";
	
	ifstream in(filename);
	
	//мы находим две строки максимально похожие на исходное положение chi
	while (!in.eof()) {

		string line;
		getline(in, line);
		stringstream ss(line);

		double chi;
		ss >> chi;

		if (chi <= deg_opc.chi && chi > chi_down) {
			chi_down = chi; down_line = line;
		}

		if (chi >= deg_opc.chi && chi < chi_up) {
			chi_up = chi; up_line = line;
		}
	}
	in.close();
}


//из строки данных создаём осмысленный вектор данных
//пара обозначает промежуток (first, second) - в котором разрешено находиться
//если first > second - промежуток проходит через 180 градусов, на котором меняется нумерация
std::vector<std::pair<double, double>> Gonio::BOOM_predictor::create_position_data(const std::string & line)
{
	vector<pair<double, double>> data;
	stringstream ss(line);
	double chi;
	ss >> chi;
	while (!ss.eof()) {
		pair<double, double> begin_end = { 0, 0 };

		ss >> begin_end.first >> begin_end.second;

		if (begin_end != pair<double, double>(0, 0)) data.push_back(begin_end);
	}
	return data;
}

bool Gonio::BOOM_predictor::safety_prediction(const OPC & opc, const std::vector<std::pair<double, double>>& up, const std::vector<std::pair<double, double>>& down)
{
	bool contains(const Gonio::OPC & opc, const std::vector<std::pair<double, double>> & data);

	return contains(opc, up) && contains(opc, down);
}

bool Gonio::BOOM_predictor::impossible_prediction(const OPC & opc, const std::vector<std::pair<double, double>>& up, const std::vector<std::pair<double, double>>& down)
{
	bool contains(const Gonio::OPC & opc, const std::vector<std::pair<double, double>> & data);

	return !contains(opc, up) && !contains(opc, down);
}

bool contains(const Gonio::OPC & opc, const std::vector<std::pair<double, double>> & data) {
	
	for (auto range : data) {
		
		if ((range.first < range.second) && (opc.omega > range.first && opc.omega < range.second)) return true;
		if ((range.first > range.second) && (opc.omega > range.first || opc.omega < range.second)) return true;
		
	}
	return false;
}