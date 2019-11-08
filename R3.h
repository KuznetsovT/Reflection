//utf-8
#pragma once
#include <iostream>
#include <fstream>
#define PI 3.1415926535897932384626433832795    //Википедия, ctrl+C -> ctrl+V

//==================================================================================================================
//==================================================================================================================
//==========================ОБЪЯВЛЕНИЯ==============================================================================
//==================================================================================================================

//Класс векторов {x, y, z}
class R3 {
public:
	double x = 0, y = 0, z = 0;

	//производим перегрузку операторов для удобной арифметики векторов

	R3 operator/(const double) const;     //деление вектора на число
	R3 operator+(const R3&) const;        //сложение двух векторов
	double operator^(const R3&) const;    //скалярное произведение  
	R3 operator*(const R3&) const;        //векторное произведение
	R3 operator*(const double) const;     //домножение на число
	double length() const;                //модуль вектора

	friend std::istream& operator>> (std::istream&, R3&);                //перегружаем оператор ввода 
	friend std::ostream& operator<< (std::ostream&, const R3&);          //перегружаем оператор вывода
	double operator[] (const unsigned el_num) const {          //оператор вывода el_num - элемента вектора
		switch (el_num) {
		case 0: return x;
		case 1: return y;
		case 2: return z;
		default: return 0.0;  // catch exception ???
		}
	}
};

/////////////////////////////ТРИГОНОМЕТРИЯ\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

//Функция, возвращающая угол в РАДИАНАХ между векторами a и b
double phi(const R3& a, const R3& b);

//Функция, возвращает угол в сферических координатах равный углу поворота вектора в плоскости ху в РАДИАНАХ
double Alpha(const R3& r) {
	return (r.y >= 0) ?  // в направлении tth
		phi({ r.x, r.y, 0 }, { 1, 0, 0 }) :
		-phi({ r.x, r.y, 0 }, { 1, 0, 0 });
}

//Функция, возвращает угол в сферических координатах, равный углу между вектором и осью z в РАДИАНАХ
double Beta(const R3& r) {
	return acos(r.z / r.length());
}

//перевод радиан в градусы
double rad_to_degrees(const double);

///////////////////////////МАТРИЦЫ\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

//структура определяющая базис из трёх векторов-строк { a, b, c }
class matrix {
public:
	R3 a, b, c;
	double V() const;
	matrix _M_() const;
	friend std::istream& operator>>(std::istream&, matrix&);
	friend std::ostream& operator<<(std::ostream&, const matrix&);
	R3 operator[](const unsigned col_num) const {
		return { a[col_num], b[col_num], c[col_num] };
	}

	matrix operator*(const matrix& m) const;
};


/*===================================================================================================================
=====================================================================================================================
==============================ОПРЕДЕЛЕНИЯ============================================================================
=====================================================================================================================
*/

//Определяем деление вектора на скаляр
R3 R3::operator/(const double k) const
{
	return { x / k, y / k, z / k };
}

//Опредделяем сумму двух векторов
R3 R3::operator+(const R3& r) const
{
	return { x + r.x, y + r.y, z + r.z };
}

//Определяем скалярное произведение двух векторов.
double R3::operator^(const R3& r) const
{
	return x * r.x + y * r.y + z * r.z;
}

//Определяем векторное произведение двух векторов.
R3 R3::operator*(const R3& r) const
{
	return {
		y * r.z - z * r.y,
		z * r.x - x * r.z,
		x * r.y - y * r.x
	};
}

//Определяем домножение вектора на скаляр
R3 R3::operator*(const double k) const
{
	return { x * k, y * k, z * k };
}

//Длина вектора!
double R3::length() const
{
	return sqrt(x * x + y * y + z * z);
}

//Определяем ввод вектора через поток ввода
std::istream& operator>> (std::istream& in, R3& r)
{
	in >> r.x >> r.y >> r.z;
	return in;
}

//Определяем вывод вектора через поток вывода
std::ostream& operator<<(std::ostream& out, const R3& r)
{
	out << r.x << " " << r.y << " " << r.z;
	return out;
}

//Возвращает угол между двумя векторами в РАДИАНАХ
double phi(const R3& a, const R3& b)
{
	return acos((a ^ b) / (a.length() * b.length()));
}

//Перевод рад в градусы
double rad_to_degrees(const double rad)
{
	return (rad * 180 / PI);
}

//Оператор ввода матрицы 3х3 
std::istream& operator>>(std::istream& in, matrix& m)
{
	in >> m.a >> m.b >> m.c;
	return in;
}

//Оператор вывода матрицы 3х3
std::ostream& operator<<(std::ostream& out, const matrix& m)
{
	out << m.a << std::endl
		<< m.b << std::endl
		<< m.c << std::endl;
	return out;
}


//Функция, возвращает ориентированный обьём параллелограмма на трёх базисных векторах
double matrix::V() const
{
	return a ^ (b * c);
}

//Возвращает матрицу, обратную данной
inline matrix matrix::_M_() const
{
	return {
	   (b * c) / V(),
	   (c * a) / V(),
	   (a * b) / V()
	};
}

//Вводим умножение двух матриц, будем использовать в случае поворота кристалла
 matrix matrix::operator*(const matrix& m) const
{
	 return {
		 {a ^ m[0], a ^ m[1], a ^ m[2] },
		 {b ^ m[0], b ^ m[1], b ^ m[2] },
		 {c ^ m[0], c ^ m[1], c ^ m[2] }
	 };
}

