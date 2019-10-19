#pragma once
#include <iostream>
#include <fstream>
#define PI 3.1415926535897931

//==================================================================================================================
//==================================================================================================================
//=======ОБЪЯВЛЕНИЯ=================================================================================================
//==================================================================================================================

//Класс векторов {x, y, z}
class R3 {
public:
	double x = 0, y = 0, z = 0;

	R3 operator/(const double) const;
	R3 operator+(const R3&) const;
	double operator^(const R3&) const;    //скалярное произведение  
	R3 operator*(const R3&) const;         //векторное произведение
	R3 operator*(const double) const;
	double length() const;

	friend std::istream& operator>> (std::istream&, R3&);
	friend std::ostream& operator<< (std::ostream&, const R3&);
};

//Функция, возвращающая угол в РАДИАНАХ между векторами a и b
double phi(const R3& a, const R3& b);

//перевод радиан в градусы
double rad_to_degrees(const double);

//структура определяющая базис из трёх векторов { a, b, c }
class matrix {
public:
	R3 a, b, c;
	double V() const;
	matrix _M_() const;
	friend std::istream& operator>>(std::istream&, matrix&);
	friend std::ostream& operator<<(std::ostream&, const matrix&);
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

