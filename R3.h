#pragma once
#include <iostream>
#include <fstream>
#define PI 3.1415926535897931

// ласс векторов {x, y, z}
class R3 {
public:
	double x = 0, y = 0, z = 0;

	R3 operator/(const double) const;
	R3 operator+(const R3&) const;
	double operator^(const R3&) const;    //скал€рное произведение  
	R3 operator*(const R3&) const;         //векторное произведение
	R3 operator*(const double) const;
	double length() const;

	friend std::istream & operator>> (std::istream &, R3 & );
	friend std::ostream & operator<< (std::ostream &, const R3 &);
};

//‘ункци€, возвращающа€ угол в –јƒ»јЌј’ между векторами a и b
double phi(const R3& a, const R3& b);

//перевод радиан в градусы
double rad_to_degrees(const double);

//структура определ€юща€ базис из трЄх векторов { a, b, c }
class matrix {
public:
	R3 a, b, c;
	double V() const;
	friend std::istream& operator>>(std::istream&, matrix&);
	friend std::ostream& operator<<(std::ostream&, const matrix&);
};



//ќпредел€ем деление вектора на скал€р
R3 R3::operator/(const double k) const
{
	return { x / k, y / k, z / k };
}

//ќпреддел€ем сумму двух векторов
R3 R3::operator+(const R3& r) const
{
	return { x + r.x, y + r.y, z + r.z };
}

//ќпредел€ем скал€рное произведение двух векторов.
double R3::operator^(const R3& r) const
{
	return x * r.x + y * r.y + z * r.z;
}

//ќпредел€ем векторное произведение двух векторов.
R3 R3::operator*(const R3& r) const
{
	return {
		y * r.z - z * r.y,
		z * r.x - x * r.z,
		x * r.y - y * r.x
	};
}

//ќпредел€ем домножение вектора на скал€р
R3 R3::operator*(const double k) const
{
	return { x * k, y * k, z * k };
}

//ƒлина вектора.
double R3::length() const
{
	return sqrt(x * x + y * y + z * z);
}

//ќпредел€ем ввод вектора через поток ввода
std::istream& operator>> (std::istream& in, R3& r)
{
	in >> r.x >> r.y >> r.z;
	return in;
}

std::ostream& operator<<(std::ostream& out, const R3& r)
{
	out << r.x << " " << r.y << " " << r.z;
	return out;
}

double phi(const R3& a, const R3& b)
{
	return acos((a ^ b) / (a.length() * b.length()));
}

double rad_to_degrees(const double rad)
{
	return (rad*180/PI);
}

std::istream& operator>>(std::istream& in, matrix& m)
{
	in >> m.a >> m.b >> m.c;
	return in;
}

std::ostream& operator<<(std::ostream& out,const matrix& m)
{
	out << m.a << std::endl
		<< m.b << std::endl
		<< m.c << std::endl;
	return out;
}


//‘ункци€, возвращает ориентированный обьЄм параллелограмма на трЄх базисных векторах
double matrix::V() const
{
	return a ^ (b * c);
}
