#pragma once
#include <iostream>
#include <fstream>
#define PI 3.1415926535897931

//Êëàññ âåêòîðîâ {x, y, z}
class R3 {
public:
	double x = 0, y = 0, z = 0;

	R3 operator/(const double) const;
	R3 operator+(const R3&) const;
	double operator^(const R3&) const;    //ñêàëÿðíîå ïðîèçâåäåíèå  
	R3 operator*(const R3&) const;         //âåêòîðíîå ïðîèçâåäåíèå
	R3 operator*(const double) const;
	double length() const;

	friend std::istream & operator>> (std::istream &, R3 & );
	friend std::ostream & operator<< (std::ostream &, const R3 &);
};

//Ôóíêöèÿ, âîçâðàùàþùàÿ óãîë â ÐÀÄÈÀÍÀÕ ìåæäó âåêòîðàìè a è b
double phi(const R3& a, const R3& b);

//ïåðåâîä ðàäèàí â ãðàäóñû
double rad_to_degrees(const double);

//ñòðóêòóðà îïðåäåëÿþùàÿ áàçèñ èç òð¸õ âåêòîðîâ { a, b, c }
class matrix {
public:
	R3 a, b, c;
	double V() const;
	friend std::istream& operator>>(std::istream&, matrix&);
	friend std::ostream& operator<<(std::ostream&, const matrix&);
};



//Îïðåäåëÿåì äåëåíèå âåêòîðà íà ñêàëÿð
R3 R3::operator/(const double k) const
{
	return { x / k, y / k, z / k };
}

//Îïðåääåëÿåì ñóììó äâóõ âåêòîðîâ
R3 R3::operator+(const R3& r) const
{
	return { x + r.x, y + r.y, z + r.z };
}

//Îïðåäåëÿåì ñêàëÿðíîå ïðîèçâåäåíèå äâóõ âåêòîðîâ.
double R3::operator^(const R3& r) const
{
	return x * r.x + y * r.y + z * r.z;
}

//Îïðåäåëÿåì âåêòîðíîå ïðîèçâåäåíèå äâóõ âåêòîðîâ.
R3 R3::operator*(const R3& r) const
{
	return {
		y * r.z - z * r.y,
		z * r.x - x * r.z,
		x * r.y - y * r.x
	};
}

//Îïðåäåëÿåì äîìíîæåíèå âåêòîðà íà ñêàëÿð
R3 R3::operator*(const double k) const
{
	return { x * k, y * k, z * k };
}

//Äëèíà âåêòîðà.
double R3::length() const
{
	return sqrt(x * x + y * y + z * z);
}

//Îïðåäåëÿåì ââîä âåêòîðà ÷åðåç ïîòîê ââîäà
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


//Ôóíêöèÿ, âîçâðàùàåò îðèåíòèðîâàííûé îáü¸ì ïàðàëëåëîãðàììà íà òð¸õ áàçèñíûõ âåêòîðàõ
double matrix::V() const
{
	return a ^ (b * c);
}

