#pragma once
#include <iostream>

//����� �������� {x, y, z}
class R3 {
public:
	double x = 0, y = 0, z = 0;

	R3 operator/(const double) const;
	R3 operator+(const R3&) const;
	double operator^(const R3&) const;    //��������� ������������  
	R3 operator*(const R3&) const;         //��������� ������������
	R3 operator*(const double) const;
	double length() const;

	friend std::istream & operator>> (std::istream &, R3 & );
	friend std::ostream & operator<< (std::ostream &, const R3 &);
};


//��������� ������������ ����� �� ��� �������� { a, b, c }
class matrix {
public:
	const R3 a, b, c;
	double V() const;
};



//���������� ������� ������� �� ������
R3 R3::operator/(const double k) const
{
	return { x / k, y / k, z / k };
}

//����������� ����� ���� ��������
R3 R3::operator+(const R3& r) const
{
	return { x + r.x, y + r.y, z + r.z };
}

//���������� ��������� ������������ ���� ��������.
double R3::operator^(const R3& r) const
{
	return x * r.x + y * r.y + z * r.z;
}

//���������� ��������� ������������ ���� ��������.
R3 R3::operator*(const R3& r) const
{
	return {
		y * r.z - z * r.y,
		z * r.x - x * r.z,
		x * r.y - y * r.x
	};
}

//���������� ���������� ������� �� ������
R3 R3::operator*(const double k) const
{
	return { x * k, y * k, z * k };
}

//����� �������.
double R3::length() const
{
	return sqrt(x * x + y * y + z * z);
}

//���������� ���� ������� ����� ����� �����
inline std::istream& operator>> (std::istream& in, R3& r)
{
	in >> r.x >> r.y >> r.z;
	return in;
}

inline std::ostream& operator<<(std::ostream& out, const R3& r)
{
	out << r.x << " " << r.y << " " << r.z;
	return out;
}


//�������, ���������� ��������������� ����� ��������������� �� ��� �������� ��������
double matrix::V() const
{
	return a ^ (b * c);
}
