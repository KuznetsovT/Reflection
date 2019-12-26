//utf-8
#include "R3.h"

//файл для определений из файла R3.h


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

//Вводим оператор умножения вектора на матрицу(поворота). !! Поворот вектора !!
R3 R3::operator*(const matrix& rotation) const
{
	return {
		*this ^ rotation[0],
		*this ^ rotation[1],
		*this ^ rotation[2] };
}

//Длина вектора!
double R3::length() const
{
	return sqrt(x * x + y * y + z * z);
}

//оператор вывода el_num - элемента вектора
double R3::operator[](const unsigned el_num) const {
	switch (el_num) {
	case 0: return x;
	case 1: return y;
	case 2: return z;
	default: return 0.0;  // catch exception ???
	}
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
double ang(const R3& a, const R3& b)
{
	return acos((a ^ b) / (a.length() * b.length()));
}

//Перевод рад в градусы
double rad_to_degrees(const double rad)
{
	return (rad * 180 / PI);
}

double degrees_to_rades(const double deg)
{
	return (deg*PI)/180;
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
matrix matrix::_M_() const
{
	return {
	   (b * c) / V(),
	   (c * a) / V(),
	   (a * b) / V()
	};
}

//Оператор, возвращает СТОЛБЕЦ матрицы
R3 matrix::operator[](const unsigned col_num) const {
	return { a[col_num], b[col_num], c[col_num] };
}

//Вводим умножение двух матриц, будем использовать в случае поворота кристалла
matrix matrix::operator*(const matrix& m) const
{
	return { a * m, b * m, c * m };
}

