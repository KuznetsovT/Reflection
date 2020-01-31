//utf-8
#pragma once
#include <iostream>
#include <fstream>
#include <cmath>
#define PI 3.1415926535897932384626433832795    //Википедия, ctrl+C -> ctrl+V

//==================================================================================================================
//==================================================================================================================
//==========================ОБЪЯВЛЕНИЯ==============================================================================
//==================================================================================================================
class matrix;   //Обьявление класса матрицы для определения умножения вектора на матрицу поворота. см ниже↓

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
	R3 operator*(const matrix &) const;     //умножение вектора на матрицу(поворота)
	double length() const;                //модуль вектора
	double operator[] (const unsigned el_num) const; //вывод элемента по номеру, x=0, y=1, z=2;


	friend std::istream& operator>> (std::istream&, R3&);                //перегружаем оператор ввода 
	friend std::ostream& operator<< (std::ostream&, const R3&);          //перегружаем оператор вывода
	
};

/////////////////////////////ТРИГОНОМЕТРИЯ\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

//Функция, возвращающая угол в РАДИАНАХ между векторами a и b
double ang(const R3& a, const R3& b);

//перевод радиан в градусы
double rad_to_degrees(const double);

//перевод градусов в радианы
double degrees_to_rades(const double);

///////////////////////////МАТРИЦЫ\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

//структура определяющая базис из трёх векторов-строк { a, b, c }
class matrix {
public:
	R3 a = { 0,0,0 }, b = { 0,0,0 }, c = { 0,0,0 };
	double V() const;
	matrix _M_() const; //обрятная матрица
	matrix T() const;  //транспонирование
	friend std::istream& operator>>(std::istream&, matrix&);
	friend std::ostream& operator<<(std::ostream&, const matrix&);
	R3 operator[](const unsigned col_num) const;

	matrix operator*(const matrix& m) const;
};

