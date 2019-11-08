//utf-8
#pragma once

#include "R3.h"
#include "Diffraction.h"
#include <utility> //std::pair

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//ОБЬЯВЛЕНИЯ

/* 
                   ↑Z
                   |                               
                   |  
                ** | **
             **    |    *
            **           *
           **      @_____*________________________\X
           **     //     *                        /
           **    //     *          
            ** ////|  **                 
                ***|                    
 /*                |                   
                   |                
                   |                
                        
	omega - угол поворота кольца вокруг оси Z.
	  сhi - угол отклонения держателя кристала от вертикали.
	  phi - угол поворота кристалла в держателе относительно оси держателя.
*/


//класс, реализующий поворот кристалла в гониометре, использующем Эйлерову геометрию.
class Gonio {
private:
	double omega = 0.0, phi = 0.0, chi = 0.0;  //угловые параметры в Эйлеровой геометрии.
	const matrix _M_;    //Гониометр хранит матрицу Ориентации кристалла
	matrix _M_rotated;   //Матрица ориентации, после поворотов. По ней можно считать положения рефлексов.

	//Функция, возвращает угол в CФЕРИЧЕСКИХ координатах равный углу поворота вектора в плоскости ху в РАДИАНАХ
	static double Alpha(const R3& r);

	//Функция, возвращает угол в СФЕРИЧЕСКИХ координатах, равный углу между вектором и осью z в РАДИАНАХ
	static double Beta(const R3& r);

public:

	//Основной конструктор
	Gonio(const matrix& _M_);

	//реализуем установки сразу omega, phi, chi.
	Gonio& set(const double omega, const double phi, const double chi);
	Gonio& set_omega(const double omega);
	Gonio& set_phi(const double phi);
	Gonio& set_chi(const double chi);

	//функция возвращает повернутый базис, где X' - нормаль к оси кольца, Y' - параллельна оси кольца, Z' - ось держателя!
	matrix basic() const;


	/*Задача трёхосного отражателя.
	Если закрепить phi и chi, то есть считать что они нам заданы -> получаем задачу одноосного отражателя с вращением 'omega'.
	В нашем решении был пункт, что если угол между осью z и вектором s(hkl) рефлекса не входит в промежуток (-th; th),
	Отражающего положения добиться невозможно. Далее нужно для всех phi подобрать такие chi, чтобы кристалл был в отражающем положении.
	или выяснить такие phi, для которых это невозможно.
	*/

	//Далее запишем решение задачи одноосного отражателя. в виде функции rotate_omega(hkl)
	//функция возвращает true, если возможен поворот по omega,  чтобы кристалл оказался в отражающем положении.
	bool is_omega_rotation_available(hkl) const;

	//Запишем решение задачи одноосного отражателя. в виде функции rotate_omega(hkl).
	//Обычно существует два решения. Поэтому оба решения мы возвращаем в виде пары. 
	std::pair<double, double> rotate_omega(hkl) const;
};


