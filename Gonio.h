//utf-8
#pragma once

#include "R3.h"
#include "Diffraction.h"
#include <vector> //std::vector

////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////
//ОБЬЯВЛЕНИЯ

/* 
                   ↑Z
                   |                               
                   |
              ***  |  **
           **      |     **
          **              **
         **        @__ ___ **________________________\X
         **       //       **                        /
          **     //|      **          
           ***  // |    **                 
             ******|  *                    
 /*                |                   
                   |                
                   |                
                        
	omega - угол поворота кольца вокруг оси Z.
	  сhi - угол отклонения держателя кристала от вертикали.
	  phi - угол поворота кристалла в держателе относительно оси держателя.
*/


//класс, реализующий поворот кристалла в гониометре, использующем Эйлерову геометрию.
class Gonio {
public:

	//Omega, Phi, Chi
	struct OPC {
		double omega = 0;
		double phi = 0;
		double chi = 0;
	};

private:
	OPC opc;   //угловые параметры в Эйлеровой геометрии.
	const matrix _M_;    //Гониометр хранит матрицу Ориентации кристалла
	matrix _M_rotated;   //Матрица ориентации, после поворотов. По ней можно считать положения рефлексов.

	//Функция, возвращает угол в CФЕРИЧЕСКИХ координатах равный углу поворота вектора в плоскости ху в РАДИАНАХ
	static double Alpha(const R3& r);

	//Функция, возвращает угол в СФЕРИЧЕСКИХ координатах, равный углу между вектором и осью z в РАДИАНАХ
	static double Beta(const R3& r);


	//Возвращает матрицу поворот по оси Х
	static matrix X_rotation(double xang);
	//Возвращает матрицу поворота по оси У
	static matrix Y_rotation(double yang);
	//Возвращает матрицу поворота по оси Z
	static matrix Z_rotation(double zang);

	//Возвращает матрицу поворота в Эйлеровой геометрии
	static matrix OPC_rotation(OPC opc);

public:

	//Основной конструктор
	Gonio(const matrix& _M_);   

	//реализуем установки сразу omega, phi, chi.
	Gonio & set(const OPC& opc);
	Gonio & set_omega(const double omega);
	Gonio & set_phi(const double phi);
	Gonio & set_chi(const double chi);

	//функция возвращает повернутый базис, где X' - нормаль к оси кольца, Y' - параллельна оси кольца, Z' - ось держателя!
	matrix basic() const;
	matrix M_rotated() const;
	R3 S_rotated(const HKL &hkl) const;




	//--------------------------ОДНООСНЫЙ ОТРАЖАТЕЛЬ---------------------------------------


	//функция возвращает true, если возможен поворот по omega,  чтобы кристалл оказался в отражающем положении.
	bool is_omega_rotation_available(const HKL &hkl) const;


	//Обычно существует два решения. Поэтому оба решения мы возвращаем в виде пары. 
	//возвращает насколько нужно повернуть кристалл по omega, чтобы получить отражающее положение.
	std::vector<Gonio::OPC> delta_omega_rotation(const HKL &hkl) const;

	//возвращает пару omega, до которых нужно поверуть кристалл, чтобы быть в отражающем положении.
	std::vector<Gonio::OPC> omega_rotation(const HKL &hkl) const;







	//=========================ТРЁХОСНЫЙ ОТРАЖАТЕЛЬ=================================



	//решение задачи трёхосного отражателя
	std::vector<OPC> diff_rotation(const double psi, const double ksi, const HKL & hkl) const;


	//функция, возвращающая OPC для psi = 0 для минимального chi
	OPC psi0(const double ksi,const HKL & hkl) const; 

};



