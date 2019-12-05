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
	R3 S_rotated(hkl) const;

	//Далее запишем решение задачи одноосного отражателя. в виде функции rotate_omega(hkl)
	//функция возвращает true, если возможен поворот по omega,  чтобы кристалл оказался в отражающем положении.
	bool is_omega_rotation_available(hkl) const;

	//Запишем решение задачи одноосного отражателя. в виде функции delta_omega_rotation(hkl).
	//Обычно существует два решения. Поэтому оба решения мы возвращаем в виде пары. 
	//возвращает насколько нужно повернуть кристалл по omega, чтобы получить отражающее положение.
	std::vector<Gonio::OPC> delta_omega_rotation(hkl) const;

	//возвращает пару omega, до которых нужно поверуть кристалл, чтобы быть в отражающем положении.
	std::vector<Gonio::OPC> omega_rotation(hkl) const;


	/*Задача трёхосного отражателя.
	Если закрепить phi и chi, то есть считать что они нам заданы -> получаем задачу одноосного отражателя с вращением 'omega'.
	В нашем решении был пункт, что если угол между осью z и вектором s(hkl) рефлекса не входит в промежуток (-th; th),
	Отражающего положения добиться невозможно. Далее нужно для всех phi подобрать такие chi, чтобы кристалл был в отражающем положении.
	или выяснить такие phi, для которых это невозможно.
	*/

	OPC psi0(const double ksi, hkl) const; //функция, возвращающая OPC для psi0 для минимального chi

	std::vector<OPC> diff_rotation(const double psi, const double ksi, hkl) const;

	/*
	
	запишем условие, что |sin Bt| >= |sin th| -> |z'| <= |cos th * |s||  -> -cos th*|s| <= z'<= cos th*|s|
	Ротация по phi и по chi. omega = 0!!!
	



	Попробуем посмотреть зависимость phi от chi
	z' = -sin(phi)*sin(chi)*x +cos(phi)*sin(chi)*y+cos(chi)*z
	(cos(chi)*z - z' )/sin(chi) = sin(phi)*[ x ] - cos(phi)*[ y ]
	[ cos(chi)*z - z']/(sin(chi)*sqrt(x^2 + y^2) ) = sin ( phi - КАКОЙ-ТО УГОЛ С ИЗВЕСТНЫМИ СИНУСОМ И КОСИНУСОМ )

	phi = УГОЛ   +    asin( [ cos(chi)*z - z']/( sin(chi)*sqrt( x^2 + y^2 ) ) )
	phi = УГОЛ + PI - asin( [ cos(chi)*z - z']/( sin(chi)*sqrt( x^2 + y^2 ) ) )

	phi определён если определён АРКСИНУС, то есть

	( cos(chi)*z - z' ) ^2 <= sin(chi)^2*(x^2 + y^2)
	(cos(chi)*z)^2 + ( z')^2 - 2cos(chi)*z*z' <= (x^2+y^2)( 1 - cos^2(chi))

	cos(chi)^2 * (s^2) - 2cos(chi)*(z*z') - (x^2+y^2 - (z')^2) <= 0    --- РЕШАЕМ

	cos(chi) ={ (z*z') +- sqrt( (z*z')^2 + (s^2)*(x^2+y*2 - (z')^2 ) }/s^2

	cos(chi) ={ (z*z') +- sqrt( (s^2 - (z')^2)*(x^2+y^2) ) }/ s^2

	ПОЛУЧАЕМ, что cos(chi) на пересечении интервалов [-1; 1] и 
	[ { (z*z') - sqrt( (s^2 - (z')^2)*(x^2+y^2) ) }/ s^2 ; { (z*z') + sqrt( (s^2 - (z')^2)*(x^2+y^2) ) }/ s^2 ]

	получаем, что chi = +-acos(  { (z*z') + sqrt( (s^2 - (z')^2)*(x^2+y^2) ) }/ s^2 ),
	либо 0, если аргумент аркосинуса больще 1!

	Далее находим phi и omega
	*/


	/*
	Далее нужно организовать поворот по psi
	
	Для этого повернём Сначала сделаем поворот на найденные нами (phi, chi, omega). 
	Получим рефлекс в отражающем положении при минимальном chi (psi = 0).
	Далее делаем поворот на Альфа(s') по оси z, на Бета(s') по оси y.
	Теперь делаем поворот на psi по оси х.
	Обратными преобразованиями на -Бета и -Альфа Получаем конечное положение вектора s''.

	Тогда конечную матрицу поворота от  (phi, chi, omega) можно записать как произведение матриц 
	(psi0, ch0, omega0)*Альфа*Бета*Пси*(-Бета)*(-Альфа) = RM (Rotation_Matrix) = (phi, chi, omega)
	
	Можно заметить, что cos(chi) = RM[3][3]
						tg(phi)  = - RM[1][3]/RM[2][3]
						th(omega) = RM[3][1]/RM[3][2]

						Тогда можно задать chi = +- acos(RM[3][3])
						Тогда задаётся phi и omega по знакам синуса и косинуса.
						Остальные уравнения можно использовать для проверки.
	
	*/
};



