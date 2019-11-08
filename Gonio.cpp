//utf-8
#include "Gonio.h" 
#include "Diffraction.h"
#include "R3.h"


//Функция, возвращает угол в CФЕРИЧЕСКИХ координатах равный углу поворота вектора в плоскости ху в РАДИАНАХ
double Gonio::Alpha(const R3& r)
{
	return (r.y >= 0) ?  // в направлении tth
		ang({ r.x, r.y, 0 }, { 1, 0, 0 }) :
	    -ang({ r.x, r.y, 0 }, { 1, 0, 0 });
}


//Функция, возвращает угол в СФЕРИЧЕСКИХ координатах, равный углу между вектором и осью z в РАДИАНАХ
double Gonio::Beta(const R3& r)
{
	return acos(r.z / r.length());
}

//Основной Конструктор
Gonio::Gonio(const matrix& _M_) : _M_(_M_), _M_rotated(_M_) {}

//установить omega, phi, chi
Gonio& Gonio::set(const double omega, const double phi, const double chi) 
{
	//Порядок действий: устанавливаем phi, далее устанавливаем chi, далее устанавливаем omega!
	//Так как ориентация кристалла никак не зависит от порядка действий, мы выполняем поворот самым удобным для нас способом.

	//phi
	_M_rotated = _M_ * matrix({ { cos(phi), -sin(phi), 0.0 },
								{ sin(phi),  cos(phi), 0.0 },
								{      0.0,       0.0, 1.0 } });

	//chi
	_M_rotated = _M_rotated * matrix({ { 1.0,       0.0,      0.0 },
								       { 0.0,  cos(chi), sin(chi) },
								       { 0.0, -sin(chi), cos(chi) } });

	//omega
	_M_rotated = _M_rotated * matrix({ { cos(omega), -sin(omega), 0.0 },
							 	       { sin(omega),  cos(omega), 0.0 },
								       {        0.0,         0.0,  1.0} });

	this->phi = phi;
	this->chi = chi;
	this->omega = omega;

	return *this;
}

//устанавливаем угол omega, не меняя остальные
inline Gonio& Gonio::set_omega(const double omega)
{
	return set(omega, this->phi,this->chi);
}

//устанавливаем угол phi, не меняя остальные
inline Gonio& Gonio::set_phi(const double phi)
{
	return set(this->omega, phi, this->chi);
}

//устанавливаем угол chi, не меняя остальные
inline Gonio& Gonio::set_chi(const double chi)
{
	return set(this->omega, this->phi, chi);
}

//Возвращает базисные вектора после поворота. Так X' - нормаль к плоскости кольца, Y' - параллельна плоскости кольца, а Z' - параллельна оси держателя!
matrix Gonio::basic() const 
{   //в новом экземпляре класса вращаем ортонормированный базис
	return Gonio({ { 1,0,0 }, { 0,1,0 }, { 0,0,1 } }).set(omega, phi, chi)._M_rotated;
}


//Task 1

	/*Используя стереометрию, получаем, что если в сферических координатах s(hkl) задаётся углами "Al" и "Bt"
		и если угол между k_inc и s(hkl) равен (PI-'Gm'), то
                    ↑s(hkl)
                    /
             |   * /|.
             |*   / | .
             |Bt /  |  . Gm
             |  /  * \  .
             | / * Al \ .
____k_inc___\|/________\.______\   x
            /                  /

		cos Gm = sin Bt *( 1-2sin^2(Al/2) )

		//видно, что "Gm" меняется от PI/2 - "Bt" до PI/2 + "Bt"
		//если (th')  не входит в промежуток (PI/2-Bt; Pi/2+Bt) - отражающего положения не добиться изменением "alpha"
		//th принадлежит (-Bt; +Bt);
	*/

//Возвращает true, усли поворотом omega можно достичь отражающего положения
bool Gonio::is_omega_rotation_available(hkl) const
{
	return Diff::sin_th(Diff::S(_M_rotated, h, k, l)) <= sin(Beta(Diff::S(_M_, h, k, l)));
}

//Возвращает пару omega в РАДИАНАХ таких, что при их установке кристал находится в отражающем положении.
std::pair<double, double> Gonio::rotate_omega(hkl) const
{
	const R3 s = Diff::S(_M_rotated, h, k, l);
	double sin_Bt = sin(Beta(s));
	double sin_th = Diff::lam * s.length() / 2;
	//	s находится в отражающем положении, если Gm =PI-th', cos th' = sin th

	double Al_требуемое =/*+-*/ 2 * asin(sqrt(0.5 * (1 + sin_th / sin_Bt)));   //sin th = λ|s|/2

	double Al_начальное = Gonio::Alpha(s);

	//тогда определяем кратчайший угол поворота 

	return { Al_начальное - Al_требуемое, Al_начальное + Al_требуемое };
}




