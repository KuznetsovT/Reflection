//utf-8
#include "Gonio.h" 
#include "Diffraction.h"
#include "R3.h"

#include <iostream>

//Функция возвращает угол от -Pi до Pi по сторонам прямоугольного треугольника
double arc(const double cos, const double sin);

//Функция, возвращает угол в CФЕРИЧЕСКИХ координатах равный углу поворота вектора в плоскости ху в РАДИАНАХ по правилу правой руки
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

//Поворот производим по правилу правой руки, в Эйлеровой геометрии поворот по часовой стрелке.
/*
{ { 1.0,        0.0,       0.0 },
  { 0.0,  cos(xang), sin(xang) },
  { 0.0, -sin(xang), cos(xang) } });
*/
matrix Gonio::X_rotation(double xang)
{
	return matrix({ { 1.0,        0.0,       0.0 },
		    		{ 0.0,  cos(xang), sin(xang) },
					{ 0.0, -sin(xang), cos(xang) } });
}

/*
{ {  cos(yang), 0.0,-sin(yang) },
  {        0.0, 1.0,       0.0 },
  {  sin(yang), 0.0, cos(yang) } });
*/
matrix Gonio::Y_rotation(double yang)
{
	return matrix({ {  cos(yang), 0.0,-sin(yang) },
					{        0.0, 1.0,       0.0 },
					{  sin(yang), 0.0, cos(yang) } });
}

/*
{ {  cos(zang),-sin(zang), 0.0 },
  { -sin(zang), cos(zang), 0.0 },
  {        0.0,       0.0, 1.0 } });
*/
matrix Gonio::Z_rotation(double zang)
{
	return matrix({ {  cos(zang), sin(zang), 0.0 },
					{ -sin(zang), cos(zang), 0.0 },
					{        0.0,       0.0, 1.0 } });
}


matrix Gonio::OPC_rotation(OPC opc)
{
	return Z_rotation(-opc.phi) * X_rotation(-opc.chi) * Z_rotation(-opc.omega);
}

//Основной Конструктор
Gonio::Gonio(const matrix& _M_) : _M_(_M_), _M_rotated(_M_) {}

//установить omega, phi, chi
Gonio& Gonio::set(const Gonio::OPC& opc)
{
	//Порядок действий: устанавливаем phi, далее устанавливаем chi, далее устанавливаем omega!
	//Так как ориентация кристалла никак не зависит от порядка действий, мы выполняем поворот самым удобным для нас способом.

	_M_rotated = _M_ * OPC_rotation(opc);
	this->opc = opc;
	return *this;
}

//устанавливаем угол omega, не меняя остальные
Gonio& Gonio::set_omega(const double omega)
{
	return set({ omega, this->opc.phi,this->opc.chi });
}

//устанавливаем угол phi, не меняя остальные
Gonio& Gonio::set_phi(const double phi)
{
	return set({this->opc.omega, phi, this->opc.chi});
}

//устанавливаем угол chi, не меняя остальные
Gonio& Gonio::set_chi(const double chi)
{
	return set({ this->opc.omega, this->opc.phi, chi });
}

//Возвращает базисные вектора после поворота. Так X' - нормаль к плоскости кольца, Y' - параллельна плоскости кольца, а Z' - параллельна оси держателя!
matrix Gonio::basic() const 
{   //в новом экземпляре класса вращаем ортонормированный базис
	return Gonio({ { 1,0,0 }, { 0,1,0 }, { 0,0,1 } }).set(opc)._M_rotated;
}

matrix Gonio::M_rotated() const
{
	return this->_M_rotated;
}

R3 Gonio::S_rotated(hkl) const
{
	return Diff::S(_M_rotated, h,k,l);
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

		cos Gm = sin Bt * сos Al

		//видно, что "Gm" меняется от PI/2 - "Bt" до PI/2 + "Bt"
		//если (th')  не входит в промежуток (PI/2-Bt; Pi/2+Bt) - отражающего положения не добиться изменением "alpha"
		//th принадлежит (-Bt; +Bt);
	*/

//Возвращает true, eсли поворотом omega можно достичь отражающего положения
bool Gonio::is_omega_rotation_available(hkl) const
{
	return fabs(Diff::sin_th(Diff::S(_M_rotated, h, k, l))) <= fabs(sin(Beta(Diff::S(_M_, h, k, l))));
}



//Возвращает пару delta_omega в РАДИАНАХ таких, что при повороте на delta при данных phi и chi кристал находится в отражающем положении,
std::vector<Gonio::OPC> Gonio::delta_omega_rotation(hkl) const
{
	const R3 s = Diff::S(_M_rotated, h, k, l);
	double sin_Bt = sin(Beta(s));
	double sin_th = Diff::sin_th(_M_rotated, h, k, l);
	//	s находится в отражающем положении, если Gm =PI-th', cos th' = sin th
	// cos Al = (- sin th/sin Bt)

	//Al = +- acos(-sin_th/sin_Bt);
	double Al_требуемое =/*+-*/ acos(- sin_th / sin_Bt);   //sin th = λ|s|/2

	double Al_начальное = Gonio::Alpha(s);
               
	return { { (Al_начальное + (Al_требуемое)), opc.phi, opc.chi },
			 { (Al_начальное - (Al_требуемое)), opc.phi, opc.chi } };
}




//Возвращает, какой omega нужно установить, не меняя phi и chi, чтобы кристал был в отражающем положении.
std::vector<Gonio::OPC> Gonio::omega_rotation(hkl) const
{
	auto delta = delta_omega_rotation(h, k, l);
	for (unsigned i = 0; i < delta.size(); i++) delta[i].omega += opc.omega;
	return delta;
}


//Трёхмерная ротация, по данным углам psi, ksi и hkl, выдаём пару троек opc.
std::vector<Gonio::OPC> Gonio::diff_rotation(const double psi, const double ksi, hkl) const
{
	
	R3 s = Diff::S(_M_,h, k, l);
	OPC opc0 = psi0(ksi, h, k, l);
	R3 s1 = Gonio(_M_).set(opc0).S_rotated(h,k,l);

	//Матрица поворота Матрицы ориентации из начального положения.
	matrix RM = 
		OPC_rotation(opc0) * 
		Z_rotation(-Alpha(s1)) * Y_rotation(-Beta(s1)) * 
		Z_rotation(-psi) *                                      //вращение psi по часовой стрелке, как и phi
		Y_rotation(Beta(s1)) * Z_rotation(Alpha(s1));


	std::cout << std::endl << "|RM| : \n" << RM << std::endl;

	/*
	Можно заметить, что cos(chi) = RM[3][3]
		tg(phi) = -RM[1][3] / RM[2][3]
		th(omega) = RM[3][1] / RM[3][2]

		Тогда можно задать chi = +-acos(RM[3][3])
		Тогда задаётся phi и omega по знакам синуса и косинуса.
		Остальные уравнения можно использовать для проверки.
	*/
	std::vector<OPC> sol;
	for (auto chi : { acos(RM.c.z), -acos(RM.c.z) }) {

		double phi=arc(-RM.b.z/chi, RM.a.z/chi);

		double omega = arc(RM.c.y/chi, RM.c.x/chi);

		sol.push_back({ omega, phi, chi });

		matrix opc = OPC_rotation({ omega, phi, chi });

		std::cout << "\n |OPC| : \n" << opc << std::endl;

		double diff = 0;
		for (unsigned i = 0; i < 3; i++) {
			for (unsigned j = 0; j < 3; j++) {
				diff += abs(opc[i][j] - RM[i][j]);
			}
		}
		std::cerr << "\nDifferense between Rot Matrixes is : " << diff << std::endl << std::endl;
	}
	return sol;

}









//Однозначно восстанавливаем возможные phi по сhi и ksi.
std::vector<double> Phi_(const double chi,const double ksi, const matrix& m, hkl);




/*

........................ВЫВОД......................................


Попробуем посмотреть зависимость phi от chi
z' = sin(phi)*sin(chi)*x - cos(phi)*sin(chi)*y + cos(chi)*z
( z' - cos(chi)* z ) / sin(chi) = sin(phi)*[ x ] - cos(phi)*[ y ]
	[z' - cos(chi) * z]/sin(chi)*sqrt(x^2 + y^2 ) = sin ( phi - КАКОЙ-ТО УГОЛ С ИЗВЕСТНЫМИ СИНУСОМ И КОСИНУСОМ )

	phi = УГОЛ + asin([ z' - cos(chi) * z ]/( sin(chi)*sqrt( x^2 + y^2 ) ) )
	phi = УГОЛ + PI - asin([z' - cos(chi) * z]/( sin(chi)*sqrt( x^2 + y^2 ) ) )

	phi определён если определён АРКСИНУС, то есть

	(z' - cos(chi) * z ) ^2 <= sin(chi)^2*(x^2 + y^2)
	(cos(chi) * z) ^ 2 + (z')^2 - 2cos(chi)*z*z' <= (x ^ 2 + y ^ 2)(1 - cos ^ 2(chi))

	cos(chi) ^ 2 * (s ^ 2) - 2cos(chi) * (z * z') - (x^2+y^2 - (z') ^ 2) <= 0    -- - РЕШАЕМ

	cos(chi) = { (z * z') +- sqrt( (z*z') ^ 2 + (s ^ 2) * (x ^ 2 + y * 2 - (z')^2 ) }/s^2

	cos(chi) = { (z * z') +- sqrt( (s^2 - (z') ^ 2) * (x ^ 2 + y ^ 2)) } / s ^ 2

	ПОЛУЧАЕМ, что cos(chi) на пересечении интервалов[-1; 1] и
	[{ (z* z') - sqrt( (s^2 - (z') ^ 2)* (x ^ 2 + y ^ 2)) } / s ^ 2; { (z * z') + sqrt( (s^2 - (z') ^ 2)* (x ^ 2 + y ^ 2)) } / s ^ 2]

	получаем, что chi = +-acos({ (z * z') + sqrt( (s^2 - (z') ^ 2) * (x ^ 2 + y ^ 2) ) } / s ^ 2 ),
	либо 0, если аргумент аркосинуса больще 1!



	СМОТРИТЕ Phi_(...) и Psi0(...)
*/









//Находим тройку omega, phi, chi, соответствующую минимальному сhi/
Gonio::OPC Gonio::psi0(const double ksi, hkl) const {
	R3 s = Diff::S(_M_, h, k, l);
	R3 xy = { s.x, s.y, 0 };
	
	//debug
	//std::cout << "S_rotated should be : " << R3({ -s.length() * Diff::sin_th(s), s.length() * cos(Diff::th(s)) * cos(ksi), s.length() * cos(Diff::th(s)) * sin(ksi) }) << std::endl;

	/*
	Найдём z' от ksi.
	
			| z'| = | sin ksi * cos th * |s| |
	*/

	//найдём требуемый Z, он не зависит от Omega, только от phi и chi
	double z1 = sin(ksi) * abs( cos(Diff::th(s)) * s.length() );

	//return!
	OPC opc0;

	//Расскрывая матрицу преобразования, получаем
	//z' = sin(phi)*sin(chi)*x -cos(phi)*sin(chi)*y+cos(chi)*z
	//
	//Выводим зависимость phi от chi, смотрим ОДЗ,
	//
	//

	//Условие, что полученный максимальный сos(chi) не пинает chi в комплексную плоскость
	if (abs((s.z * z1) + sqrt(((s ^ s) - z1 * z1) * (xy ^ xy))) > (s ^ s)) opc0.chi = 0;
	else opc0.chi = acos(((s.z * z1) + sqrt(((s ^ s) - z1*z1) * (xy ^ xy))) / (s ^ s));

	//std::cout << "OPC0 chi min : " << opc0.chi << std::endl;
	opc0.phi = Phi_(opc0.chi, ksi, _M_, h, k, l)[0];

	//std::cout << "z1 should be : " << z1 << " and z1 is " << s.x * sin(opc0.phi) * sin(opc0.chi) - s.y * cos(opc0.phi) * sin(opc0.chi) + s.z * cos(opc0.chi) << std::endl;

	auto opc1 =  Gonio(_M_).set(opc0).omega_rotation(h, k, l)[0];
	//std::cout << " OPC0 [ omega: " << opc1.omega << " phi: " << opc1.phi << " chi: " << opc1.chi << " ]\n";


	//std::cout << "Calculated S1 for PSI 0 : " << Gonio(_M_).set(opc1).S_rotated(h, k, l) << std::endl;

	return opc1;
}






//Восстанавливаем phi по chi: phi = УГОЛ + PI - asin([z' - cos(chi) * z]/( sin(chi)*sqrt( x^2 + y^2 ) ) )
std::vector<double> Phi_(const double chi,const double ksi, const matrix& m, hkl)
{
	R3 s = Diff::S(m, h, k, l);
	R3 xy = { s.x, s.y, 0 };
	double z1 = sin(ksi) * abs(cos(Diff::th(m, h, k, l)) * s.length());

	double ang = arc(s.x, s.y);


	//Возникла проблема с обработкой значений очень близких к Pi/2
	double arcsin = (abs((z1 - cos(chi) * s.z) / (sin(chi) * xy.length()) < 1))?
		asin((z1 - cos(chi) * s.z) / (sin(chi) * xy.length())):
		PI/2;

	return { ang + arcsin, ang + PI - arcsin };
}



//Возвращает угол по его катетам.
double arc(const double cos, const double sin)
{
	return (sin > 0) ? acos(cos / sqrt(sin * sin + cos * cos)) : -acos(cos / sqrt(sin * sin + cos * cos));
}
