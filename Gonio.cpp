//utf-8
#include "Gonio.h" 
#include "Diffraction.h"
#include "R3.h"

#include <iostream>

//Функция возвращает угол от -Pi до Pi по сторонам прямоугольного треугольника
double arc(const double x, const double y);

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
	if (r.length() == 0) return 0;
	else return acos(r.z / r.length());
}

//Поворот производим по правилу правой руки, в Эйлеровой геометрии поворот по часовой стрелке.
/*
     x↓          y↓         z↓
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

/*      x↓       y↓      z↓
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

/*      x↓         y↓       z↓
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

/*
{ cos phi, -sin phi, 0}   { 1,       0,        0}   {  cos omega, sin omega, 0}
{ sin phi,  cos phi, 0} * { 0, cos chi, -sin chi} * { -sin omega, cos omega, 0}
{       0,        0, 1}   { 0, sin chi,  cos chi}   {          0,         0, 1}
*/
matrix Gonio::OPC_matrix(OPC opc)
{
	return Z_rotation(-opc.phi) * X_rotation(-opc.chi) * Z_rotation(opc.omega);
}

//Основной Конструктор
Gonio::Gonio(const matrix& _M_) : _M_(_M_), _M_rotated(_M_) {}

//установить omega, phi, chi
Gonio& Gonio::set(const Gonio::OPC& opc)
{
	//Порядок действий: устанавливаем phi, далее устанавливаем chi, далее устанавливаем omega!
	//Так как ориентация кристалла никак не зависит от порядка действий, мы выполняем поворот самым удобным для нас способом.

	_M_rotated = _M_ * OPC_matrix(opc);
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

R3 Gonio::S_rotated(const HKL &hkl) const
{
	return Diff::S(_M_rotated, hkl);
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
		Gm = 90 deg + th
		//видно, что "Gm" меняется от PI/2 - "Bt" до PI/2 + "Bt"
		//если (th')  не входит в промежуток (PI/2-Bt; Pi/2+Bt) - отражающего положения не добиться изменением "alpha"
		//th принадлежит (-Bt; +Bt);
	*/

//Возвращает true, eсли поворотом omega можно достичь отражающего положения
bool Gonio::is_omega_rotation_available(const HKL &hkl) const
{
	return abs(Diff::sin_th(Diff::S(_M_rotated, hkl))) <= abs(sin(Beta(Diff::S(_M_, hkl))));
}



//Возвращает пару delta_omega в РАДИАНАХ таких, что при повороте на delta при данных phi и chi кристал находится в отражающем положении,
std::vector<Gonio::OPC> Gonio::delta_omega_rotation(const HKL &hkl) const
{
	const R3 s = Diff::S(_M_rotated, hkl);
	double sin_th = Diff::sin_th(_M_rotated, hkl);
	double sin_Bt = sin(Beta(s));
	//	s находится в отражающем положении, если Gm =PI-th', cos th' = sin th
	// cos Al = (- sin th/sin Bt)

	//Al = +- acos(-sin_th/sin_Bt);
	double Al_требуемое =/*+-*/ acos(- sin_th / sin_Bt);   //sin th = λ|s|/2
	if (abs(sin_th) > abs(sin_Bt) && abs(sin_th) - abs(sin_Bt) <= 20 * DBL_EPSILON) {
		Al_требуемое = acos(-1);
	}

	double Al_начальное = Gonio::Alpha(s);
               
	return { { (Al_требуемое - Al_начальное), opc.phi, opc.chi },
			 {-(Al_требуемое + Al_начальное), opc.phi, opc.chi } };
}




//Возвращает, какой omega нужно установить, не меняя phi и chi, чтобы кристал был в отражающем положении.
std::vector<Gonio::OPC> Gonio::omega_rotation(const HKL &hkl) const
{
	auto delta = delta_omega_rotation(hkl);
	for (unsigned i = 0; i < delta.size(); i++) delta[i].omega += opc.omega;
	return delta;
}










//??????????????????????????????????????????????????????????????????????????????????????????
//??????????????????????????????????????????????????????????????????????????????????????????
//??????????????????????????????????????????????????????????????????????????????????????????
//                                ТРЁХМЕРНЫЙ ОТРАЖАТЕЛЬ

/*
Данная задача решается в несколько шагов.

1. При данном ksi однозначно восстанавливается z'. Находим z'. 
	Далее хотим найти тройку opc, чтобы кристалл был в отражающем положении,
	и чтобы chi был минимален - это соответствует psi = 0.
	Расскрываем матрицу поворота, видим, что z' зависит только от phi и chi. 
	Выражаем phi через chi, и определяем минимально возможный chi.

........................ВЫВОД......................................
					R' = R * OPC_matrix

	z' = sin(phi)*sin(chi)*x - cos(phi)*sin(chi)*y + cos(chi)*z
	( z' - cos(chi)* z ) / sin(chi) = sin(phi)*[ x ] - cos(phi)*[ y ]
	[z' - cos(chi) * z]/sin(chi)*sqrt(x^2 + y^2 ) = sin ( phi - КАКОЙ-ТО УГОЛ С ИЗВЕСТНЫМИ СИНУСОМ И КОСИНУСОМ )

	phi = УГОЛ + asin([ z' - cos(chi) * z ]/( sin(chi)*sqrt( x^2 + y^2 ) ) )                            (1)
	phi = УГОЛ + PI - asin([z' - cos(chi) * z]/( sin(chi)*sqrt( x^2 + y^2 ) ) )

	phi определён если определён АРКСИНУС, то есть

	(z' - cos(chi) * z ) ^2 <= sin(chi)^2*(x^2 + y^2)
	(cos(chi) * z) ^ 2 + (z')^2 - 2cos(chi)*z*z' <= (x ^ 2 + y ^ 2)(1 - cos ^ 2(chi))

	cos(chi) ^ 2 * (s ^ 2) - 2cos(chi) * (z * z') - (x^2+y^2 - (z') ^ 2) <= 0    ---- РЕШАЕМ КВАДРАТНОЕ УРАВНЕНИЕ
	cos(chi) = { (z * z') +- sqrt( (z*z') ^ 2 + (s ^ 2) * (x ^ 2 + y * 2 - (z')^2 ) }/s^2
	cos(chi) = { (z * z') +- sqrt( (s^2 - (z') ^ 2) * (x ^ 2 + y ^ 2)) } / s ^ 2

	получаем, что chi = +-acos({ (z * z') + sqrt( (s^2 - (z') ^ 2) * (x ^ 2 + y ^ 2) ) } / s ^ 2 ),    (2)
	либо 0, если аргумент аркосинуса больще 1!

	......................................................................


2. По формуле (2) в функции Psi0(...) вычисляем chi0 - минимально возможный chi.

3. Далее по формуле (1) находим phi, соответствующий минимально возможному chi.

4. В конце определяется omega, соответствующая phi и chi, данное определяем при помощи одноосного отражателя.

5. Мы нашли opc, соответствующее нулевому psi. Далее записываем матрицу поворота(RM), если бы мы повернули s до совпадения с осью z,
	повернули s на угол psi, вернули обратно, и далее повернули в отражающее положение. 

	Данная матрица поворота должна соответствовать матрице поворота, записанной в Эйлеровой геометрии, находим omega, phi, chi через матрицу поворота.
	
	Например,
	Можно заметить, что cos(chi) = RM[3][3]
		tg(phi) = -RM[1][3] / RM[2][3]
		th(omega) = -RM[3][1] / RM[3][2]

		Тогда можно задать chi = +-acos(RM[3][3])
		Тогда задаётся phi и omega по знакам синуса и косинуса.
		Остальные уравнения можно использовать для проверки.


	Полученные углы - искомые.

*/




//Трёхмерная ротация, по данным углам psi, ksi и hkl, выдаём пару троек opc.
std::vector<Gonio::OPC> Gonio::diff_rotation(const double psi, const double ksi, const HKL&hkl) const
{
	
	R3 s = Diff::S(_M_,hkl);
	OPC opc0 = psi0(ksi, hkl);

	//Матрица поворота Матрицы ориентации из начального положения.
	matrix RM =
		Z_rotation(-Alpha(s)) * Y_rotation(-Beta(s)) *
		Z_rotation(-psi) *                                      //вращение psi по часовой стрелке, как и phi
		Y_rotation(Beta(s)) * Z_rotation(Alpha(s)) *
		OPC_matrix(opc0);

	
	std::vector<OPC> sol;
	for (auto chi : { acos(RM.c.z), -acos(RM.c.z) }) {
		
		double phi=0;
		double omega=0;

		//дополнительно нужно рассмотреть случай, когда sin(chi) == 0
		
		if (abs(sin(chi)) <= 20*DBL_EPSILON) {
			phi = /*ANY NUMBER, но мы возьмём */0;
		}
		else {
			phi = (sin(chi) >= 0) ? arc(-RM.b.z, RM.a.z) : arc(+RM.b.z, -RM.a.z);
		}





		if (abs(sin(chi)) <= 20*DBL_EPSILON) {
			omega = arc(RM.a.x, RM.a.y);
		}
		else {
			omega = (sin(chi) >= 0) ? arc(RM.c.y, -RM.c.x) : arc(-RM.c.y, +RM.c.x);
		}

		

		sol.push_back({ omega, phi, chi });
		/*
		//При проверке показывается, что матрицы действительно совпадают с точностью до точности вычислений
		using namespace std;
		cout << "\n sin( chi ) : " << sin(chi) << endl;
		cerr << "__RM__:\n" << RM << endl << endl;
		cerr << "__OPC_:\n" << OPC_matrix({ omega, phi, chi }) << endl << endl;
		*/
	}
	return sol;

}


//Однозначно восстанавливаем возможные phi по сhi и ksi.
std::vector<double> Phi_(const double chi,const double ksi, const R3 & s);



//Находим тройку omega, phi, chi, соответствующую минимальному сhi/
Gonio::OPC Gonio::psi0(const double ksi,const HKL& hkl) const {
	R3 s = Diff::S(_M_, hkl);
	R3 xy = { s.x, s.y, 0 };
	
	//Найдём z' от ksi.
	//		| z'| = | sin ksi * cos th * |s| |

	//найдём требуемый Z, он не зависит от Omega, только от phi и chi
	double z1 = sin(ksi) * abs( cos(Diff::th(s)) * s.length() );

	//return!
	OPC opc0 = { 0, 0, 0 };


	//Дальнейшее решение рассматривает случай, когда x != 0 || y != 0. Запишем отдельно решение для этого случая

	if (s.x == 0 && s.y == 0) {
		opc0.chi = acos(z1 / s.z);

		//phi в данном случае можно брать любой!                            
		opc0.phi = 0;                                        //           ←---------ANY NUMBER
		return Gonio(_M_).set(opc0).delta_omega_rotation(hkl)[(cos(ksi) >=0)?0:1];
	}

	//Расскрывая матрицу преобразования, получаем
	//z' = sin(phi)*sin(chi)*x -cos(phi)*sin(chi)*y+cos(chi)*z
	//Выводим зависимость phi от chi, смотрим ОДЗ,

	//Условие, что полученный максимальный сos(chi) не пинает chi в комплексную плоскость
	if (abs(((s.z * z1) + sqrt(((s ^ s) - z1 * z1) * (xy ^ xy))) / (s ^ s)) >= 1) {

		opc0.chi = 0;

		//phi в данном случае можно брать любой!
		opc0.phi = 0;                                       //           ←---------ANY NUMBER
	}
	else {

		opc0.chi = acos(((s.z * z1) + sqrt(((s ^ s) - z1 * z1) * (xy ^ xy))) / (s ^ s));
		opc0.phi = Phi_(opc0.chi, ksi, s)[0];
	}

	

	return Gonio(_M_).set(opc0).delta_omega_rotation(hkl)[(cos(ksi) >= 0) ? 0 : 1];
	
}






//Восстанавливаем phi по chi: phi = УГОЛ + PI - asin([z' - cos(chi) * z]/( sin(chi)*sqrt( x^2 + y^2 ) ) )
std::vector<double> Phi_(const double chi,const double ksi, const R3 & s)
{
	R3 xy = { s.x, s.y, 0 };
	//if (s.x == 0 && s.y == 0) return { 0 }; //-------!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!ПРОРАБОТАТЬ!!!

	double z1 = sin(ksi) * abs(cos(Diff::th(s)) * s.length());

	double ang = arc(s.x, s.y);

	//Возникла проблема с обработкой значений очень близких к Pi/2
	double _sin_ = (z1 - cos(chi) * s.z) / (sin(chi) * xy.length());
	double arcsin = (abs(_sin_) < 1) ?
		asin(_sin_) :
		((_sin_ > 0) ? PI / 2 : -PI / 2);

	return { ang + arcsin, ang + PI - arcsin };
}



//Возвращает угол по его катетам.
double arc(const double x, const double y)
{
	return (y >= 0) ? acos(x / sqrt(y*y + x*x)) : -acos(x / sqrt(y*y + x*x));
}
