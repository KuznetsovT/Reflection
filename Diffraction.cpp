#include "Diffraction.h"


const double Diff::lam = 0.71073;      //Длина вволны рентгеновского источника, взвешанная MoKa = 2/3 MoKa1 + 1/3 MoKa2, в Ангстремах

const R3 Diff::k_inc = { 1 / lam, 0.0, 0.0 };    // волновой вектор падающего излучения

R3 Diff::S(const matrix& _M_, hkl)
{
	return _M_.a * h + _M_.b * k + _M_.c * l;
}

R3 Diff::k_diff(const R3& s)
{
	return k_inc+s;
}

double Diff::sin_th(const matrix& _M_, hkl)
{
	return sin_th(S(_M_, h, k, l));
}

double Diff::sin_th(const R3& s)
{
	return lam * s.length() / 2;
}


