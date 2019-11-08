﻿#pragma once


#include "R3.h"

/*
	Вывод дифракции в кристаллах
              ↑
               \         !!! принятое обозначение 2θ - угол между падающим и дифрагировавшим лучем
            θ / \
    ================================
            θ \  A  / θ    |
                /|\        | d
               / | \       |
    ==========/==+==\===============
             /   A   \
            /   /|\   \
           /   /   \   \
          /   / 2•θ'\   \  θ' = 90deg - θ
         ↓   ↓   |   ↑   ↑
        k_diff   ↓s  k_inc (волновой вектор падающего излучения)

	s - вектор обратного пространства, s = _a * h + _b * k + _c * l и |s| = 1/d

	разность хода лучей должна быть кратна/равна λ = 2d•cos θ'

	 -> d = λ/2cos θ', |s| = 1/d,
	 cos θ' = λ|s|/2

	 sin θ  = λ|s|/2
*/

//для более удобной и быстрой сокращенной записи 
#define hkl const double h, const double k, const double l

class Diff {
public:

	const static double lam;

	const static R3 k_inc;

	static R3 S(const matrix& _M_, hkl);

	static R3 k_diff(const R3& s);     //волновой вектор дифрагировавшего излучения.

	static double sin_th(const matrix& _M_, hkl);   //sin половинного угла между падающим и дифрагировавшим лучом

	static double sin_th(const R3& s);             //sin половинного угла между падающим и дифрагировавшем лучом

};



