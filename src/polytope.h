/*  polyvest.h
 *
 *  Copyright (C) 2015-2024 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of PolyVest.
 *  See COPYING for more information on using this software.
 */


#include <iostream>
#include <algorithm>
#include "armadillo"
#include "math.h"
#include "time.h"
#include "memory.h"

#ifndef POLYVOL_H
#define POLYVOL_H

namespace polyvest{

class polytope{
public:
	polytope(int rows, int cols);
	~polytope();
	double 	matA(double val, int i, int j){ return (A(i, j) = val); }
	double 	matA(int i, int j){ return A(i, j); }
	double 	vecb(double val, int i){ return (b(i) = val); }
	double 	vecb(int i){ return b(i); }
	
	void	GetBounds();	// results in xmax and xmin
	void	GetBounds(int index);
	bool	Enlarge();
	bool 	AffineTrans();
	double 	EstimateVol(double epsilon, double delta, double coef);
	double	ExactVol();		// call vinci

	void	PrepForWalk();
	void	Walk();
	arma::vec GetInvPoint(arma::vec point);
	bool	isInside(arma::vec point);
	bool	CUonBoundary(arma::vec point);
	
	double 	Volume() const { return vol; }
	void 	Print();
	polytope* Clone();

	bool 	msg_off;

	//reciprocal of beta, beta-cut
	double beta_r;
	
	//polytope denoted by: Ax<=b, A is an (m x n) matrix.
	int 	m, n;
	arma::mat A;
	arma::vec b;
	int 	*var_flag;	// 0 = Real; 1 = Integer

	// random walk point
	arma::vec x;
	
	// upper & lower bounds for each xi
	double 	*xmax;
	double 	*xmin;
	double 	*hpoffset;
	
private:
	double 	walk(int k);
	void 	genInitE(double &R2, arma::vec &Ori);

	double 	randd(double u){ return rand() * u / RAND_MAX; }
	int 	randi(int u){ return rand() % u; }

	arma::mat invT;		// to compute invert points
	arma::vec invOri;
	
	double 	vol;		// volume and det(A)
	double	determinant;
	int 	l;
	double 	*r2;

	arma::vec *B;
	arma::mat *Ai;
};

inline polytope::polytope(int rows, int cols) :
	msg_off(true),
	m(rows),
	n(cols),
	A(rows, cols),
	b(rows),
	x(cols),
	invT(cols, cols),
	invOri(cols), 
	vol(0),
	determinant(0)
{
	srand((int)time(0));

	beta_r = 2 * n; //2 * n;

	l = (int)(n * log((double)beta_r) / log((double)2)) + 2;
	r2 = new double[l];
	for (int i = 0; i < l; i++) 
		r2[i] = pow((double)2, (double)(2 * i) / n);

	var_flag = new int[n] ();

	xmax = new double[n] ();
	xmin = new double[n] ();
	hpoffset = new double[m] ();

	B = new arma::vec[n];
	Ai = new arma::mat[n];
}

inline polytope::~polytope(){
	delete []r2;
	delete []var_flag;
	delete []xmax;
	delete []xmin;
	delete []hpoffset;
	delete []B;
	delete []Ai;
}

}

#endif
