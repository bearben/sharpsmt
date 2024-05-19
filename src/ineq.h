/*  ineq.h
 *
 *  Copyright (C) 2016-2024 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of sharpSMT.
 *  See COPYING for more information on using this software.
 */


#include <global.h>


#ifndef INEQ_HEADER
#define INEQ_HEADER

namespace volce{

class term{
public:
	//<index, multiplier>
	unsigned int 	id;
	double 			m;
	
	term() : id(0), m(1) {};
	term(const unsigned int index, const double multiplier) : id(index), m(multiplier) {};
	~term(){};
	
	// terms order
	const bool operator>(const term &elem) const{ return id > elem.id; };
	const bool operator<(const term &elem) const{ return id < elem.id; };
	
private:
};

//inequality core
class ineqc{
public:
	
	ineqc(const bool isequality, const double constant, std::vector<term> &terms) :
			eq(isequality), cst(constant), tm(terms) { unify(); };
	~ineqc(){};
	
	const bool iseq() const { return eq; };
	const bool isle() const { return !eq; };
	
	const double get_const() const { return cst; };
	const double get_const_r() const { return -cst; };
	
	const term operator[](unsigned int index) const { return tm[index]; };
	const unsigned int size() const { return tm.size(); };
	
	const std::vector<double> get_key() const;
	
	// sort & merge and keep the first term positive
	void unify();
	
private:

	// <const> + <term_0> + <term_1> + ... = (<=) 0
	bool eq;
	double cst;
	std::vector<term> tm;
	
	//sort and merge same id terms
	void merge_terms();
	
	//reverse all multipliers and the constant
	void reverse_mults();
	
};

}

#endif
