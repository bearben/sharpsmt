/*  print.cpp
 *
 *  Copyright (C) 2016-2024 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of sharpSMT.
 *  See COPYING for more information on using this software.
 */
 

#include <solver.h>


/*
	Print in treeview
*/
void volce::solver::print_ast(const dagc root) {

	if (root.iscbool()) {
		//boolean constant
		if (root.bval()) std::cout << "true";
		else std::cout << "false";
	} else if (root.iscnum()) {
		//constant number
		std::cout << root.nval();
	} else if (root.isvbool()) {
		//boolean variable
		if (root.isneg()) std::cout << "[N]";
		std::cout << vbool_list.name(root.id);
	} else if (root.isineq()) {
		//inequality
		if (root.isneg()) std::cout << "[N]";
		std::cout << "_ie" << root.id;
	} else if (root.isvnum()) {
		//numeric variable
		if (root.m != 1) std::cout << '[' << root.m << ']';
		std::cout << vnum_list.name(root.id);
	} else {

		std::cout << '(';

		if (root.isbool() && root.isneg())
			std::cout << "[N]";
		else if (!root.isbool()) {
			if (root.m != 1) std::cout << "[m:" << root.m << ']';
			if (root.v != 0) std::cout << "[v:" << root.v << ']';
		}
			

		if (root.isand()) 
			std::cout << "and";
		else if (root.isor()) 
			std::cout << "or";
		else if (root.iscomp()) 
			std::cout << "=";
		else if (root.isite()) 
			std::cout << "ite";
		else if (root.isadd()) 
			std::cout << "+";
		else if (root.ismul()) 
			std::cout << "*";
		else if (root.isdiv()) 
			std::cout << "/";
		else assert(false);
	
		const std::vector<dagc> c = bop_list[root.id];
		for (unsigned int i = 0; i < c.size(); i++) {
			std::cout << ' ';
			print_ast(c[i]);
		}
	
		std::cout << ')';
		
	}

}


/*
	print model (bunch)
*/
void volce::solver::print_model() {

	std::cout << "(and" << std::endl;

	for (unsigned int i = 0; i < vbool_list.size(); i++) {
		if (vbool_list(i).is_true())
			std::cout << "  (= " << vbool_list.name(i) << " true)" << std::endl;
		else if (vbool_list(i).is_false())
			std::cout << "  (= " << vbool_list.name(i) << " false)" << std::endl;
	}
	for (unsigned int i = 0; i < ineq_list.size(); i++) {
		// print ineq with assignment (true or false)
		if (ineq_list(i).is_true()) {
			// is true
			std::cout << "  (= " << ineq_list.name(i) << " true)" << std::endl;
		} else if (ineq_list(i).is_false()) {
			// is false
			std::cout << "  (= " << ineq_list.name(i) << " false)" << std::endl;
		}
	}
	std::cout << ')';

}


/*
	print inequality
*/
void volce::solver::print_ineq(unsigned int index) {

	volce::ineqc ie = ineq_list[index];
	
	std::cout << "(= " << ineq_list.name(index) << ' ';
	
	if (ie.iseq()) std::cout << "(=";
	else std::cout << "(<=";
	
	if (ie.size() > 1) std::cout << " (+";
	
	for (unsigned int j = 0; j < ie.size(); j++)
		if (ie[j].m == 1) std::cout << ' ' << vnum_list.name(ie[j].id);
		else {
			std::cout << " (* ";
			if (ie[j].m < 0) std::cout << "(- " << -ie[j].m << ')';
			else if (ie[j].m == 0) std::cout << 0;
			else std::cout << ie[j].m;
			std::cout << ' ' << vnum_list.name(ie[j].id) << ')';
		}
		
	if (ie.size() > 1) std::cout << ')';
	
	if (ie.get_const_r() < 0) std::cout << " (- " << -ie.get_const_r() << ")))";
	else if (ie.get_const_r() == 0) std::cout << " 0))";
	else std::cout << ' ' << ie.get_const_r() << "))";

}



