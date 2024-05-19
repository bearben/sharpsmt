/*  solve.cpp
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
	SMT solving
*/

//initialize z3 solver
void volce::solver::z3_init() {

	//mark solver initialized, cannot add constraints any more
	solving_initialized = true;
	
	//init volume routine
	vol_init();
	
	//for (unsigned int i = 0; i < ineq_list.size(); i++) {
		//print_ineq(i);
		//std::cout << std::endl;
	//}
	
	//make z3 variables, bool, numeric and ineqs
	for (unsigned int i = 0; i < vbool_list.size(); i++)
		vbool_expr.push_back(z3context.bool_const(vbool_list.name(i).c_str()));

	for (unsigned int i = 0; i < vnum_list.size(); i++)
		if (islia())
			vnum_expr.push_back(z3context.int_const(vnum_list.name(i).c_str()));
		else
			vnum_expr.push_back(z3context.real_const(vnum_list.name(i).c_str()));
			
	for (unsigned int i = 0; i < ineq_list.size(); i++)
		ineq_expr.push_back(z3context.bool_const(ineq_list.name(i).c_str()));
	
	//make z3 operators
	for (unsigned int i = 0; i < bop_list.size(); i++)
		op_expr.push_back(z3_mk_op(i));
	
	//add assertions
	for (unsigned int i = 0; i < assert_list.size(); i++)
		z3solver.add(z3_mk_expr(assert_list[i]));
	
	//add inequalities
	for (unsigned int i = 0; i < ineq_list.size(); i++)
		z3solver.add(z3_mk_ineq(i));
		
	z3_init_bounds(wordlength);
		
}

//more initializing operations
void volce::solver::z3_init_bounds(const int wordlength) {

	// disable bounds
	if (wordlength <= 0) return;

	//add bounds for variables
	for (unsigned int i = 0; i < vnum_list.size(); i++) {
		z3solver.add(vnum_expr[i] <= (int)pow(2, wordlength - 1) - 1);
		z3solver.add(vnum_expr[i] >= -(int)pow(2, wordlength - 1));
	}
	
}


z3::expr volce::solver::z3_mk_ineq(const unsigned int index) {
	
	ineqc ie = ineq_list[index];
	assert(ie.size() > 0);
	
	//make first term
	z3::expr expr = z3_mk_term(ie[0]);
	
	//add all terms
	for (unsigned int i = 1; i < ie.size(); i++)
		expr = expr + z3_mk_term(ie[i]);
	
	if (ie.iseq()) 
		expr = (expr == z3_mk_nconst(ie.get_const_r()));
	else
		expr = (expr <= z3_mk_nconst(ie.get_const_r()));
	
	return (ineq_expr[index] == expr);

}


z3::expr volce::solver::z3_mk_term(const term t) {
	if (t.m == 1) 
		// 1 * x = x
		return vnum_expr[t.id];
	else 
		return z3_mk_nconst(t.m) * vnum_expr[t.id];
}


z3::expr volce::solver::z3_mk_nconst(const double val) {

//	if (islia()) {
		//int const
//		return z3context.int_val((int)val);
//	} else {
		//real const
		char num[STRLEN];
		sprintf(num, "%.8g", val);
		//remove char '+'
		for (unsigned int i = 0; num[i] != 0; i++)
			if (num[i] == '+')
				while (num[i] != 0)
					num[i] = num[i + 1], i++;
		return z3context.real_val(num);
//	}

}


z3::expr volce::solver::z3_mk_op(const unsigned int index) {

	const dagc node(bop_list.type(index), index);
	const std::vector<dagc> c = bop_list[index];

	if (node.isand()) {
		//operator AND
		z3::expr expr = z3_mk_expr(c[0]);
		for (unsigned int i = 1; i < c.size(); i++)
			expr = (expr && z3_mk_expr(c[i]));
		return expr;
	} else if (node.isor()) {
		//operator OR
		z3::expr expr = z3_mk_expr(c[0]);
		for (unsigned int i = 1; i < c.size(); i++)
			expr = (expr || z3_mk_expr(c[i]));
		return expr;
	} else if (node.iscomp()) {
		//operator EQ
		return z3_mk_expr(c[0]) == z3_mk_expr(c[1]);
	} else if (node.isitebool()) {
		//operator ITE
		return z3::ite(z3_mk_expr(c[0]), z3_mk_expr(c[1]), z3_mk_expr(c[2]));
	} else {
		//others
		assert(false);
	}
}


z3::expr volce::solver::z3_mk_expr(const dagc node) {

	if (node.iscbool()) {
		//constant bool
		return z3context.bool_val(node.bval());
	} else if (node.isvbool()) {
		//boolean variable
		return (node.isneg() ? !vbool_expr[node.id] : vbool_expr[node.id]);
	} else if (node.isineq()) {
		//introduced ineq variable
		return (node.isneg() ? !ineq_expr[node.id] : ineq_expr[node.id]);
	} else if (node.isboolop()) {
		//operator
		return (node.isneg() ? !op_expr[node.id] : op_expr[node.id]);
	} else assert(false);

}


// repeatable smt solving
// true: achieve a feasible bunch
// false: unsat
const bool volce::solver::solve() {

	// check, return false if unsat
	if (z3solver.check() != z3::sat) return false;
	
	// extract the model
	z3::model z3model = z3solver.get_model();
	
	//init with -1
	ineq_list.init_vals();
	vbool_list.init_vals();
	
	// scan and extract the partial assignments
	for (unsigned int i = 0; i < z3model.size(); i++){
		z3::func_decl var = z3model[i];
		
		// model presented by z3 doesnt contain unknown assignments
		// first, check whether it is an inequality by name
		unsigned int id = ineq_list.find(var.name().str());
		if (id != ineq_list.size()) {
			// inequalities
			// set ineq val
			ineq_list(id) = dagv(eq(z3model.get_const_interp(var), z3context.bool_val(true)));
		} else {
			// then, check whether it is a pure bool variable by name
			id = vbool_list.find(var.name().str());
			if (id != vbool_list.size())
				vbool_list(id) = dagv(eq(z3model.get_const_interp(var), z3context.bool_val(true)));
		}
	}

	if (enable_bunch) {

		// try to reduce the solution into a bunch
		std::vector<bool> ineq_flip(ineq_list.size(), true);
		for (unsigned int i = 0; i < ineq_list.size(); i++)
			if (ineq_list(i).is_unknown()) 
				ineq_flip[i] = false;
		std::vector<bool> vbool_flip(vbool_list.size(), true);
		for (unsigned int i = 0; i < vbool_list.size(); i++)
			if (vbool_list(i).is_unknown()) 
				vbool_flip[i] = false;
		
GOTO_BUNCH_CYCLE:
	
		get_flip_list(ineq_flip, vbool_flip);

		for (unsigned int i = 0; i < vbool_list.size(); i++) {
			// skip vbool should not be flipped
			if (!vbool_flip[i]) continue;

			vbool_list(i).negate();
			if (get_result()) {
				// reduce i-th vbool and start next round
				vbool_list(i) = dagv();
				vbool_flip[i] = false;
				goto GOTO_BUNCH_CYCLE;
			} else {
				vbool_list(i).negate();
			}
		}
		
		for (unsigned int i = 0; i < ineq_list.size(); i++) {
			// skip ineq should not be flipped
			if (!ineq_flip[i]) continue;
		
			ineq_list(i).negate();
			if (get_result()) {
				// reduce i-th ineq and start next round
				ineq_list(i) = dagv();
				ineq_flip[i] = false;
				goto GOTO_BUNCH_CYCLE;
			} else {
				ineq_list(i).negate();
			}
		}

		// no more var can be flipped
	
	}
	
/*
	// since the variables in the condition of ITEs are in front of the variables in branch statements
	// we select variables from the back and give the priority to the inequality-vairables
	int iid = ineq_list.size() - 1;
	int bid = vbool_list.size() - 1;
	std::vector<bool> ineq_negated(ineq_list.size(), false);
	std::vector<bool> vbool_negated(vbool_list.size(), false);
	while (true)
		if (iid >= 0) {
			if (ineq_list(iid).is_unknown() || ineq_negated[iid]) {
				//skip unknown and negated var
				iid--;
			} else {
				ineq_list(iid).negate();
				if (get_result()) {
					//reduce ineq(iid) and start next round
					ineq_negated[iid] = true;
					ineq_list(iid) = dagv();
					iid = ineq_list.size() - 1;
					bid = vbool_list.size() - 1;
				} else {
					ineq_list(iid).negate();
					iid--;
				}
			}
		} else if (bid >= 0) {
			if (vbool_list(bid).is_unknown() || vbool_negated[bid]) {
				//skip unknown and negated var
				bid--;
			} else {
				vbool_list(bid).negate();
				if (get_result()) {
					//reduce vbool(bid) and start next round
					vbool_negated[bid] = true;
					vbool_list(bid) = dagv();
					iid = ineq_list.size() - 1;
					bid = vbool_list.size() - 1;
				} else { 
					vbool_list(bid).negate();
					bid--;
				}
			}
		} else break;
*/
/*
	// select variables with standard sequence
	unsigned int iid = 0;
	unsigned int bid = 0;
	std::vector<bool> ineq_negated(ineq_list.size(), false);
	std::vector<bool> vbool_negated(vbool_list.size(), false);
	while (true)
		if (bid < vbool_list.size()) {
			if (vbool_list(bid).is_unknown() || vbool_negated[bid]) {
				//skip unknown and negated var
				bid++;
			} else {
				vbool_list(bid).negate();
				if (get_result()) {
					//reduce vbool(bid) and start next round
					vbool_negated[bid] = true;
					vbool_list(bid) = dagv();
					iid = bid = 0;
				} else { 
					vbool_list(bid).negate();
					bid++;
				}
			}
		} else if (iid < ineq_list.size()) {
			if (ineq_list(iid).is_unknown() || ineq_negated[iid]) {
				//skip unknown and negated var
				iid++;
			} else {
				ineq_list(iid).negate();
				if (get_result()) {
					//reduce ineq(iid) and start next round
					ineq_negated[iid] = true;
					ineq_list(iid) = dagv();
					iid = bid = 0;
				} else {
					ineq_list(iid).negate();
					iid++;
				}
			}
		} else break;
*/
/*
	print_model();
	std::cout << std::endl;
	for (unsigned int i = 0; i < assert_list.size(); i++) {
		print_ast(assert_list[i]);
		std::cout << std::endl;
	}
	for (unsigned int i = 0; i < ineq_list.size(); i++) {
		print_ineq(i);
		std::cout << std::endl;
	}
	std::cout << std::endl << std::endl;
*/	
	// store bunches
	bunch_elem bunch(ineq_list, vbool_list);
	bunch_list.push_back(bunch);
	
	// multipliers = 2^(the number of unassigned bools)
	unsigned int m = 1;
	for (unsigned int i = 0; i < vbool_list.size(); i++) {
		if (vbool_list(i).is_unknown()) m *= 2;
	}
	multiplier.push_back(m);
	
	// note: nFormulas = ineq_list.size()
	int *sol = new int[nFormulas];
	for (unsigned int i = 0; i < nFormulas; i++) {
		if (ineq_list(i).is_true()) sol[i] = 1;
		else if (ineq_list(i).is_false()) sol[i] = 0;
		else sol[i] = -1;
		//std::cout << sol[i] << " ";
	}
	//std::cout << std::endl;
	bsols.push_back(sol);

	//print_model(); std::cout << std::endl;

	//negate the model and add into solver
	z3::expr new_assert = z3context.bool_val(false);
	for (unsigned int i = 0; i < vbool_list.size(); i++) {
		const dagv val = vbool_list(i);
		if (val.is_unknown()) continue;
		else if (val.is_true()) new_assert = (new_assert || !vbool_expr[i]);
		else new_assert = (new_assert || vbool_expr[i]);
	}
	for (unsigned int i = 0; i < ineq_list.size(); i++) {
		const dagv val = ineq_list(i);
		if (val.is_unknown()) continue;
		else if (val.is_true()) new_assert = (new_assert || !ineq_expr[i]);
		else new_assert = (new_assert || ineq_expr[i]);
	}
	z3solver.add(new_assert);
	
	return true;
	
}

const bool volce::solver::get_result() {
	bop_list.init_vals();
	for (unsigned int i = 0; i < bop_list.size(); i++)
		bop_list(i) = eval_op(i);
	for (unsigned int i = 0; i < assert_list.size(); i++)
		if (!get_val(assert_list[i]).is_true()) 
			return false;
	return true;
}

const volce::solver::dagv volce::solver::eval_op(const unsigned int index) {
	// may return unknown even if the formula can be evaluated
	// a: false, b: false, c: unknown
	// (and (or -c a) (or c b)) = false
	// (or -c a) = unknown, (or c b) = unknown, (and unknown unknown) = unknown

	const dagc node(bop_list.type(index), index);
	const std::vector<dagc> c = bop_list[index];
	
	//nodes of propositional operators
	if (node.isand()) {
		//AND
		bool is_unkwn = false;
 		for (unsigned int i = 0; i < c.size(); i++) {
 			dagv val = get_val(c[i]);
 			if (val.is_false()) return dagv(false);
 			else if (val.is_unknown()) is_unkwn = true;
 		}
 		if (is_unkwn) return dagv();
 		else return dagv(true);
	} else if (node.isor()) {
		//OR
		bool is_unkwn = false;
		for (unsigned int i = 0; i < c.size(); i++) {
			dagv val = get_val(c[i]);
			if (val.is_true()) return dagv(true);
			else if (val.is_unknown()) is_unkwn = true;
		}
		if (is_unkwn) return dagv();
		else return dagv(false);
	} else if (node.iscomp()) {
		//EQUAL
		dagv val_l = get_val(c[0]);
		dagv val_r = get_val(c[1]);
		if (val_l.is_unknown() || val_r.is_unknown()) return dagv();
		else return dagv(val_l.is_true() == val_r.is_true());
	} else if (node.isitebool()) {
		//ITE
		dagv val_c = get_val(c[0]);
		dagv val_l = get_val(c[1]);
		dagv val_r = get_val(c[2]);
		if (val_c.is_true()) {
			//(ite true l r) -> l
			return val_l;
		} else if (val_c.is_false()) {
			//(ite false l r) -> r
			return val_r;
		} else if (val_l.is_true() && val_r.is_true()) {
			//(ite unkwn true true) -> true
			return dagv(true);
		} else if (val_l.is_false() && val_r.is_false()) {
			//(ite unkwn false false) -> false
			return dagv(false);
		} else return dagv();
	} else assert(false);

}

const volce::solver::dagv volce::solver::get_val(const dagc node) {

	if (node.iscbool()) {
		//constant: true or false
		return dagv(node.bval());
	} else if (node.isvbool()) {
		//boolean variable
		return (node.isneg()) ? !vbool_list(node.id) : vbool_list(node.id);
	} else if (node.isineq()) {
		//ineq variable
		return (node.isneg()) ? !ineq_list(node.id) : ineq_list(node.id);
	} else if (node.isboolop()) {
		//boolean operator
		return (node.isneg()) ? !bop_list(node.id) : bop_list(node.id);
	} else assert(false);

}

// method for bunch strategy
// input: ineq_list and vbool_list
// output: ineqs and vbools which should be flipped
void volce::solver::get_flip_list(std::vector<bool> &ineq_flip, std::vector<bool> &vbool_flip)
{
	
	//ineq_flip.clear();
	//ineq_flip.resize(ineq_list.size(), true);
	
	//vbool_flip.clear();
	//vbool_flip.resize(vbool_list.size(), true);
	
	//check all previous bunches
	for (unsigned int i = 0; i < bunch_list.size(); i++) {
		
		bunch_elem &prev_bunch = bunch_list[i];
		
		unsigned int differ_count = 0;
		int ineq_id = -1;
		int vbool_id = -1;
		
		// check ineq_val
		// ineq_list.size() == target.ineq_val.size() == prev_bunch.ineq_val.size()
		for (unsigned int j = 0; j < ineq_list.size() && differ_count <= 1; j++) {
		
			if (ineq_list(j).is_unknown() || prev_bunch.ineq_vals[j].is_unknown()) {
				// skip unknown, treat as equal value
				continue;
			} else if (ineq_list(j).is_true() != prev_bunch.ineq_vals[j].is_true()) {
				// differ value found
				if (differ_count == 1) 
					goto GOTO_NEXT_BUNCH;
				else {
					differ_count++;
					ineq_id = j;
				}
			}
		}
		
		// check vbool_val
		// vbool_list.size() == target.vbool_val.size() == prev_bunch.vbool_val.size()
		for (unsigned int j = 0; j < vbool_list.size() && differ_count <= 1; j++) {
		
			if (vbool_list(j).is_unknown() || prev_bunch.vbool_vals[j].is_unknown()) {
				// skip unknown, treat as equal value
				continue;
			} else if (vbool_list(j).is_true() != prev_bunch.vbool_vals[j].is_true()) {
				// differ value found
				if (differ_count == 1) 
					goto GOTO_NEXT_BUNCH;
				else {
					differ_count++;
					vbool_id = j;
				}
			}
		}
		
		assert(differ_count == 1);
		assert(ineq_id >= 0 || vbool_id >= 0);
		assert(ineq_id < 0 || vbool_id < 0);
		
		// a similar bunch found, set dangerous flip
		if (ineq_id >= 0) {
			ineq_flip[ineq_id] = false;
		} else {
			vbool_flip[vbool_id] = false;
		} 
		
GOTO_NEXT_BUNCH:;
		
	}
}
					
					
					
