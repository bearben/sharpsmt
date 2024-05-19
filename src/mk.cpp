/*  mk.cpp
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
	make AND
*/
const volce::solver::dagc volce::solver::mk_and(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);

	std::vector<dagc> new_params;
	
	for (unsigned int i = 0; i < params.size(); i++) {
	
		//check returned type
		if (!params[i].isbool()) return mk_err(ERR_PARAM_NBOOL);
	
		if (params[i].isconst()) {
			if (!params[i].bval())
				//false: one of params is false
				return params[i];
		} else if (params[i].isand() && !params[i].isneg()) {
			//reduce nested AND
			const std::vector<dagc> c = bop_list[params[i].id];
			new_params.insert(new_params.end(), c.begin(), c.end());
		} else
			//insert uncertain param
			new_params.push_back(params[i]);

	}
	
	if (new_params.size() == 0) {
		//all true constant
		return mk_true();
	} else if (new_params.size() == 1) {
		//only one uncertain param
		return new_params[0];
	} else {
		//make new AND operator
		return mk_oper(NT_AND, new_params);
	}
}


/*
	make OR
*/
const volce::solver::dagc	volce::solver::mk_or(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);

	std::vector<dagc> new_params;
	
	for (unsigned int i = 0; i < params.size(); i++) {
	
		//check returned type
		if (!params[i].isbool()) return mk_err(ERR_PARAM_NBOOL);
	
		if (params[i].isconst()) {
			if (params[i].bval())
				//true: one of params is false
				return params[i];
		} else if (params[i].isor() && !params[i].isneg()) {
			//reduce nested OR
			const std::vector<dagc> c = bop_list[params[i].id];
			new_params.insert(new_params.end(), c.begin(), c.end());
		} else
			//insert uncertain param
			new_params.push_back(params[i]);

	}
	
	if (new_params.size() == 0) {
		//all false constant
		return mk_false();
	} else if (new_params.size() == 1) {
		//only one uncertain param
		return new_params[0];
	} else {
		//make new OR operator
		return mk_oper(NT_OR, new_params);
	}
}


/*
	make IMPLY
*/
const volce::solver::dagc volce::solver::mk_imply(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	std::vector<dagc> new_params(params);

	//(=> a b c d) <=> (or -a -b -c d)
	for (unsigned int i = 0; i < new_params.size() - 1; i++)
		new_params[i].negate();
	
	//all syntax checks rely on mk_or func	
	return mk_or(new_params);

}


/*
	make XOR
*/
const volce::solver::dagc volce::solver::mk_xor(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);
	
	//check returned type
	if (!params[0].isbool()) return mk_err(ERR_PARAM_NBOOL);
	
	dagc res = params[0];
	
	//(xor a b c d) <=> (neq (neq (neq a b) c) d)
	for (unsigned int i = 1; i < params.size(); i++) {
	
		//check returned type
		if (!params[i].isbool()) return mk_err(ERR_PARAM_NBOOL);
	
		//(xor pi-1 pi) <=> (neq pi-1 pi)
		res = mk_not(mk_eq_bool(res, params[i]));
		
	}
	
	return res;

}


/*
	make EQUAL
*/
//binary operation version
const volce::solver::dagc volce::solver::mk_eq(const dagc l, const dagc r) {

	if (solving_initialized) err_solving_initialized();

	if (l.isbool() != r.isbool()) return mk_err(ERR_PARAM_NSAME);
	else if (l.isbool()) return mk_eq_bool(l, r);
	else {
		// l = r ==> (l <= r) and (l >= r)
		std::vector<dagc> new_ineqs;
		new_ineqs.push_back(mk_le(l, r));
		new_ineqs.push_back(mk_ge(l, r));
		return mk_and(new_ineqs);
	}

}

const volce::solver::dagc volce::solver::mk_eq(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);
	
	//call binary operation version
	if (params.size() == 2) return mk_eq(params[0], params[1]);
	
	//(= a pivot c d) <=> (and (= pivot a) (= pivot c) (= pivot d))
	std::vector<dagc> new_params;
	unsigned int pivot = 0;
	while (!params[pivot].isconst()) pivot++;
	
	bool isbool = params[pivot].isbool();	
	for (unsigned int i = 0; i < params.size(); i++) {
		if (i == pivot) continue;

		//check returned type
		if (params[i].isbool() != isbool) return mk_err(ERR_PARAM_NSAME);
		
		//(= pivot p[i])
		new_params.push_back(mk_eq(params[pivot], params[i]));
	}
	
	return mk_and(new_params);
	
}

//core of equal, inner use
//assume l and r are boolean-elements
const volce::solver::dagc volce::solver::mk_eq_bool(const dagc l, const dagc r) {

	if (l.isconst()) {
		if (r.isconst()) {
			//left and right are either constants
			if (l.bval() == r.bval()) 
				return mk_true();
			else 
				return mk_false();
		} else if (l.bval()) {
			//(= true r) <=> r
			return r;
		} else {
			//(= false r) <=> (not r)
			return mk_not(r);
		}
	} else if (r.isconst()) {
		if (r.bval()) {
			//(= l true) <=> l
			return l;
		} else {
			//(= l false) <=> (not l)
			return mk_not(l);
		}
	} else {
		//make new EQ operator
		std::vector<dagc> new_params{l, r};
		return mk_oper(NT_EQ, new_params);
	}

}


/*
	make DISTINCT
*/
//binary operation version
//distinct <=> not equal
const volce::solver::dagc volce::solver::mk_distinct(const dagc l, const dagc r) {

	if (solving_initialized) err_solving_initialized();

	if (l.isbool() != r.isbool()) return mk_err(ERR_PARAM_NSAME);
	else if (l.isbool()) return mk_not(mk_eq_bool(l, r));
	else {
		// l != r ==> (l < r) or (l > r)
		std::vector<dagc> new_ineqs;
		new_ineqs.push_back(mk_lt(l, r));
		new_ineqs.push_back(mk_gt(l, r));
		return mk_or(new_ineqs);
	}
}

const volce::solver::dagc volce::solver::mk_distinct(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);
	
	//call binary operation version
	if (params.size() == 2) return mk_distinct(params[0], params[1]);
	
	//pairwise: (distinct a b c) <=> (and (neq a b) (neq a c) (neq b c))
	std::vector<dagc> new_params;
	bool isbool = params.back().isbool();
	
	for (unsigned int i = 0; i < params.size() - 1; i++) {
		//check returned type
		if (params[i].isbool() != isbool) return mk_err(ERR_PARAM_NSAME);
		
		for (unsigned int j = i + 1; j < params.size(); j++)
			new_params.push_back(mk_distinct(params[i], params[j]));
	}
	
	return mk_and(new_params);

}


/*
	make ITE
*/
const volce::solver::dagc volce::solver::mk_ite(const dagc c, const dagc l, const dagc r) {

	if (solving_initialized) err_solving_initialized();

	bool isbool = l.isbool();

	//check condition type
	if (!c.isbool()) return mk_err(ERR_PARAM_NBOOL);
	
	//check returned type
	if (isbool != r.isbool()) return mk_err(ERR_PARAM_NSAME);
	
	//condition is constant
	if (c.isconst()) {
		if (c.bval()) return l;
		else return r;
	}
	
	//make new ITE operator
	std::vector<dagc> new_params{c, l, r};
	if (isbool) 
		return mk_oper(NT_ITE_BOOL, new_params);
	else 
		return mk_num_oper(NT_ITE_NUM, new_params);

}


/*
	make ADD
*/
const volce::solver::dagc volce::solver::mk_add(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);

	std::vector<dagc> new_params;
	double value = 0;
	
	for (unsigned int i = 0; i < params.size(); i++) {
		
		//check returned type
		if (params[i].isbool()) return mk_err(ERR_PARAM_NNUM);
		
		if (params[i].isconst()) {
			//merge constants
			value += params[i].nval();
		} else if (params[i].isadd()) {
			//reduce nested addition
			//(+ a (+ b c) d) <=> (+ a b c d)
			value += params[i].nval();
			//assert(params[i].m == 1);
			//value += params[i].v;
			const std::vector<dagc> c = nop_list[params[i].id];
			new_params.insert(new_params.end(), c.begin(), c.end());
		} else
			new_params.push_back(params[i]);
		
	}
	
	if (new_params.size() == 0) {
		//all constants
		return mk_const(value);
	} else if (new_params.size() == 1 && value == 0) {
		return new_params[0];
	} else {
		//push the costant value at last
		return mk_num_oper(NT_ADD, new_params, value);
	}


}


/*
	make MINUS
*/
const volce::solver::dagc volce::solver::mk_minus(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	std::vector<dagc> new_params(params);

	//(- a b c d) <=> (+ a -b -c -d)
	for (unsigned int i = 1; i < new_params.size(); i++)
		new_params[i].negate();
	
	//all syntax checks rely on mk_add func	
	return mk_add(new_params);

}


/*
	make MUL
*/
const volce::solver::dagc volce::solver::mk_mul(const std::vector<dagc> &params) {

	if (solving_initialized) err_solving_initialized();

	//check the number of params
	if (params.size() < 2) return mk_err(ERR_PARAM_MIS);

	std::vector<dagc> new_params;
	double value = 1;
	
	for (unsigned int i = 0; i < params.size(); i++) {

		//check returned type
		if (params[i].isbool()) return mk_err(ERR_PARAM_NNUM);
		
		if (params[i].isconst()) {
			//merge constants
			value *= params[i].nval();
			if (value == 0) return mk_const(0);
		} else
			new_params.push_back(params[i]);

	}
	
	if (new_params.size() == 0) {
		//all constants
		return mk_const(value);
	} else if (new_params.size() == 1) {
		return mk_mul(new_params[0], value);		
	} else {
		//make new MUL operator
		return mk_num_oper(NT_MUL, new_params, 0, value);
	}
}

//(* param multiplier)
const volce::solver::dagc volce::solver::mk_mul(const dagc p, const double m) {

	if (solving_initialized) err_solving_initialized();

	if (m == 1) {
		// 1 * p = p
		return p;
	} else if (m == 0) {
		// 0 * p = 0
		return mk_const(0);
	} else if (p.isadd()) {
		// '+' expr, distributive law
		for (unsigned int i = 0; i < nop_list[p.id].size(); i++) 
			nop_list[p.id][i].m *= m;
		return dagc(p.t, p.id, p.v * m);
	} else
		return dagc(p.t, p.id, p.v, p.m * m);

}


/*
	make DIV
*/
const volce::solver::dagc volce::solver::mk_div(const dagc l, const dagc r) {

	if (solving_initialized) err_solving_initialized();

	//check logic, enabled for LRA
	//if (!islra()) return mk_err(ERR_LOGIC);

	//check returned type
	if (l.isbool() || r.isbool()) return mk_err(ERR_PARAM_NNUM);

	if (r.isconst()) {
		if (l.isconst()) {
			//all constants
			return mk_const(l.nval() / r.nval());
		} else {
			//divisor is zero
			std::vector<dagc> new_params{l, mk_const(1.0 / r.nval())};
			return mk_mul(new_params);
		}
	} else {
		//make new DIV operator
		std::vector<dagc> new_params{l, r};
		return mk_num_oper(NT_DIV, new_params);
	}

}


/*
	make INEQ
*/
const volce::solver::dagc volce::solver::mk_ineq(const bool iseq, const dagc l, const dagc r) {

	//check returned type
	if (l.isbool() || r.isbool()) return mk_err(ERR_PARAM_NNUM);

	dagc res = mk_ineq_core(iseq, l, r);
	
	if (res.iserr()) {

		//exist ITEs
		//enumerate paths with different assignments of conditions
		//(and path1 path2 ...)
		std::vector<dagc> aparams;
		
		unsigned int size = nop_list.size();
		
		while (true) {
			
			//provide assignments for conditions in ITEs and remake AST
			dagc lhs = remk_ineq(l);
			dagc rhs = remk_ineq(r);
		
			res = mk_ineq_core(iseq, lhs, rhs);
			if (res.iserr()) return res;
			
			//(or -cond1 ... -condn ineq)
			std::vector<dagc> oparams(cond_stack);
			for (unsigned int i = 0; i < cond_stack.size(); i++)
				oparams[i].negate();
			oparams.push_back(res);
			aparams.push_back(mk_or(oparams));
			
			//pop conditions with negative assignments
			while (cond_stack.size() > 0 && !cond_stack.back().bval()) 
				cond_stack.pop_back();
			
			//negate the last condition or break out
			if (cond_stack.size() > 0) cond_stack.back().negate();
			else break;
			
		}
		
		nop_list.resize(size);
		
		return mk_and(aparams);
	
	} else return res;
	
}

const volce::solver::dagc volce::solver::mk_ineq_core(const bool iseq, const dagc lhs, const dagc rhs) {

	double constant = 0;
	std::vector<term> terms;

	if (lhs.iscnum()) 
		//<const>
		constant += lhs.nval();
	else if (lhs.isvnum()) 
		//<term>
		terms.push_back(term(lhs.id, lhs.m));
	else if (lhs.isadd()) {
		//(+ <const> <term>+)
		constant += lhs.nval();
		const std::vector<dagc> c = nop_list[lhs.id];
		for (unsigned int i = 0; i < c.size(); i++)
			if (c[i].isvnum()) {
				terms.push_back(term(c[i].id, c[i].m));
			} else return mk_err(ERR_NLINEAR);
	} else return mk_err(ERR_NLINEAR);

	if (rhs.iscnum()) 
		//<const>
		constant -= rhs.nval();
	else if (rhs.isvnum()) 
		//<term>
		terms.push_back(term(rhs.id, -rhs.m));
	else if (rhs.isadd()) {
		//(+ <const> <term>+)
		constant -= rhs.nval();
		const std::vector<dagc> c = nop_list[rhs.id];
		for (unsigned int i = 0; i < c.size(); i++)
			if (c[i].isvnum()) {
				terms.push_back(term(c[i].id, -c[i].m));
			} else return mk_err(ERR_NLINEAR);
	} else return mk_err(ERR_NLINEAR);
	
	//only have constants
	if (terms.size() == 0) {
		if (iseq) {
			if (constant == 0) return mk_true();
			else return mk_false();
		} else {
			if (constant <= 0) return mk_true();
			else return mk_false();
		}
	}
	
	//make inequality
	ineqc ie(iseq, constant, terms);
	return dagc(NT_INEQ, ineq_list.push_back(ie));

}

/*
	remake INEQ
*/
const volce::solver::dagc volce::solver::remk_ineq(const dagc root) {

	if (root.iscnum() || root.isvnum()) {
		//constant or numeric variable
		return root;
	}
	
	const std::vector<dagc> c = nop_list[root.id];
	
	if (root.isitenum()) {
		std::vector<dagc>::iterator it = find(cond_stack.begin(), cond_stack.end(), c[0]);
		bool cond = true;
		if (it == cond_stack.end()) {
			//new condition
			cond_stack.push_back(c[0]);
			if (!(cond = c[0].bval())) 
				cond_stack.back().negate();
		} else {
			//existed
			cond = (c[0].bval() == it->bval());
		}
		if (cond) return mk_mul(remk_ineq(c[1]), root.m);
		else return mk_mul(remk_ineq(c[2]), root.m);
		
	} else if (root.isadd()) {
		//remake add
		std::vector<dagc> new_params{mk_const(root.v)};
		for (unsigned int i = 0; i < c.size(); i++)
			new_params.push_back(remk_ineq(c[i]));
		return mk_mul(mk_add(new_params), root.m);
		
	} else if (root.ismul()) {
		//remake mul
		std::vector<dagc> new_params{mk_const(root.m)};
		for (unsigned int i = 0; i < c.size(); i++)
			new_params.push_back(remk_ineq(c[i]));
		return mk_mul(new_params);
		
	} else if (root.isdiv()) {
		//remake div
		return mk_mul(mk_div(remk_ineq(c[0]), remk_ineq(c[1])), root.m);
		
	} else assert(false);

}


/*
	make declared bool
*/
const volce::solver::dagc volce::solver::mk_bool_decl(const std::string &name) {

	if (solving_initialized) err_solving_initialized();

	dagc expr(NT_VBOOL, vbool_list.size());
	
	std::pair<boost::unordered_map<std::string, dagc>::iterator, bool>
			p = key_map.insert(std::pair<std::string, dagc>(name, expr));
			
	if (!p.second) {
		//multiple declarations
		return mk_err(ERR_MUL_DECL);
	} else {
		//insertion took place
		if (!vbool_list.push_back(name)) assert(false);
		return expr;
	}
}


/*
	make declared var
*/
const volce::solver::dagc volce::solver::mk_var_decl(const std::string &name) {

	if (solving_initialized) err_solving_initialized();

	dagc expr(NT_VNUM, vnum_list.size());

	std::pair<boost::unordered_map<std::string, dagc>::iterator, bool>
			p = key_map.insert(std::pair<std::string, dagc>(name, expr));
			
	if (!p.second) {
		//multiple declarations
		return mk_err(ERR_MUL_DECL);
	} else {
		//insertion took place
		if (!vnum_list.push_back(name)) assert(false);
		return expr;
	}
}


/*
	make key binding
*/
const volce::solver::dagc volce::solver::mk_key_bind(const std::string &key, const dagc expr) {

	if (solving_initialized) err_solving_initialized();

	std::pair<boost::unordered_map<std::string, dagc>::iterator, bool>
			p = key_map.insert(std::pair<std::string, dagc>(key, expr));
			
	if (!p.second) {
		//multiple key bindings
		return mk_err(ERR_MUL_DECL);
	} else {
		//insertion took place
		return expr;
	}
}



