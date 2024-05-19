/*  error.cpp
 *
 *  Copyright (C) 2016-2024 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of sharpSMT.
 *  See COPYING for more information on using this software.
 */



#include <solver.h>

void volce::solver::err_all(const ERROR_TYPE e, const std::string s, const unsigned int ln) const {

	switch (e){
	case ERR_UNEXP_EOF:
		err_unexp_eof();
		break;
	case ERR_SYM_MIS:
		err_sym_mis(s, ln);
		break;
	case ERR_UNKWN_SYM:
		err_unkwn_sym(s, ln);
		break;
	case ERR_PARAM_MIS:
		err_param_mis(s, ln);
		break;
	case ERR_PARAM_NBOOL:
		err_param_nbool(s, ln);
		break;
	case ERR_PARAM_NNUM:
		err_param_nnum(s, ln);
		break;
	case ERR_PARAM_NSAME:
		err_param_nsame(s, ln);
		break;
	case ERR_LOGIC:
		err_logic(s, ln);
		break;
	case ERR_MUL_DECL:
		err_mul_decl(s, ln);
		break;
	case ERR_MUL_DEF:
		err_mul_def(s, ln);
		break;
	case ERR_NLINEAR:
		err_nlinear(s, ln);
		break;
	case ERR_ZERO_DIVISOR:
		err_zero_divisor(ln);
	}

}

void volce::solver::err_all(const dagc e, const std::string s, const unsigned int ln) const {
	err_all((ERROR_TYPE)e.id, s, ln);
}

//unexpected end of file
void volce::solver::err_unexp_eof() const {
	std::cout << "error: Unexpected end of file found." << std::endl;
	exit(0);
}

//symbol missing
void volce::solver::err_sym_mis(const std::string mis, const unsigned int ln) const {
	std::cout << "error: \"" << mis << "\" missing in line " << ln << '.' << std::endl;
	exit(0);
}

void volce::solver::err_sym_mis(const std::string mis, const std::string nm, const unsigned int ln) const {
	std::cout << "error: \"" << mis << "\" missing before \"" << nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//unknown symbol
void volce::solver::err_unkwn_sym(const std::string nm, const unsigned int ln) const {
	if (nm == "") err_unexp_eof();
	std::cout << "error: Unknown or unexptected symbol \"" << nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//wrong number of parameters
void volce::solver::err_param_mis(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Wrong number of parameters of \"" << nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//paramerter type error
void volce::solver::err_param_nbool(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Invalid command \"" << nm << "\" in line " 
			<< ln << ", paramerter is not a boolean." << std::endl;
	exit(0);
}

void volce::solver::err_param_nnum(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Invalid command \"" << nm << "\" in line " 
			<< ln << ", paramerter is not an integer or a real." << std::endl;
	exit(0);
}

//paramerters are not in same type
void volce::solver::err_param_nsame(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Invalid command \"" << nm << "\" in line " 
				<< ln << ", paramerters are not in same type." << std::endl;
	exit(0);
}

//logic doesnt support
void volce::solver::err_logic(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Logic does not support \"" << nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//multiple declaration
void volce::solver::err_mul_decl(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Multiple declarations of \"" << nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//multiple definition
void volce::solver::err_mul_def(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Multiple definitions or keybindings of \"" 
		 << nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//nonlinear arithmetic
void volce::solver::err_nlinear(const std::string nm, const unsigned int ln) const {
	std::cout << "error: Logic does not support nonlinear arithmetic of command \""
				<< nm << "\" in line " << ln << '.' << std::endl;
	exit(0);
}

//divisor is zero
void volce::solver::err_zero_divisor(const unsigned int ln) const {
	std::cout << "error: Divisor is zero in line " << ln << '.' << std::endl;
	exit(0);
}


//global errors
//cannot open file
void volce::solver::err_open_file(const std::string filename) const {
	std::cout << "error: Cannot open file \"" << filename << "\"." << std::endl;
	exit(0);
}

//error while adding new constraints
void volce::solver::err_solving_initialized() const {
	std::cout << "error: Cannot add new constraints, since smt solving has been initialized." << std::endl;
	exit(0);	
}

//unbounded polytope
void volce::solver::err_unbounded_polytope() const {
	std::cout << "\nThe problem is unbounded.\n\n";
	std::cout << "Hint: VolCE provides wordlength parameter (-w) to quickly set bound \n"
			  << "      to each variable with bit-wise domain. For details, check the \n"
			  << "      help menu with '-h' or '--help'.\n";
	exit(1);
}

//logic not support latte
void volce::solver::err_logic_latte() const {
	//std::cout << "warning: Logic does not match LattE." << std::endl;
	//std::cout << "Use '-h' or '--help' for help." << std::endl;
	//exit(0);
}

//logic not support vinci
void volce::solver::err_logic_vinci() const {
	//std::cout << "warning: Logic does not match Vinci." << std::endl;
	//std::cout << "Use '-h' or '--help' for help." << std::endl;
	//exit(0);
}

//logic not support polyvest
void volce::solver::err_logic_polyvest() const {
	//std::cout << "warning: Logic does not match PolyVest." << std::endl;
	//std::cout << "Use '-h' or '--help' for help." << std::endl;
	//exit(0);
}

/*
	warning
*/
//command not support
void volce::solver::warn_cmd_nsup(const std::string nm, const unsigned int ln) const {
	//std::cout << "warning: \"" << nm << "\" command is safely ignored in line " << ln << "." << std::endl;
}

