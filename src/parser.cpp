/*  parser.cpp
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
	get symbol from buffer
*/
std::string volce::solver::get_symbol() {

	char *beg = bufptr;
	
	//first char already scanned
	bufptr++;

	//while not eof	
	while (*bufptr != 0) {
		
		switch (scan_mode) {
		case SM_SYMBOL:
			if (isblank(*bufptr)) {
			
				//out of symbol mode by ' ' and \t
				std::string tmp_s(beg, bufptr - beg);
				
				//skip space
				bufptr++;
				scan_to_next_symbol();
				
				return tmp_s;
				
			} else if(*bufptr == '\n' || *bufptr == '\r' || *bufptr == '\v' || *bufptr == '\f') {
				line_number++;
				
				//out of symbol mode by '\n', '\r', '\v' and '\f'
				std::string tmp_s(beg, bufptr - beg);
				
				//skip space
				bufptr++;
				scan_to_next_symbol();
				
				return tmp_s;
				
			} else if (*bufptr == ';' || *bufptr == '|' || *bufptr == '"' || *bufptr == '(' || *bufptr == ')') {
			
				//out of symbol mode bu ';', '|', '(', and ')'
				std::string tmp_s(beg, bufptr - beg);
				return tmp_s;
				
			}
			break;
			
		case SM_COMP_SYM:
		
			if(*bufptr == '\n' || *bufptr == '\r' || *bufptr == '\v' || *bufptr == '\f') {
				line_number++;
			} else if (*bufptr == '|') {
			
				//out of complicated symbol mode
				bufptr++;
				std::string tmp_s(beg, bufptr - beg);
				
				//skip space
				scan_to_next_symbol();
				
				return tmp_s;
				
			}
			break;
			
		case SM_STRING:
			
			if(*bufptr == '\n' || *bufptr == '\r' || *bufptr == '\v' || *bufptr == '\f') {
				line_number++;
			} else if (*bufptr == '"') {
			
				//out of string mode
				bufptr++;
				std::string tmp_s(beg, bufptr - beg);
				
				//skip space
				scan_to_next_symbol();
				
				return tmp_s;
				
			}
			break;
			
		default:
			assert(false);
		}
		
		//go next char
		bufptr++;
	}
	
	err_unexp_eof();
	
	return NULL;
}

void volce::solver::scan_to_next_symbol() {

	//init scan mode
	scan_mode = SM_COMMON;
	
	//while not eof
	while (*bufptr != 0) {
		
		if (*bufptr == '\n' || *bufptr == '\r' || *bufptr == '\v' || *bufptr == '\f') {
		
			line_number++;
			
			//out of comment mode
			if (scan_mode == SM_COMMENT) scan_mode = SM_COMMON;
			
		} else if (scan_mode == SM_COMMON && !isblank(*bufptr)) {
			
			switch (*bufptr) {
			case ';':
				//encounter comment
				scan_mode = SM_COMMENT;
				break;
			case '|':
				//encounter next complicated symbol
				scan_mode = SM_COMP_SYM;
				return;
			case '"':
				//encounter next string
				scan_mode = SM_STRING;
				return;
			default:
				//encounter next symbol
				scan_mode = SM_SYMBOL;
				return;
			}
			
		}
		
		//go next char
		bufptr++;
	}

}

void volce::solver::parse_lpar() {
	if (*bufptr == '(') {
		bufptr++;
		scan_to_next_symbol();
	} else {
		err_sym_mis("(", line_number);
	}
}

void volce::solver::parse_rpar() {
	if (*bufptr == ')') {
		bufptr++;
		scan_to_next_symbol();
	} else {
		err_sym_mis(")", line_number);
	}
}

void volce::solver::skip_to_rpar() {

	//skip to next right parenthesis with same depth	
	scan_mode = SM_COMMON;
	unsigned int level = 0;
	
	while (*bufptr != 0) {
	
		if (*bufptr == '\n' || *bufptr == '\r' || *bufptr == '\v' || *bufptr == '\f') {
			//new line
			line_number++;
			if (scan_mode == SM_COMMENT) 
				scan_mode = SM_COMMON;
		} else if (scan_mode == SM_COMMON) {
		
			if (*bufptr == '(') level++;
			else if (*bufptr == ')') {
				if (level == 0) return;
				else level--;
			} else if (*bufptr == ';')
				scan_mode = SM_COMMENT; 
			else if (*bufptr == '|')
				scan_mode = SM_COMP_SYM;
			else if (*bufptr == '"')
				scan_mode = SM_STRING;
		
		} else if (scan_mode == SM_COMP_SYM && *bufptr == '|')
			scan_mode = SM_COMMON;
		else if (scan_mode == SM_STRING && *bufptr == '"')
			scan_mode = SM_COMMON;
		
		//go to next char
		bufptr++;
	}

}


void volce::solver::parse_smtlib2_file(const std::string filename) {

	if (solving_initialized) {
		err_solving_initialized();
	}

	/*
		load file
	*/
	std::ifstream fin(filename, std::ifstream::binary);

	if (!fin) {
		err_open_file(filename);
	}

	fin.seekg(0, std::ios::end);
	buflen = (long)fin.tellg();
	fin.seekg(0, std::ios::beg);

	buffer = new char[buflen + 1];
	fin.read(buffer, buflen);
	buffer[buflen] = 0;

	fin.close();


	/*
		parse command
	*/
	bufptr = buffer;
	if (buflen > 0) line_number = 1;

	//skip to first symbol;
	scan_to_next_symbol();
	
	while (*bufptr) {
		parse_lpar();
		if (parse_command() == CT_EXIT) break;
		parse_rpar();
	}

	//parse finished
	key_map.clear();
	delete []buffer;
	
}

volce::CMD_TYPE volce::solver::parse_command() {
	
	unsigned int command_ln = line_number;
	std::string command = get_symbol();

	//(assert <expr>)
	if (command == "assert") {
	
		//parse expression
		dagc assert_expr = parse_expr();
		
		//expression return boolean
		if (!assert_expr.isbool()) {
			err_param_nbool(command, command_ln);
		}
		
		//insert into dag
		assert_list.push_back(assert_expr);
		
		return CT_ASSERT;
		
	}
	
	//(check-sat)
	if (command == "check-sat") {
		skip_to_rpar();
		return CT_CHECK_SAT;
	}
	
	if (command == "check-sat-assuming") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_CHECK_SAT_ASSUMING;
	}

	//(declare-const <symbol> <sort>)
	if (command == "declare-const") {
	
		//get name
		unsigned int name_ln = line_number;
		std::string name = get_symbol();
		
		//get returned type
		dagc res;
		unsigned int type_ln = line_number;
		std::string type = get_symbol();
		if (type == "Bool") {
			res = mk_bool_decl(name);
		} else if (type == "Int") {
			if (logic_not_set()) logic = QF_LIA;
			if (!islia()) err_logic(type, type_ln);
			res = mk_var_decl(name);
		} else if (type == "Real") {
			if (logic_not_set()) logic = QF_LRA;
			if (!islra()) err_logic(type, type_ln);
			res = mk_var_decl(name);
		} else {
			err_unkwn_sym(type, type_ln);
		}
		
		//multiple declarations
		if (res.iserr()) err_all(res, name, name_ln);
	
		return CT_DECLARE_CONST;
		
	}
	
	//(declare <symbol> () <sort>)
	if (command == "declare-fun") {

		//get name
		unsigned int name_ln = line_number;
		std::string name = get_symbol();

		parse_lpar();
		parse_rpar();

		//get returned type
		dagc res;
		unsigned int type_ln = line_number;
		std::string type = get_symbol();
		if (type == "Bool") {
			res = mk_bool_decl(name);
		} else if (type == "Int") {
			if (logic_not_set()) logic = QF_LIA;
			if (!islia()) err_logic(type, type_ln);
			res = mk_var_decl(name);
		} else if (type == "Real") {
			if (logic_not_set()) logic = QF_LRA;
			if (!islra()) err_logic(type, type_ln);
			res = mk_var_decl(name);
		} else {
			err_unkwn_sym(type, type_ln);
		}

		//multiple declarations
		if (res.iserr()) err_all(res, name, name_ln);

		return CT_DECLARE_FUN;

	}

	if (command == "declare-sort") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_DECLARE_SORT;
	}
	
	//(define-fun <symbol> () <sort> <expr>)
	if (command == "define-fun") {

		//get name
		unsigned int name_ln = line_number;
		std::string name = get_symbol();

		parse_lpar();
		parse_rpar();

		//get returned type
		bool isbool = true;
		unsigned int type_ln = line_number;
		std::string type = get_symbol();
		if (type == "Bool") {
			isbool = true;
		} else if (type == "Int") {
			if (logic_not_set()) logic = QF_LIA;
			if (!islia()) err_logic(type, type_ln);
			else isbool = false;
		} else if (type == "Real") {
			if (!islra()) err_logic(type, type_ln);
			else isbool = false;
		} else {
			err_unkwn_sym(type, type_ln);
		}

		//parse expression
		dagc expr = parse_expr();
		
		//expression should return same type with defined func
		if (expr.isbool() != isbool) {
			if (isbool) err_param_nbool(command, command_ln);
			else err_param_nnum(command, command_ln);
		}
		
		//make func
		dagc res = mk_key_bind(name, expr);
		if (res.iserr()) err_all(res, name, name_ln);

		return CT_DEFINE_FUN;
		
	}

	if (command == "define-fun-rec") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_DEFINE_FUN_REC;
	}
	
	if (command == "define-funs-rec") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_DEFINE_FUNS_REC;
	}

	if (command == "define-sort") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_DEFINE_SORT;
	}
	
	if (command == "echo") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_ECHO;
	}

	//(exit)
	if (command == "exit") {
		skip_to_rpar();
		return CT_EXIT;
	}
	
	if (command == "get-assertions") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_ASSERTIONS;
	}
	
	if (command == "get-assignment") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_ASSIGNMENT;
	}
	
	if (command == "get-info") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_INFO;
	}
	
	if (command == "get-option") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_OPTION;
	}
	
	if (command == "get-model") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_MODEL;
	}
	
	if (command == "get-option") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_OPTION;
	}
	
	if (command == "get-proof") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_PROOF;
	}

	if (command =="get-unsat-assumptions") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_UNSAT_ASSUMPTIONS;
	}
	
	if (command == "get-unsat-core") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_UNSAT_CORE;
	}
	
	if (command == "get-value") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_GET_VALUE;
	}
	
	if (command == "pop") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_POP;
	}

	if (command == "push") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_PUSH;
	}
	
	if (command == "reset") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_RESET;
	}

	if (command == "reset-assertions") {
		//ignore
		warn_cmd_nsup(command, command_ln);
		skip_to_rpar();
		return CT_RESET_ASSERTIONS;
	}

	//<attribute ::= <keyword> | <keyword> <attribute_value>
	//(set-info <attribute>)
	if (command == "set-info") {
		skip_to_rpar();
		return CT_SET_INFO;
		
	}
	
	//(set-logic <symbol>)
	if (command == "set-logic") {
		
		unsigned int type_ln = line_number;
		std::string type = get_symbol();
		if (type == "QF_LIA") {
			logic = QF_LIA;
		} else if (type == "QF_LRA") {
			logic = QF_LRA;
		} else {
			err_unkwn_sym(type, type_ln);
		}
		
		return CT_SET_LOGIC;
		
	}
	
	//<option ::= <attribute>
	//(set-option <option>)
	if (command == "set-option") {
		skip_to_rpar();
		return CT_SET_OPTION;
	}

	err_unkwn_sym(command, command_ln);

	return CT_UNKNOWN;
	
}


//expr ::= const | func | (<identifier> <expr>+)
volce::solver::dagc volce::solver::parse_expr() {

	if (*bufptr != '(') {
		//const | func

		unsigned int ln = line_number;
		std::string s = get_symbol();

		if (isdigit(s[0])) {
			//constant number
			return mk_const(s);
		} else {
			boost::unordered_map<std::string, dagc>::iterator kmap_iter = key_map.find(s);
			if (kmap_iter != key_map.end()) {
				//func found
				return kmap_iter->second;
			} else if (s == "true") {
				return mk_true();
			} else if (s == "false") {
				return mk_false();
			} else {
				//unknown symbol
				err_unkwn_sym(s, ln);
			}
		}

	}
	
	//(<identifier> <expr>+)
	parse_lpar();
	
	unsigned int ln = line_number;
	std::string s = get_symbol();
	
	//parse identifier and get params
	dagc expr;
	if (s == "and") expr = parse_and(ln);
	else if (s == "or") expr = parse_or(ln);
	else if (s == "not") expr = parse_not(ln);
	else if (s == "=>") expr = parse_imply(ln);
	else if (s == "xor") expr = parse_xor(ln);
	else if (s == "=") expr = parse_eq(ln);
	else if (s == "distinct") expr = parse_distinct(ln);
	else if (s == "ite") expr = parse_ite(ln);
	else if (s == "+") expr = parse_add(ln);
	else if (s == "-") expr = parse_neg(ln);
	else if (s == "*") expr = parse_mul(ln);
	else if (s == "/") expr = parse_div(ln);
	else if (s == "<=") expr = parse_le(ln);
	else if (s == "<") expr = parse_lt(ln);
	else if (s == ">=") expr = parse_ge(ln);
	else if (s == ">") expr = parse_gt(ln);
	else if (s == "let") expr = parse_let(ln);
	else err_unkwn_sym(s, ln);
	
	parse_rpar();
	
	return expr;

}

/*
	propositional logic
*/

//(and Bool Bool+ :left-assoc), return Bool
volce::solver::dagc volce::solver::parse_and(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_and(params);
	if (res.iserr()) err_all(res, "and", ln);

	return res;
}

//(or Bool Bool+ :left-assoc), return Bool
volce::solver::dagc volce::solver::parse_or(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_or(params);
	if (res.iserr()) err_all(res, "or", ln);

	return res;

}

//(not Bool), return Bool
volce::solver::dagc volce::solver::parse_not(const unsigned int ln) {

	//unary operation
	if (*bufptr == ')') err_param_mis("not", ln);
	dagc expr = parse_expr();
	if (*bufptr != ')') err_param_mis("not", ln);
	if (!expr.isbool()) err_param_nbool("not", ln);
	return mk_not(expr);
}

//(=> Bool Bool+ :right-assoc), return Bool
//(=> a b c d) <=> (or -a -b -c d)
volce::solver::dagc volce::solver::parse_imply(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_imply(params);
	if (res.iserr()) err_all(res, "=>", ln);

	return res;

}

//(xor Bool Bool+ :left-assoc), return Bool
volce::solver::dagc volce::solver::parse_xor(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_xor(params);
	if (res.iserr()) err_all(res, "xor", ln);

	return res;

}

/*
	eq, distinct
*/
	
//(= A A+ :chainable), return Bool
volce::solver::dagc volce::solver::parse_eq(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_eq(params);
	if (res.iserr()) err_all(res, "=", ln);

	return res;

}

//(distinct A A+ :std::pairwise), return Bool
volce::solver::dagc volce::solver::parse_distinct(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_distinct(params);
	if (res.iserr()) err_all(res, "distinct", ln);

	return res;

}

/*
	(ite Bool A A), return A
*/
volce::solver::dagc volce::solver::parse_ite(const unsigned int ln) {

	//parse condition, if-statement and else-statement
	if (*bufptr == ')') err_param_mis("ite", ln);
	dagc condition = parse_expr();
	if (*bufptr == ')') err_param_mis("ite", ln);
	dagc if_statement = parse_expr();
	if (*bufptr == ')') err_param_mis("ite", ln);
	dagc else_statement = parse_expr();
	if (*bufptr != ')') err_param_mis("ite", ln);
	
	dagc res = mk_ite(condition, if_statement, else_statement);
	if (res.iserr()) err_all(res, "ite", ln);

	return res;

}

/*
	arithmetic operators
*/
	
//(+ rt rt+), return rt
volce::solver::dagc volce::solver::parse_add(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_add(params);
	if (res.iserr()) err_all(res, "+", ln);

	return res;

}

//(- rt), return rt
//(- rt rt+), return rt
volce::solver::dagc volce::solver::parse_neg(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res;
	if (params.size() == 1) {
		if (params[0].isbool()) err_param_nnum("-", ln);
		res = mk_neg(params[0]);
	} else {
		res = mk_minus(params);
		if (res.iserr()) err_all(res, "-", ln);
	}

	return res;

}

//(* rt rt+), return rt
volce::solver::dagc volce::solver::parse_mul(const unsigned int ln) {

	std::vector<dagc> params;

	while (*bufptr != ')') params.push_back(parse_expr());
	
	dagc res = mk_mul(params);
	if (res.iserr()) err_all(res, "*", ln);

	return res;

}

//div is an LRA operator
//(/ Real Real), return Real
volce::solver::dagc volce::solver::parse_div(const unsigned int ln) {

	if (*bufptr == ')') err_param_mis("/", ln);
	dagc dividend = parse_expr();
	if (*bufptr == ')') err_param_mis("/", ln);
	dagc divisor = parse_expr();
	if (*bufptr != ')') err_param_mis("/", ln);
	
	dagc res = mk_div(dividend, divisor);
	if (res.iserr()) err_all(res, "/", ln);

	return res;

}

//(<= rt rt), return rt
volce::solver::dagc volce::solver::parse_le(const unsigned int ln) {

	if (*bufptr == ')') err_param_mis("<=", ln);
	dagc lhs = parse_expr();
	if (*bufptr == ')') err_param_mis("<=", ln);
	dagc rhs = parse_expr();
	if (*bufptr != ')') err_param_mis("<=", ln);
	
	dagc res = mk_le(lhs, rhs);
	if (res.iserr()) err_all(res, "<=", ln);

	return res;
	
}

//(< rt rt), return rt
volce::solver::dagc volce::solver::parse_lt(const unsigned int ln) {

	if (*bufptr == ')') err_param_mis("<", ln);
	dagc lhs = parse_expr();
	if (*bufptr == ')') err_param_mis("<", ln);
	dagc rhs = parse_expr();
	if (*bufptr != ')') err_param_mis("<", ln);
	
	dagc res = mk_lt(lhs, rhs);
	if (res.iserr()) err_all(res, "<", ln);

	return res;

}

//(>= rt rt), return rt
volce::solver::dagc volce::solver::parse_ge(const unsigned int ln) {

	if (*bufptr == ')') err_param_mis(">=", ln);
	dagc lhs = parse_expr();
	if (*bufptr == ')') err_param_mis(">=", ln);
	dagc rhs = parse_expr();
	if (*bufptr != ')') err_param_mis(">=", ln);
	
	dagc res = mk_ge(lhs, rhs);
	if (res.iserr()) err_all(res, ">=", ln);

	return res;

}

//(> rt rt), return rt
volce::solver::dagc volce::solver::parse_gt(const unsigned int ln) {

	if (*bufptr == ')') err_param_mis(">", ln);
	dagc lhs = parse_expr();
	if (*bufptr == ')') err_param_mis(">", ln);
	dagc rhs = parse_expr();
	if (*bufptr != ')') err_param_mis(">", ln);
	
	dagc res = mk_gt(lhs, rhs);
	if (res.iserr()) err_all(res, ">", ln);

	return res;

}


/*
	keybinding ::= (<symbol> expr)
	(let (<keybinding>+) expr), return expr
*/
volce::solver::dagc volce::solver::parse_let(const unsigned int ln) {
	
	parse_lpar();
	
	//parse key bindings
	std::vector<std::string> key_list;
	while (*bufptr != ')') {
	
		//(<symbol> expr)
		parse_lpar();	
		
		unsigned int name_ln = line_number;
		std::string name = get_symbol();
		
		dagc res = parse_expr();
		res = mk_key_bind(name, res);
		if (res.iserr()) err_all(res, name, name_ln);
		
		parse_rpar();	
	
		//new key
		key_list.push_back(name);

	}
	
	parse_rpar();
		
	//parse let right expression
	dagc expr = parse_expr();

	//remove key bindings
	while (key_list.size() > 0) {
		key_map.erase(key_list.back());
		key_list.pop_back();
	}

	return expr;

}



