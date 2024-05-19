/*  ineq.cpp
 *
 *  Copyright (C) 2016-2024 Cunjing Ge.
 *
 *  All rights reserved.
 *
 *  This file is part of sharpSMT.
 *  See COPYING for more information on using this software.
 */


#include <ineq.h>

const std::vector<double> volce::ineqc::get_key() const {

	std::vector<double> key{(double)eq, cst};
	for (unsigned int i = 0; i < tm.size(); i++){
		key.push_back(tm[i].id);
		key.push_back(tm[i].m);
	}
	return key;

}

void volce::ineqc::merge_terms(){
	
	sort(tm.begin(), tm.end());
	std::vector<term>::iterator p, q;
	for (p = tm.begin(), q = tm.begin() + 1; q != tm.end(); q++)
		if (p->id == q->id){
			// combine term with same index
			p->m += q->m;
		}else{
			// go next if mult is not zero
			// eleminate terms with zero multiplier
			if (p->m != 0) p++;
			*p = *q;
		}
	tm.resize(p - tm.begin() + 1);		
}

void volce::ineqc::unify(){

	merge_terms();
	
	// make first term postive for equations
	if (!eq) return;
	else if (tm.size() == 0) cst = abs(cst);
	else if (tm[0].m < 0) reverse_mults();
	
}

void volce::ineqc::reverse_mults(){

	cst = -cst;
	
	for (unsigned int i = 0; i < tm.size(); i++)
		tm[i].m = -tm[i].m;

}




