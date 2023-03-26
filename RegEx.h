#ifndef CPP_CODE_REGEX_H
#define CPP_CODE_REGEX_H

#include "template.h"
using namespace std;

class RegEx {
public:
	const static int ZERO = -1, ONE = -2, PLUS = -3, DOT = -4, STAR = -5;
	RegEx* L;
	RegEx* R;
	int eId;
	void copy(RegEx* other){
		L = other -> L;
		R = other -> R;
		eId = other -> eId;
	}
	RegEx() {}
	RegEx(int _eId) : eId(_eId), L(nullptr), R(nullptr) {}
	RegEx(int type, RegEx* _L, RegEx* _R){
		eId = type;
		L = _L;
		R = _R;
		if ((L == NULL or R == NULL) and eId != STAR)
			return;
		if (eId == PLUS){
			if (L -> to_string() == R -> to_string()){
				copy(L);
			}
			else if (L -> eId == ZERO)	
				copy(R);
			else if (R -> eId == ZERO)
				copy(L);
		}
		else if (eId == DOT){
			if (L -> eId == ZERO or R -> eId == ZERO)
				eId = ZERO, L = R = NULL;
			else if(L -> eId == ONE)
				copy(R);
			else if(R -> eId == ONE)
				copy(L);
		}
		else if (eId == STAR){
			if (L -> eId == ONE or L -> eId == ZERO)
				eId = ONE, L = R = NULL;
		}
	}
	// to create a regular expression 0, RegEx(RegEx::ZERO)
	// to create a regular expression 1, RegEx(RegEx::ONE)
	// to create a regular expression (re)* [where re is another RegEx], RegEx(RegEx::STAR, re, nullptr)
	// to create a regular expression re1 + re2 [where re1, re2 are another RegExs], RegEx(RegEx::PLUS, re1, re2)
	// to create a regular expression re1 . re2 [where re1, re2 are another RegExs], RegEx(RegEx::DOT, re1, re2)
	friend auto operator<<(ostream& os, RegEx const& m) -> ostream& {
		if (m.eId == ZERO)
			os << "0";
		else if (m.eId == ONE)
			os << "1";
		else if (m.eId == PLUS){
			if (m.L == NULL and m.R != NULL)
				os << (*(m.R));
			else if (m.L != NULL and m.R == NULL)
				os << (*(m.L));
			else 
				os << "(" << (*(m.L)) << " + " << (*(m.R)) << ")";
		}
			
		else if (m.eId == DOT){
			if (m.L == NULL and m.R != NULL)
				os << (*(m.R));
			else if (m.L != NULL and m.R == NULL)
				os << (*(m.L));
			else{
				assert(m.L -> eId != ONE);
				os << "(" << (*(m.L)) << " . " << (*(m.R)) << ")";
			}
		}
		else if (m.eId == STAR){
			if (m.L != NULL){
				if ((m.L -> eId == PLUS or m.L -> eId == DOT) 
					and m.L -> L != NULL and m.L -> R != NULL)
					os << (*(m.L)) << "*";
				else
					os << "(" << (*(m.L)) << ")*";
			}
		}
		else
			os << "e" << m.eId;

		return os;
	}

	string eId_to_string(){
		int x = eId;
		assert (x >= 0);
		if (x == 0)
			return "e0";
		string res;
		while (x > 0){
			res += (x%10) + '0';
			x /= 10;
		}
		reverse(all(res));
		return "e" + res;
	}
	string to_string(){
		if (eId == ZERO)
			return "0";
		else if (eId == ONE)
			return "1";
		else if (eId == PLUS){
			if (L == NULL and R != NULL)
				return R -> to_string();
			else if (L != NULL and R == NULL)
				return L -> to_string();
			else 
				return "(" + L -> to_string() + " + " + R -> to_string() + ")";
		}
			
		else if (eId == DOT){
			if (L == NULL and R != NULL)
				return R -> to_string();
			else if (L != NULL and R == NULL)
				return L -> to_string();
			else{
				assert(L -> eId != ONE);
				return  "(" + L -> to_string() + " . " + R -> to_string() + ")";
			}
		}
		else if (eId == STAR){
			if (L != NULL){
				if ((L -> eId == PLUS or L -> eId == DOT) 
					and L -> L != NULL and L -> R != NULL)
					return L -> to_string() + "*";
				else
					return "(" + L -> to_string() + ")*";
			}
		}
		else
			return eId_to_string();
		return "";
	}
};


#endif //CPP_CODE_REGEX_H