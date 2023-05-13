/**
* @mainpage
*
* Space-efficient sparse variant of an RNA (loop-based) free energy
* minimization algorithm (RNA folding equivalent to the Zuker
* algorithm).
*
* The results are equivalent to RNAfold -d0.
*
* Demonstration of space-efficient sparsification with trace-back.
*
* Since many matrix entries can not be efficiently recomputed in
* trace back, we store trace arrows to such entries. To save space,
* trace arrows are gc'ed and trace arrows to candidates are omitted
* and reconstructed in trace back.
*
* ----------------------------------------
* Specific recursions:

W(i,j) = min { W(i,j-1);
				min_i<k<j  W(i,k-1) + V(k,j) <-- W (same i), CLW;
				V(i,j);
		0 if i>=j-m
		}

V(i,j) = min { HairpinE(i,j);
		min_kl V(i,j)+ILoopE(i,j,k,l) <-- TAs;
		WM2(i+1,j-1) + a <-- WM2, no TAs;
		}

WM(i,j) = min { V(i,j)+b      <-- candidate in recomp;
		WM(i,j-1) + c <-- ! not via candidate list;
				min_i<k<j (k-i)*c + V(k,j) + b <-- CLWM   ( trick to save trace arrows );
				min_i<k<j  WM(i,k-1) + V(k,j) + b  <-- WM, CLWM;
		INF if i>=j-m
		}

WM2(i,j) = min{ WM2(i,j-1) + c;
				min_i<k<j  WM(i,k-1) + V(k,j) + b }  <-- WM, CLWM, no TAs;

* ----------------------------------------
* Candidate criteria:
*
(i,j) is a candidate for the split in W if
V(i,j)      < min {
					W(i,j-1);
			min_i<k<j  W(i,k-1) + V(k,j)
					}

(i,j) is a candidate for the split in WM if
V(i,j) + b  < min {
					WM(i,j-1)+c;
			min_i<k<j (k-i)*c + V(k,j) + b;
			min_i<k<j  WM(i,k-1) + V(k,j) + b
					}

*
* For simplicity and space savings, we store all candidates that
* meet either criterion in the same list.
*/

#include <iostream>
#include <iomanip>

#include <sstream>

#include <matrix.hh>

#include <limits>

#include <vector>
#include <iterator>

#include <cstring>
#include <string>
#include <cassert>
#include <numeric>

#include "base.hh"
#include "trace_arrows.hh"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}

#include "cmdline.hh"
#include <omp.h>

#include <fstream>


typedef unsigned short int cand_pos_t;


struct triplet
{
    cand_pos_t first; 
    energy_t second;
    energy_t third;
	triplet(){
		first = 1;
		second = 2;
		third = 3;
	}
	triplet(cand_pos_t x, energy_t y , energy_t z){
		first = x;
		second = y;
		third = z;
	}
};


// typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
// typedef std::vector< cand_entry_t > cand_list_t;

typedef triplet cand_entry_td1;
typedef std::vector< cand_entry_td1 > cand_list_td1;

class SparseMFEFold;

energy_t ILoopE(auto const& S_,auto const& S1_, auto const& params_,const int& ptype_closing, const size_t& i, const size_t& j, const size_t& k,  const size_t& l);
energy_t MbLoopE(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
energy_t Mlstem(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,int32_t* p_table,int32_t *last_j_array, int32_t *in_pair_array);
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,int32_t* p_table,int32_t *last_j_array, int32_t *in_pair_array);
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,int32_t* p_table,int32_t *last_j_array, int32_t *in_pair_array) ;
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,int32_t* p_table,int32_t *last_j_array, int32_t *in_pair_array);

bool evaluate_restriction(int i, int j, int32_t *last_j_array, int32_t *in_pair_array);

/**
* Space efficient sparsification of Zuker-type RNA folding with
* trace-back. Provides methods for the evaluation of dynamic
* programming recursions and the trace-back.
*/
class SparseMFEFold {

public:
	std::string seq_;
	size_t n_;

	short *S_;
	short *S1_;

	paramT *params_;

	std::string structure_;
	std::string restricted_;
	

	bool garbage_collect_;

	LocARNA::Matrix<energy_t> V_; // store V[i..i+MAXLOOP-1][1..n]
	std::vector<energy_t> W_;
	std::vector<energy_t> WM_;
	std::vector<energy_t> WM2_;

	std::vector<energy_t> dmli1_; // WM2 from 1 iteration ago
	std::vector<energy_t> dmli2_; // WM2 from 2 iterations ago

	bool mark_candidates_;


	TraceArrows ta_;

	// TraceArrows ta_dangle_;
	
	std::vector< cand_list_td1 > CL_;

	
	

	/**
	candidate list for decomposition in W or WM

	@note Avoid separate candidate lists CLW and CLWM for split cases in W and
	WM to save even more space; here, this works after
	reformulating the recursions such that both split-cases recurse to
	V-entries. (compare OCTs)
	*/
	

	// compare candidate list entries by keys (left index i) in descending order
	struct {
	bool operator ()(const cand_entry_td1 &x, size_t y) const {
		return x.first > y;
	}
	}
	cand_comp;

	


	SparseMFEFold(const std::string &seq, bool garbage_collect, std::string restricted)
	: seq_(seq),
	n_(seq.length()),
	params_(scale_parameters()),
	ta_(n_),
	// ta_dangle_(n_),
		garbage_collect_(garbage_collect)
	{
	make_pair_matrix();

	S_ = encode_sequence(seq.c_str(),0);
	S1_ = encode_sequence(seq.c_str(),1);

	V_.resize(MAXLOOP+1,n_+1);
	W_.resize(n_+1,0);

	WM_.resize(n_+1,INF);

	WM2_.resize(n_+1,INF);

	dmli1_.resize(n_+1,INF);

	dmli2_.resize(n_+1,INF);

	// init candidate lists
	CL_.resize(n_+1);

	resize(ta_,n_+1);

	// resize(ta_dangle_,n_+1);

	restricted_ = restricted;
	
	}

	

	~SparseMFEFold() {
	free(params_);
	free(S_);
	free(S1_);
	}
};


// ! TRANSLATED: -----------------------------------------------------------------------------------

energy_t HairpinE(auto const& seq, auto const& S, auto const& S1, auto const& params, size_t i, size_t j) {

	assert(1<=i);
	assert(i<j);
	
	//assert(j<=len); // don't know len here

	const int ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
}

energy_t encode(energy_t v, int d){
    return (v << 2) | d;
}

std::pair<energy_t,int> decode(energy_t v){
    return std::make_pair((v >> 2), (v & 3));
}


energy_t E_ext_Stem(auto const& vkj,auto const& vk1j,auto const& vkj1,auto const& vk1j1,auto const& S, auto const& params, const size_t i,const size_t j, int &d, size_t n, auto const& p_table){

	int e = INF;
	int en = INF;
  	unsigned int tt  = pair[S[i]][S[j]];

    if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j && p_table[j] == i)) {
				en = vkj; // i j

                int si1 = i>1 ? S[i-1] : -1;
                int sj1 = j<n ? S[j+1] : -1;
				if (en != INF) {
					if (params->model_details.dangles == 2)
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
                    else
                        en += vrna_E_ext_stem(tt, -1, -1, params);

                    e = MIN2(e, en);
					
				}

	}

	if(params->model_details.dangles == 1){
        tt  = pair[S[i+1]][S[j]];
        if (((p_table[i+1] <-1 && p_table[j] <-1) || (p_table[i+1] == j)) && p_table[i]<0) {
            en = (j-i-1>TURN) ? vk1j : INF; //i+1 j

            if (en != INF) {

                int si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = MIN2(e,en);
            if(e == en){
                d=1;
            }

        }
        tt  = pair[S[i]][S[j-1]];
        if (((p_table[i] <-1 && p_table[j-1] <-1) || (p_table[i] == j-1)) && p_table[j]<0) {
            en = (j-1-i>TURN) ? vkj1 : INF; // i j-1
            if (en != INF) {

                int sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = MIN2(e,en);
            if(e == en){
                d=2;
            }

        }
        tt  = pair[S[i+1]][S[j-1]];
        if (((p_table[i+1] <-1 && p_table[j-1] <-1) || (p_table[i+1] == j-1)) && p_table[i] < 0 && p_table[j]<0) {
            en = (j-1-i-1>TURN) ? vk1j1 : INF; // i+1 j-1

            if (en != INF) {

                int si1 = S[i];
                int sj1 = S[j];

                en += vrna_E_ext_stem(tt, si1, sj1, params);
            }
            e = MIN2(e,en);
            if(e == en){
                d=3;
            }
        }
	}
	return e;
}


/**
* @brief Rotate WM2 arrays to store the previous and previous previous iterations
* @param WM2 WM2 array
* @param dmli1 WM2 from one iteration ago
* @param dmli2 WM2 from two iterations ago
*/
void rotate_arrays(auto &WM2, auto &dmli1, auto &dmli2){
	dmli2.swap(dmli1);
    dmli1.swap(WM2);
}


/**
* @brief Computes the multiloop V contribution (in essence VM)
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago
* @param S Sequence Encoding
* @param params Parameters
* @param i Current i
* @param j Current j
* @param p_table Restricted Array
* @return energy_t 
*/
energy_t E_MbLoop(auto const& dmli1, auto const& dmli2, auto const& S, auto const& params, size_t i, size_t j, auto const& p_table){

	int e = INF;
	int en = INF;
  	unsigned int tt  = pair[S[j]][S[i]];
	bool pairable = (p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j);
	
	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			if (pairable) {
			e = dmli1[j - 1];

			if (e != INF) {

				int si1 = S[i + 1];
				int sj1 = S[j - 1];

				e += E_MLstem(tt, sj1, si1, params) + params->MLclosing;
			}

			}
			break;

		case 1:
			/**
			* ML pair D0
			*  new closing pair (i,j) with mb part [i+1,j-1]  
			*/
			
			if (pairable) {
        		e = dmli1[j - 1];

        		if (e != INF) {

          			e += E_MLstem(tt, -1, -1, params) + params->MLclosing;

        		}
      		}
      		/** 
			* ML pair 5
			* new closing pair (i,j) with mb part [i+2,j-1] 
			*/
      		if (pairable && p_table[i+1] < 0) {
        		en = dmli2[j - 1];

        		if (en != INF) {

          			int si1 =  S[i + 1];

          			en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;
      
        		}
      		}
      		e   = MIN2(e, en);

			/** 
			* ML pair 3
			* new closing pair (i,j) with mb part [i+1, j-2] 
			*/
			if (pairable && p_table[j-1] < 0) {
				en = dmli1[j - 2];

				if (en != INF) {
					int sj1 = S[j - 1];

					en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase; 
				}
			}
			e   = MIN2(e, en);
			/** 
			* ML pair 53
			* new closing pair (i,j) with mb part [i+2.j-2]
			*/
			if (pairable && p_table[i+1] < 0 && p_table[j-1] <0) {
				en = dmli2[j - 2];

				if (en != INF) {

					int si1 = S[i + 1];
					int sj1 = S[j - 1];

					en += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
				}
			}
			e   = MIN2(e, en);
      		break;
		case 0:
			if (pairable) {
				e = dmli1[j - 1];

				if (e != INF) {
					e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
				}
			}
			break; 
	}


	return e;
}
/**
* @brief Computes the Multiloop WM contribution 
* 
* @param vkj V at k and j
* @param vk1j V at k+1 and j
* @param vkj1 V at k and j-1
* @param vk1j1 V at k+1 and j-1
* @param WM WM array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params Parameters
* @param i Current i
* @param j Current j
* @param n Length
* @param p_table Restricted array
* @return energy_t 
*/
energy_t E_MLStem(auto const& vkj,auto const& vk1j,auto const& vkj1,auto const& vk1j1, auto const& WM, auto const& CL,auto const& S, auto const& params,size_t i, size_t j,int &d, auto const& n, auto const& p_table){

	int e = INF,en=INF;

	int type = pair[S[i]][S[j]];


	if ((p_table[i] < -1 && p_table[j] < -1) || (p_table[i] == j)) {
		en = vkj; // i j
		if (en != INF) {
            int si1 = i>1 ? S[i-1] : -1;
            int sj1 = j<n ? S[j+1] : -1;
			if (params->model_details.dangles == 2)
				en += E_MLstem(type, si1, sj1, params);
			else
				en += E_MLstem(type, -1, -1, params);

			e = MIN2(e, en);
		}
	}
	if(params->model_details.dangles == 1){
		int mm5 = S[i], mm3 = S[j];

		if (((p_table[i+1] < -1 && p_table[j] < -1) || (p_table[i+1] == j)) && p_table[i] < 0) {
      		en = (j-i-1 >TURN) ? vk1j : INF; // i+1 j
      		if (en != INF) {
        		en += params->MLbase;

            	type = pair[S[i+1]][S[j]];
            	en += E_MLstem(type, mm5, -1, params);

        		e = MIN2(e, en);
				if(e == en){
					d=1;
				}
      		}
    	}

		if (((p_table[i] < -1 && p_table[j-1] < -1) || (p_table[i] == j-1)) && p_table[j] < 0) {
      		en = (j-1-i>TURN) ? vkj1 : INF; // i j-1
      		if (en != INF) {
       			en += params->MLbase;

            	type = pair[S[i]][S[j-1]];
            	en += E_MLstem(type, -1, mm3, params);
 
        		e = MIN2(e, en);
				if(e == en){
					d=2;
				}
      		}
    	}
    	if (((p_table[i+1] < -1 && p_table[j-1] < -1) || (p_table[i+1] == j-1)) && p_table[i] < 0 && p_table[j]<0) {
      		en = (j-1-i-1>TURN) ? vk1j1 : INF; // i+1 j-1
      		if (en != INF) {
        		en += 2 * params->MLbase;

        		type = pair[S[i+1]][S[j-1]];
        		en += E_MLstem(type, mm5, mm3, params);
        
				e = MIN2(e, en);
				if(e == en){
					d=3;
				}
      		}
    	} 
		
	}


    return e;
}



/**
* @brief Recompute row of W
* 
* @param W W array
* @param CL Candidate List
* @param i Current i
* @param max_j Current j
* @param p_table Restricted array
* @return auto const 
*/
auto const recompute_W(auto const &W, auto const& CL, size_t i, size_t max_j, int32_t* p_table) {
	
	std::vector<energy_t> temp = W;
	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=0; }

	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {

		energy_t w = INF;

		// note: the loop covers the case W(i,j)=V(i,j),
		// since this is in the candidate list (TBS)
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			w = std::min( w, temp[it->first-1] + it->second );
		}
		// case "j unpaired" is not in the CL (anymore)
		if(p_table[j-1]<-1) w = std::min(w,temp[j-1]);

		temp[j] = w;
	}
	return temp;
}



/**
* @brief Recompute row of WM 
* 
* @param WM WM array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params Parameters
* @param n length
* @param i Current i
* @param max_j Current j
* @param p_table Restricted array
* @return auto const 
*/
auto const recompute_WM(auto const& WM, auto const &CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, int32_t* p_table, int32_t* up_array) {
	

	assert(i>=1);
	assert(max_j<=n);

	std::vector<energy_t> temp = WM;

	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }
	
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;
		bool paired = 0;
		int ptype = 0;		
		// #pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
			paired += (p_table[k] == j);
			auto const [vd,d] = decode(it->third);
			energy_t vkjwm = vd + params->MLintern[2];
			bool can_pair = up_array[k-1] >= (k-i) ? true : false;
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + vkjwm );
			wm = std::min( wm, temp[k-1]  + vkjwm );
			// if(j==354) printf("k is %lu and j is %lu and WM[k-1] is %d and vkj is %d and wm is %d\n",k,j,WM[k-1],v_kj,wm);
		}
		if(p_table[j]<0 && !paired) wm = std::min(wm, temp[j-1] + params->MLbase);
		temp[j] = wm;
	}
	return temp;
}

/**
* @brief Recompute row of WM2 
* 
* @param WM WM array
* @param WM2 WM2 array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params parameters
* @param n length
* @param i current i
* @param max_j current j
* @param p_table restricted array
* @param last_j_array restricted array
* @param in_pair_array restricted array
* @return auto const 
*/
auto const recompute_WM2(auto const& WM, auto const& WM2, auto const CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, int32_t* p_table) {
	

	assert(i>=1);
	//assert(i+2*TURN+3<=max_j);
	assert(max_j<= n);

	std::vector<energy_t> temp = WM2;

	for ( size_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	// #pragma omp parallel for num_threads(6);
	for ( size_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;
		bool paired = 0;
		int ptype = 0;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			size_t k = it->first;
			paired += (p_table[k] == j && p_table[j] == k);
			auto const [vd,d] = decode(it->third);
			energy_t vkjwm = vd + params->MLintern[2];
			wm2 = std::min( wm2, WM[k-1]  + vkjwm );
		}
		if(p_table[j]<0 && !paired) wm2 = std::min(wm2, temp[j-1] + params->MLbase);
		temp[j] = wm2;
	}
	return temp;
}

/**
 * @brief Test existence of candidate
 * 
 * @param CL Candidate List
 * @param cand_comp 
 * @param i start
 * @param j end
 * @return true 
 * @return whether (i,j) is candidate for W/WM splits 
 */
bool is_candidate(auto const& CL,auto const& cand_comp,size_t i, size_t j) {
	const cand_list_td1 &list = CL[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,cand_comp);

	return it!=list.end() && it->first==i;
}
/**
 * Within the traceback, this function finds whether the multiloop has a dangle on the closing pair
*/
void find_mb_dangle(const energy_t &vkj,const energy_t &vk1j,const energy_t &vkj1,const energy_t &vk1j1,auto const &params, auto const& S, const size_t &i, const size_t &j, size_t &k, size_t &l,const int32_t* p_table){
	if(params->model_details.dangles == 2) return;

	int tt = pair[S[j]][S[i]];
	energy_t e1 = vkj +  E_MLstem(tt, -1, -1, params);
	energy_t e2 = vk1j +  E_MLstem(tt, -1, S[i+1], params);
	energy_t e3 = vkj1 +  E_MLstem(tt, S[j-1], -1, params);
	energy_t e4 = vk1j1 +  E_MLstem(tt, S[j-1], S[i+1], params);
	energy_t e = e1;
	if(e2<e && p_table[i+1]< 0){
		e = e2;
		k = i+2;
		l = j-1;
	}
	if(e3<e && p_table[j-1]< 0){
		e = e3;
		k = i+1; 
		l = j-2;
	}
	if(e4<e && p_table[i+1]< 0 && p_table[j-1]< 0){
		e = e4;
		k = i+2;
		l = j-2;
	}
 }

/**
 * @brief Trace from W entry
 * 
 * @param seq Sequence
 * @param CL Candidate List
 * @param cand_comp Candidate Comparator
 * @param structure Final structure
 * @param params Parameters
 * @param S Sequence Encoding
 * @param S1 Sequence Encoding
 * @param ta trace arrows
 * @param W W array
 * @param WM WM array
 * @param WM2 WM2 array
 * @param n Length
 * @param mark_candidates Whether candidates are marked as [ ]
 * @param i row index
 * @param j column index
 * @param p_table Restricted Array
 * @param last_j_array Restricted Array
 * @param in_pair_array Restricted Array
 * pre: W contains values of row i in interval i..j
 */
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,int32_t* p_table, int32_t* up_array) {

	// std:: cout << "W at  " << i << " and " << j << std::endl;
	if (i+TURN+1>=j) return;
	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,j-1,p_table,up_array);
		return;
	}
	size_t m = i;
	size_t l = j; 
	energy_t v=INF;
	int dangle = 3;
	// energy_t vv = INF;
	energy_t w;
    // int d =4;
	int ptype = 0;
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i;++it ) {
		size_t k = it->first;
		auto const[vd,d] = decode(it->third);
		energy_t vkjw = vd;
		w = W[k-1] + vkjw;
		if (W[j] == w) {
		v =vd;
		dangle = d;
		m = k;
		break;
		}
	}
	size_t k = m;
	switch(dangle){
			case 1:
				ptype = pair[S[m+1]][S[l]];
				v = v - E_ExtLoop(ptype,S[m],-1,params);
				++m;
				break;
			case 2:
				ptype = pair[S[m]][S[j-1]];
				v = v - E_ExtLoop(ptype,-1,S[j],params);
				l = j-1;
				break;
			case 3:
				if(params->model_details.dangles == 1){
					ptype = pair[S[m]][S[j-1]];
					v = v - E_ExtLoop(ptype,S[m],S[j],params);
					++m;
					l = j-1;
				} else{
					ptype = pair[S[m]][S[j]];
					int sk1 = m>1 ? S[m] : -1;
            		int sj1 = j<n ? S[j+1] : -1;
					v = v - E_ExtLoop(ptype,sk1,sj1,params);
				}
				break;  
    	}
    
	assert(i<=k && k<j);
	assert(v<INF);
	// don't recompute W, since i is not changed
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,k-1,p_table,up_array);
	trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,m,l,v,p_table,up_array);
}

/**
* @brief Trace from V entry
* 
* @param seq Sequence
* @param CL Candidate List
* @param cand_comp Candidate Comparator
* @param structure Final Structure
* @param params Parameters
* @param S Sequence Encoding
* @param S1 Sequence Encoding
* @param ta Trace Arrows
* @param WM WM array
* @param WM2 WM2 array
* @param n Length
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
* @param e energy in V[i,j]
* @param p_table Restricted Array
* @param last_j_array Restricted Array
* @param in_pair_array Restricted Array
* pre: structure is string of size (n+1)
*/
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,int32_t* p_table,int32_t* up_array) {
	// std::cout << "V at " << i << " and " << j << " with " << e << std::endl;

	assert( i+TURN+1<=j );
	assert( j<=n );
	
	
	if (mark_candidates && is_candidate(CL,cand_comp,i,j)) {
		structure[i]='{';
		structure[j]='}';
	} else {
		structure[i]='(';
		structure[j]=')';
	}
	const int32_t ptype_closing = pair[S[i]][S[j]];
	if (exists_trace_arrow_from(ta,i,j)) {
		
		const TraceArrow &arrow = trace_arrow_from(ta,i,j);
		const size_t k=arrow.k(i,j);
		const size_t l=arrow.l(i,j);
		assert(i<k);
		assert(l<j);
		
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l, arrow.target_energy(),p_table,up_array);
		return;

	}
	else {

		// assert(ptype_closing>0);
		// try to trace back to a candidate: (still) interior loop case
		for ( size_t l=i; l<j; l++) {
			// Break if it's an assured dangle case
			for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
				const size_t k=it->first;
				
				if (e == it->second + ILoopE(S,S1,params,ptype_closing,i,j,k,l) ) {
			// if (  e == it->second + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)) ) {
				trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,it->second,p_table,up_array);
				return;
				}
				
			}
		}
	}
	// is this a hairpin?
	if ( e == HairpinE(seq,S,S1,params,i,j) ) {
		return;
	}
	
	// if we are still here, trace to wm2 (split case);
	// in this case, we know the 'trace arrow'; the next row has to be recomputed
	auto const temp = recompute_WM(WM,CL,S,params,n,i+1,j-1,p_table,up_array);
	WM = temp;
	auto const temp2 = recompute_WM2(WM,WM2,CL,S,params,n,i+1,j-1,p_table);
	WM2 = temp2;
	
	// Dangle for Multiloop
	size_t k = i+1;
	size_t l = j-1;
	if(params->model_details.dangles == 1){
		auto const temp3 = recompute_WM(WM,CL,S,params,n,i+2,j-1,p_table,up_array);
		auto const temp4 = recompute_WM2(temp3,WM2,CL,S,params,n,i+2,j-1,p_table);
		find_mb_dangle(temp2[j-1],temp4[j-1],temp2[j-2],temp4[j-2],params,S,i,j,k,l,p_table);
		if (k>i+1){
			WM = temp3;
			WM2 = temp4;
		}
	}
	
	
	trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,p_table,up_array);
}

/**
* @brief Trace from WM
* 
* @param seq Sequence
* @param CL Candidate List
* @param cand_comp Candidate Comparator
* @param structure Final Structure
* @param params Parameters
* @param S Sequence Encoding
* @param S1 Sequence Encoding
* @param ta Trace Arrows
* @param WM WM array
* @param WM2 Wm2 array
* @param n Length
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index 
* @param e energy in WM[i,j] 
* @param p_table Restricted array
* @param last_j_array Restricted array
* @param in_pair_array Restricted array
* @param dangles Determines Multiloop Contribution
* pre: vector WM is recomputed for row i
*/
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j, energy_t e, int32_t* p_table, int32_t* up_array) {
	// std::cout << "WM at " << i << " and " << j << " with " << e << std::endl;

	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,WM[j-1],p_table,up_array);
		return;
	}
	energy_t v = INF;
	energy_t wm_split = INF;
	energy_t c_split = INF;
	size_t m = i;
	size_t l = j;
	int ptype = 0;
    int dangle = 3;
	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		size_t k = it->first;
		auto const[vd,d] = decode(it->third);
		energy_t vkjwm = vd + params->MLintern[2];
		
		if ( e == WM[k-1] + vkjwm ) {
		// no recomp, same i
			wm_split = WM[k-1] + vkjwm;
			m = k;
			v = vkjwm;
			dangle = d;
			break;
		} else if ( e == static_cast<energy_t>((k-i)*params->MLbase) + vkjwm ) {
			c_split = static_cast<energy_t>((k-i)*params->MLbase) + vkjwm;
			m = k;
			v= vkjwm;
			dangle = d;
			break;
		}
	}
	size_t k = m;
	switch(dangle){
			case 0:
				ptype = pair[S[m]][S[j]];
				v = v - E_MLstem(ptype,-1,-1,params);
				break;
			case 1:
				ptype = pair[S[m+1]][S[l]];
				v = v - E_MLstem(ptype,S[m],-1,params);
				++m;
				break;
			case 2:
				ptype = pair[S[m]][S[j-1]];
				v = v - E_MLstem(ptype,-1,S[j],params);
				l = j-1;
				break;
			case 3:
				if(params->model_details.dangles == 1){
					ptype = pair[S[m+1]][S[j-1]];
					v = v - E_MLstem(ptype,S[m],S[j],params);
					++m;
					l = j-1;
				} else{
					ptype = pair[S[m]][S[j]];
					v = v - E_MLstem(ptype,S[m-1],S[j+1],params);
				}
				break;  
    	}
	if ( e == wm_split ) {
		// no recomp, same i
			trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1],p_table,up_array);
			trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,m,l,v,p_table,up_array);
		return;
		} else if ( e == c_split) {
			trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,m,l,v,p_table,up_array);
			return;
		}
	
	assert(false);
}

/**
* @brief Trace from WM2
* 
* @param seq Sequence
* @param CL Candidate List
* @param cand_comp Candidate Comparator
* @param structure Final Structure
* @param params Parameters
* @param S Sequence Encoding
* @param S1 Sequence Encoding
* @param ta Trace Arrows
* @param WM WM array
* @param WM2 Wm2 array
* @param n Length
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
* @param p_table Restricted array
* @param last_j_array Restricted array
* @param in_pair_array Restricted array
* pre: vectors WM and WM2 are recomputed for row i
 */
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,int32_t* p_table, int32_t* up_array) {
	// std::cout << "WM2 at " << i << " and " << j << " with " << WM2[j] << std::endl;

	if (i+2*TURN+3>j) {return;}

	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		
		// same i, no recomputation
		trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,p_table,up_array);
		return;
	}
	size_t m = i;
	size_t l = j;
	int ptype = 0;
    energy_t v = INF;
	energy_t wm2_split = INF;
    int dangle = 3;
	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		size_t k = it->first;
		auto const[vd,d] = decode(it->third);
		energy_t vkjwm = vd + params->MLintern[2];
		// printf("k is %lu and j is %lu and e is %d and WM[k-1] is %d and vkjwm is %d and wm2 is %d\n",k,j,e,WM[k-1],vkjwm,WM[k-1] + vkjwm);
		if ( e == WM[k-1] + vkjwm ) {
			dangle = d;
			m = k;
			wm2_split = WM[k-1] + vkjwm;
			v = vkjwm;
			break;
		}
	}
	size_t k = m;
	switch(dangle){
		case 0:
			ptype = pair[S[m]][S[j]];
			v = v - E_MLstem(ptype,-1,-1,params);
			break;
		case 1:
			ptype = pair[S[m+1]][S[l]];
			v = v - E_MLstem(ptype,S[m],-1,params);
			++m;
			break;
		case 2:
			ptype = pair[S[m]][S[j-1]];
			v = v - E_MLstem(ptype,-1,S[j],params);
			l = j-1;
			break;
		case 3:
			if(params->model_details.dangles == 1){
				ptype = pair[S[m+1]][S[j-1]];
				v = v - E_MLstem(ptype,S[m],S[j],params);
				++m;
				l = j-1;
			} else{
				ptype = pair[S[m]][S[j]];
				v = v - E_MLstem(ptype,S[m-1],S[j+1],params);
			}
			break;  
    }
	if ( e == wm2_split ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1],p_table,up_array);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,m,l,v,p_table,up_array);
		return;
	}
	assert(false);
}
/**
* @brief Trace back
* pre: row 1 of matrix W is computed
* @return mfe structure (reference)
*/
const std::string & trace_back(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n,int32_t* p_table, int32_t *up_array,auto const& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,1,n,p_table,up_array);
	structure = structure.substr(1,n);

	return structure;
}

/* pre: ptype_closing>0 */
energy_t ILoopE(auto const& S, auto const& S1, auto const& params, const int32_t& ptype_closing,const size_t &i,const size_t &j,const size_t &k,const size_t &l)  {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);
	//assert(l<=len); // don't know len here

	// note: enclosed bp type 'turned around' for lib call
	const int32_t ptype_enclosed = rtype[pair[S[k]][S[l]]];

	if (ptype_enclosed==0) return INF;

	return E_IntLoop(k-i-1,j-l-1,ptype_closing,ptype_enclosed,S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
}


/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
*/
void register_candidate(auto &CL, size_t const& i, size_t const& j, energy_t const& e,energy_t const& ed) {
	assert(i<=j+TURN+1);
	
	CL[j].push_back( cand_entry_td1(i,e,ed) );
}


// /**
// * @brief Register a candidate
// * @param i start
// * @param j end
// * @param e energy of candidate "V(i,j)"
// */
// void register_candidate(auto &CL, size_t i, size_t j, energy_t e) {
// 	assert(i<=j+TURN+1);
// 	CL[j].push_back( cand_entry_t(i, e) );
// }

/**
 * @brief Evaluates whether a pairing can occur based on the restriction
 * 
 * @param i Current i
 * @param j Current j
 * @param p_table Restricted array
 * @param last_j_array Restricted array
 * @param in_pair_array Restricted array
 * @param multiloop Boolean to check if we are looking at WM and WM2
 * @return whether i and j can be non INF 
 */
bool evaluate_restriction(int i, int j, int32_t const*last_j_array, int32_t const*in_pair_array){
	bool evaluate = 1;
	if(in_pair_array[i]>in_pair_array[j]) evaluate = 0;

	if(in_pair_array[i]<in_pair_array[j]) evaluate = 0;

	if(in_pair_array[i]==in_pair_array[j]){
		if(j>last_j_array[i]) evaluate = 0;
	}
	return evaluate;
}

energy_t fold(auto const& seq, auto &V, auto const& cand_comp, auto &CL, auto const& S, auto const& S1, auto const& params, auto &ta, auto &W, auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto const& n, auto const& garbage_collect,int const*p_table,int const*last_j_array,int const*in_pair_array,int const*up_array) {
	int d = 3;
    if(params->model_details.dangles == 0) d = 0;
    
    for (size_t i=n; i>0; --i) {
		int si1 = (i>1) ? S[i-1] : -1;
		for ( size_t j=i+TURN+1; j<=n; j++ ) {

			bool evaluate = evaluate_restriction(i,j,last_j_array,in_pair_array);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			energy_t wm_split = INF;
			energy_t wm2_split = INF;
			for ( auto const [key,val,vald] : CL[j] ) {
				size_t k=key;
				// bool unpairedkj = (p_table[k]<-1 && p_table[j]<-1);
				auto const [v,d] = decode(vald);
				
				
				
				// auto const [v_kjw,dw] = decode(val_w);
				bool can_pair = up_array[k-1] >= (k-i) ? true: false;
				const energy_t v_kj = v + params->MLintern[2];
				wm_split = std::min( wm_split, WM[k-1] + v_kj );
				if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
				wm2_split = std::min( wm2_split, WM[k-1] + v_kj);
				// if(i==1 && j==82) printf("k is %lu and j is %lu and wsplit is %d and wk-1 is %d and vkjw is %d and d is %d\n",k,j,W[k-1] + vkjw,W[k-1],vkjw,d);
				w_split = std::min( w_split, W[k-1] + v );
		
				
			}
			if(p_table[j]<0) w_split = std::min(w_split,W[j-1]);
			if(p_table[j]<0) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			if(p_table[j]<0) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
			
			
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V


			size_t i_mod=i%(MAXLOOP+1);

			const int32_t ptype_closing = pair[S[i]][S[j]];
			const bool restricted = p_table[i] == -1 || p_table[j] == -1;

			const bool unpaired = (p_table[i]<-1 && p_table[j]<-1);
			const bool paired = (p_table[i] == j && p_table[j] == i);
			energy_t v = INF;
			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair
				bool canH = (paired || unpaired);
				if(up_array[j-1]<(j-i-1)) canH=false;
				
				energy_t v_h = canH ? HairpinE(seq,S,S1,params,i,j) : INF;
				// info of best interior loop decomposition (if better than hairpin)
				size_t best_l=0;
				size_t best_k=0;
				energy_t best_e;

				energy_t v_iloop=INF;

				// constraints for interior loops
				// i<k; l<j
				// k-i+j-l-2<=MAXLOOP  ==> k <= MAXLOOP+i+1
				//            ==> l >= k+j-i-MAXLOOP-2
				// l-k>=TURN+1         ==> k <= j-TURN-2
				//            ==> l >= k+TURN+1
				// j-i>=TURN+3
				//
				size_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
				// #pragma omp parallel for
				if((p_table[i]<-1 && p_table[j] < -1) || p_table[i] == j) { 
					for ( size_t k=i+1; k<=max_k; k++) {
						size_t k_mod=k%(MAXLOOP+1);
						
						bool cank = up_array[k-1]>=(k-i-1) ? true : false;
						size_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;
						for (size_t l=j-1; l>=min_l; --l) {
							assert(k-i+j-l-2<=MAXLOOP);
							if(V(k_mod,l) == INF) continue;
							
							const energy_t v_iloop_kl = cank && up_array[j-1]>=(j-l-1) ? V(k_mod,l) + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)) : INF;
							if ( v_iloop_kl < v_iloop) {
								v_iloop = v_iloop_kl;
								best_l=l;
								best_k=k;
								best_e=V(k_mod,l);
							}
						}
					}
				}
				
				
				const energy_t v_split = E_MbLoop(dmli1,dmli2,S,params,i,j,p_table);

				v = std::min(v_h,std::min(v_iloop,v_split));
				// register required trace arrows from (i,j)
				if ( v_iloop < std::min(v_h,v_split) ) {
					if ( is_candidate(CL,cand_comp,best_k,best_l) ) {
						//std::cout << "Avoid TA "<<best_k<<" "<<best_l<<std::endl;
						avoid_trace_arrow(ta);
					} else {
						//std::cout<<"Reg TA "<<i<<","<<j<<":"<<best_k<<","<<best_l<<std::endl;
						
						register_trace_arrow(ta,i,j,best_k,best_l,best_e);
					}
				}

				V(i_mod,j) = v;
			} else {
				V(i_mod,j) = INF;
			} // end if (i,j form a canonical base pair)

			
			
            
 			size_t ip1_mod = (i+1)%(MAXLOOP+1);
			energy_t vi1j = V(ip1_mod,j);
			energy_t vij1 = V(i_mod,j-1);
			energy_t vi1j1 = V(ip1_mod,j-1);	
			
			// Checking the dangle positions for W
			if(params->model_details.dangles == 1) d =0;
			
			energy_t w_v  = E_ext_Stem(v,vi1j,vij1,vi1j1,S,params,i,j,d,n,p_table);
			// Checking the dangle positions for W
			const energy_t wm_v = E_MLStem(v,vi1j,vij1,vi1j1,WM,CL,S,params,i,j,d,n,p_table);
			int k = i;
            int l = j;
			if(params->model_details.dangles == 1){
                if(d>0){
                    switch(d){
                        case 1:
                            k = i+1;
							break;
                        case 2:
                            l = j-1;
							break;
                        case 3: 
                            k = i+1;
                            l = j-1;
							break;
                    }
                    if(exists_trace_arrow_from(ta,k,l) && (wm_v < wm_split || w_v < w_split)) inc_source_ref_count(ta,k,l);	
                }
			}
			
			w  = std::min(w_v, w_split);
			wm = std::min(wm_v, wm_split);
			
			if ( w_v < w_split || wm_v < wm_split || paired) {
				int k_mod = k%(MAXLOOP+1);
				register_candidate(CL, i, j,V(i_mod,j),encode(w_v,d));
				// always keep arrows starting from candidates
				inc_source_ref_count(ta,i,j);
			}		

			W[j]       = w;
			WM[j]      = wm;
			WM2[j]     = wm2_split;

		} // end loop j
		rotate_arrays(WM2,dmli1,dmli2);

		// Clean up trace arrows in i+MAXLOOP+1
		if (garbage_collect && i+MAXLOOP+1 <= n) {
            gc_row(ta,i + MAXLOOP + 1 );
		}
		// Reallocate candidate lists in i
		for ( auto &x: CL ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_td1 vec(x.size());
				copy(x.begin(),x.end(),vec.begin());
				vec.swap(x);
			}
		}

		compactify(ta);
	}
	return W[n];
}

/**
 * @brief Fills the restriction arrays
 * p_table will contain the index of each base pair
 * X or x tells the program the base cannot pair and . sets it as unpaired but can pair
 * Pseudoknots (denoted by [ ], < >, or { } ) are filled the same way as ( )
 * That is, a structure like this (<)> is not possible.
 * @param structure Input structure
 * @param p_table Restricted array
 * @param last_j_array Restricted Array
 * @param in_pair_array Restricted Array
 */
void detect_restricted_pairs(auto const &structure, int32_t *p_table, int32_t *last_j_array, int32_t *in_pair_array){
	int i, j, count = 0, length = structure.length(),last_j=length;
	std::vector<int>  pairs;
	pairs.push_back(length);

	for (i=length; i >=1; --i){
		if ((structure[i-1] == 'x') || (structure[i-1] == 'X'))
			p_table[i] = -1;
		else if (structure[i-1] == '.')
			p_table[i] = -2;
		if (structure[i-1] == ')' || structure[i-1] == ']' || structure[i-1] == '}' || structure[i-1] == '>'){
			pairs.push_back(i);
			count++;
		}
		last_j_array[i] = pairs[pairs.size()-1];
		in_pair_array[i] = count;
		if (structure[i-1] == '(' || structure[i-1] == '[' || structure[i-1] == '{' || structure[i-1] == '<'){
			j = pairs[pairs.size()-1];
			pairs.erase(pairs.end()-1);
			p_table[i] = j;
			p_table[j] = i;
			count--;
		}
	}
	pairs.pop_back();
	if (pairs.size() != 0)
	{
		fprintf (stderr, "The given structure is not valid: more left parentheses than right parentheses: \n");
		exit (1);
	}
}

/**
 * @brief Sums the number of Candidates at each index over all indices
 * 
 * @param CL_ Candidate list
 * @return total number of candidates
 */
size_t num_of_candidates(auto const& CL_)  {
	size_t c=0;
	for ( auto const &x: CL_ ) {
		c += x.size();
	}
	return c;
}
/**
 * @brief Finds the size of allocated storage capacity across all indices
 * 
 * @param CL_ Candidate List
 * @return the amount of allocated storage 
 */
size_t capacity_of_candidates(auto const& CL_) {
	size_t c=0;
	for ( auto const &x: CL_ ) {
		c += x.capacity();
	}
	return c;
}

/**
* @brief Simple driver for @see SparseMFEFold.
*
* Reads sequence from command line or stdin and calls folding and
* trace-back methods of SparseMFEFold.
*/
int
main(int argc,char **argv) {

	args_info args_info;

	// get options (call gengetopt command line parser)
	if (cmdline_parser (argc, argv, &args_info) != 0) {
	exit(1);
	}

	std::string seq;
	if (args_info.inputs_num>0) {
	seq=args_info.inputs[0];
	} else {
	std::getline(std::cin,seq);
	}
	int n = seq.length();

	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (n,'.');

	if(restricted != "" && restricted.length() != n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		exit(0);
	}

	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		// FILE *fp;
    	// fp = fopen(parameter_file.c_str(),"r");
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}

	

	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;

	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted);

	if(args_info.dangles_given) sparsemfefold.params_->model_details.dangles = dangles;

	// Make replicate mx array in linear space
	int32_t last_j_array[n+1] = {0};
	int32_t in_pair_array[n+1] = {0};
	int32_t p_table[n+1] = {0};
	int32_t up_array[n+1] = {0};
	

	std::cout << seq << std::endl;
	
	detect_restricted_pairs(restricted,p_table,last_j_array,in_pair_array);

	int temp = 0;
	for(int i=1;i<=n;++i){
		if(p_table[i]>0) temp = 0;
		up_array[i]= temp;
		++temp;
		
	}
	
	energy_t mfe = fold(sparsemfefold.seq_,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.n_,sparsemfefold.garbage_collect_, p_table,last_j_array,in_pair_array,up_array);		
	// std::cout << mfe << std::endl;
	std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_,p_table,up_array, mark_candidates);
	
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

	std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

	// float factor=1024;
	
	// const std::string unit=" kB";
	
	
	
	if (verbose) {
		

	std::cout <<std::endl;

	std::cout << "TA cnt:\t"<<sizeT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TA max:\t"<<maxT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TA av:\t"<<avoidedT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TA rm:\t"<<erasedT(sparsemfefold.ta_)<<std::endl;

	std::cout <<std::endl;
	std::cout << "Can num:\t"<<num_of_candidates(sparsemfefold.CL_)<<std::endl;
	std::cout << "Can cap:\t"<<capacity_of_candidates(sparsemfefold.CL_)<<std::endl;
	std::cout << "TAs num:\t"<<sizeT(sparsemfefold.ta_)<<std::endl;
	std::cout << "TAs cap:\t"<<capacityT(sparsemfefold.ta_)<<std::endl;
	}

	return 0;
}
