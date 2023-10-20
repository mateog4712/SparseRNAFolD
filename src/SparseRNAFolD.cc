/**
* @mainpage
*
* Space-efficient sparse variant of an RNA (loop-based) free energy
* minimization algorithm (RNA folding equivalent to the Zuker
* algorithm).
*
* The results are equivalent to RNAfold.
*/
#define NDEBUG
#include "base_types.hh"
#include "cmdline.hh"
#include "matrix.hh"
#include "trace_arrow.hh"
#include <iostream>
#include <iomanip>
#include <vector>
#include <iterator>
#include <cstring>
#include <string>
#include <cassert>
#include <sstream>


extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
#include "ViennaRNA/params/io.h"
}

// New Candidate structure
struct quatret
{
    cand_pos_t first; 
    energy_t second;
    energy_t third;
	energy_t fourth;
	quatret(){
		first = 1;
		second = 2;
		third = 3;
		fourth = 4;
	}
	quatret(cand_pos_t x, energy_t y , energy_t z,energy_t w){
		first = x;
		second = y;
		third = z;
		fourth = w;
	}
};

typedef quatret cand_entry_t;
typedef std::vector< cand_entry_t > cand_list_t;

energy_t ILoopE(const short* S, const short* S1, const paramT* params, const pair_type& ptype_closing,const cand_pos_t &i,const cand_pos_t &j,const cand_pos_t &k,const cand_pos_t &l);

class SparseMFEFold;

/**
* Space efficient sparsification of Zuker-type RNA folding with
* trace-back. Provides methods for the evaluation of dynamic
* programming recursions and the trace-back.
*/
class SparseMFEFold {
	public:
		std::string seq_;
		cand_pos_t n_;

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
		std::vector< cand_list_t > CL_;	

		// compare candidate list entries by keys (left index i) in descending order
		struct Cand_comp{
		bool operator ()(const cand_entry_t &x, size_t y) const {
			return x.first > y;
		}
		} cand_comp;

		SparseMFEFold(const std::string &seq, bool garbage_collect, std::string restricted)
		: seq_(seq),n_(seq.length()),params_(scale_parameters()),ta_(n_),garbage_collect_(garbage_collect){
		
		make_pair_matrix();

		S_ = encode_sequence(seq.c_str(),0);
		S1_ = encode_sequence(seq.c_str(),1);

		V_.resize(MAXLOOP+1,n_+1,INF);
		W_.resize(n_+1,0);

		WM_.resize(n_+1,INF);

		WM2_.resize(n_+1,INF);

		dmli1_.resize(n_+1,INF);

		dmli2_.resize(n_+1,INF);

		// init candidate lists
		CL_.resize(n_+1);

		resize(ta_,n_+1);

		restricted_ = restricted;
		
		}

		~SparseMFEFold() {
		free(params_);
		free(S_);
		free(S1_);
		}
};

void trace_V(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array);
void trace_W(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta, const std::vector<energy_t> & W, std::vector< energy_t > &WM, std::vector< energy_t > &WM2, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array);
void trace_WM(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params,const short* S, const short* S1, TraceArrows &ta, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array) ;
void trace_WM2(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array);

bool evaluate_restriction(cand_pos_t i, cand_pos_t j, const std::vector<cand_pos_t> &last_j_array, const std::vector<cand_pos_t> &in_pair_array);


/**
 * @brief This code returns the hairpin energy for a given base pair.
 * @param i The left index in the base pair
 * @param j The right index in the base pair
*/
energy_t HairpinE(const std::string& seq, const short* S, const short* S1,  const paramT* params, cand_pos_t i, cand_pos_t j) {
	
	const pair_type ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
}

/**
 * @brief Returns the internal loop energy for a given i.j and k.l
 * 
*/
energy_t ILoopE(const short* S, const short* S1, const paramT* params, const pair_type& ptype_closing,const cand_pos_t &i,const cand_pos_t &j,const cand_pos_t &k,const cand_pos_t &l)  {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);

	// note: enclosed bp type 'turned around' for lib call
	const pair_type ptype_enclosed = rtype[pair[S[k]][S[l]]];

	// if (ptype_enclosed==0) return INF;

	return E_IntLoop(k-i-1,j-l-1,ptype_closing,ptype_enclosed,S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
}

/**
 * @brief Encodes the dangle type into the energy by shifting the energy two bits to the left and putting the dangle.
 * Dangles has 4 possible values: 00, 01, 10, 11, for no dangle, 5' dangle, 3' dangle, and 53' dangle, respectively.
 * 
 * @param e The energy being shifted. In practice, W(i,j) or WM(i,j)
 * @param d The type of dangle.
*/
// energy_t encode(energy_t e, Dangle d){
//     return (e << 2) | d;
// }
/**
 * @brief The complement to encode. This returns the energy and dangle, respectively, in a tuple format
 * 
 * @param enc The encoded energy 
*/
// std::pair<energy_t,Dangle> decode(energy_t enc){
//     return std::make_pair((enc >> 2), (enc & 3));
// }

/**
 * @brief Gives the W(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t E_ext_Stem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params, const cand_pos_t i,const cand_pos_t j, Dangle &d, cand_pos_t n, const std::vector<cand_pos_t>& p_table){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[i]][S[j]];
	
    if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j && p_table[j] == i)) {
				en = vij; // i j

				if (en != INF) {
					if (params->model_details.dangles == 2){
						base_type si1 = i>1 ? S[i-1] : -1;
                		base_type sj1 = j<n ? S[j+1] : -1;
                        en += vrna_E_ext_stem(tt, si1, sj1, params);
					}
                    else{
                        en += vrna_E_ext_stem(tt, -1, -1, params);
						d = 0;
					}

                    e = MIN2(e, en);
					
				}

	}

	if(params->model_details.dangles  == 1){
        tt  = pair[S[i+1]][S[j]];
        if (((p_table[i+1] <-1 && p_table[j] <-1) || (p_table[i+1] == j)) && p_table[i]<0) {
            en = (j-i-1>TURN) ? vi1j : INF; //i+1 j

            if (en != INF) {

                base_type si1 = S[i];
                en += vrna_E_ext_stem(tt, si1, -1, params);
            }

            e = MIN2(e,en);
            if(e == en){
                d=1;
            }

        }
        tt  = pair[S[i]][S[j-1]];
        if (((p_table[i] <-1 && p_table[j-1] <-1) || (p_table[i] == j-1)) && p_table[j]<0) {
            en = (j-1-i>TURN) ? vij1 : INF; // i j-1
            if (en != INF) {

                base_type sj1 = S[j];

                en += vrna_E_ext_stem(tt, -1, sj1, params);
            }
            e = MIN2(e,en);
            if(e == en){
                d=2;
            }

        }
        tt  = pair[S[i+1]][S[j-1]];
        if (((p_table[i+1] <-1 && p_table[j-1] <-1) || (p_table[i+1] == j-1)) && p_table[i] < 0 && p_table[j]<0) {
            en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1

            if (en != INF) {

                base_type si1 = S[i];
                base_type sj1 = S[j];

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
void rotate_arrays(std::vector<energy_t> &WM2, std::vector<energy_t> &dmli1, std::vector<energy_t> &dmli2){
	dmli2.swap(dmli1);
    dmli1.swap(WM2);
}


/**
* @brief Computes the multiloop V contribution. This gives back essentially VM(i,j).
* 
* @param dmli1 Row of WM2 from one iteration ago
* @param dmli2 Row of WM2 from two iterations ago 
*/
energy_t E_MbLoop( const std::vector<energy_t>& dmli1, const std::vector<energy_t>& dmli2, const short* S, paramT* params, cand_pos_t i, cand_pos_t j, const std::vector<cand_pos_t>& p_table){

	energy_t e = INF,en = INF;
  	pair_type tt  = pair[S[j]][S[i]];
	bool pairable = (p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j);
	
	/* double dangles */
	switch(params->model_details.dangles){
		case 2:
			if (pairable) {
			e = dmli1[j - 1];

			if (e != INF) {

				base_type si1 = S[i + 1];
				base_type sj1 = S[j - 1];

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

          			base_type si1 =  S[i + 1];

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
					base_type sj1 = S[j - 1];

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

					base_type si1 = S[i + 1];
					base_type sj1 = S[j - 1];

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
 * @brief Gives the WM(i,j) energy. The type of dangle model being used affects this energy. 
 * The type of dangle is also changed to reflect this.
 * 
 * @param vij The V(i,j) energy
 * @param vi1j The V(i+1,j) energy
 * @param vij1 The V(i,j-1) energy
 * @param vi1j1 The V(i+1,j-1) energy
*/
energy_t E_MLStem(const energy_t& vij,const energy_t& vi1j,const energy_t& vij1,const energy_t& vi1j1,const short* S, paramT* params,cand_pos_t i, cand_pos_t j,Dangle &d, const  cand_pos_t& n, const std::vector<cand_pos_t>& p_table){

	energy_t e = INF,en=INF;

	pair_type type = pair[S[i]][S[j]];


	if ((p_table[i] < -1 && p_table[j] < -1) || (p_table[i] == j)) {
		en = vij; // i j
		if (en != INF) {
			if (params->model_details.dangles == 2){
				base_type mm5 = i>1 ? S[i-1] : -1;
            	base_type mm3 = j<n ? S[j+1] : -1;
				en += E_MLstem(type, mm5, mm3, params);
			}
			else{
				en += E_MLstem(type, -1, -1, params);
				d = 0;
			}
			e = MIN2(e, en);
		}
	}
	if(params->model_details.dangles == 1){
		const base_type mm5 = S[i], mm3 = S[j];

		if (((p_table[i+1] < -1 && p_table[j] < -1) || (p_table[i+1] == j)) && p_table[i] < 0) {
      		en = (j-i-1 >TURN) ? vi1j : INF; // i+1 j
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
      		en = (j-1-i>TURN) ? vij1 : INF; // i j-1
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
      		en = (j-1-i-1>TURN) ? vi1j1 : INF; // i+1 j-1
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
* @brief Recompute row of WM. This is used in the traceback when we haved decided the current i.j pair closes a multiloop,
* and the WM energies need to be recomputed fom the candidates.
* 
* @param WM WM array
* @param CL Candidate List
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WM(const std::vector<energy_t>& WM, const std::vector< cand_list_t > &CL, const short* S, paramT* params, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array) {
	assert(i>=1);
	assert(max_j<=n_);
	std::vector<energy_t> temp = WM;

	for ( cand_pos_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }
	
	for ( cand_pos_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			cand_pos_t k = it->first;
			energy_t v_kj = it->third >> 2;

			bool can_pair = up_array[k-1] >= (k-i);
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wm = std::min( wm, temp[k-1]  + v_kj );
		}
		if(p_table[j]<0) wm = std::min(wm, temp[j-1] + params->MLbase);
		temp[j] = wm;
	}
	return temp;
}

/**
* @brief Recompute row of WM2. This is used in the traceback when we haved decided the current i.j pair closes a multiloop,
* and the WM2 energies need to be recomputed fom the candidates to get the corresponding energy for it.
* 
* @param WM WM array
* @param WM2 WM2 array
* @param CL Candidate List
* @param i Current i
* @param max_j Current j
*/
const std::vector<energy_t> recompute_WM2(const std::vector<energy_t>& WM,const std::vector<energy_t>& WM2, const std::vector< cand_list_t > &CL, const short* S, paramT* params, const cand_pos_t& n, cand_pos_t i, cand_pos_t max_j, const std::vector<cand_pos_t>& p_table) {
	assert(i>=1);
	assert(max_j<=n_);
	std::vector<energy_t> temp = WM2;

	for ( cand_pos_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	for ( cand_pos_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			cand_pos_t k = it->first;
			energy_t v_kj = it->third >> 2;
			
			wm2 = std::min( wm2, WM[k-1]  + v_kj );
		}
		if(p_table[j]<0) wm2 = std::min(wm2, temp[j-1] + params->MLbase);
		temp[j] = wm2;
	}
	return temp;
}

/**
 * @brief Test existence of candidate. Used primarily for determining whether (i,j) is candidate for W/WM splits
 * 
 * @param CL Candidate List
 * @param cand_comp Candidate Comparator
 * @param i start
 * @param j end
 * @return  
 */
bool is_candidate(const std::vector< cand_list_t >& CL,const SparseMFEFold::Cand_comp& cand_comp,cand_pos_t i, cand_pos_t j) {
	const cand_list_t &list = CL[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,cand_comp);

	return it!=list.end() && it->first==i;
}
/**
 * @brief Determines the type of dangle being used for a closing multiloop while in traceback.
 * 
 * @param WM2ij The WM2 energy for the region [i,j]
 * @param WM2i1j The WM2 energy for the region [i+1,j]
 * @param WM2ij1 The WM2 energy for the region [i,j-1]
 * @param WM2i1j1 The WM2 energy for the region [i+1,j-1]
*/
void find_mb_dangle(const energy_t WM2ij,const energy_t WM2i1j,const energy_t WM2ij1,const energy_t WM2i1j1, paramT* params, const short* S, const cand_pos_t i, const cand_pos_t j, cand_pos_t &k, cand_pos_t &l,const std::vector<cand_pos_t>& p_table){

	const pair_type tt = pair[S[j]][S[i]];
	const energy_t e1 = WM2ij +  E_MLstem(tt, -1, -1, params);
	const energy_t e2 = WM2i1j +  E_MLstem(tt, -1, S[i+1], params);
	const energy_t e3 = WM2ij1 +  E_MLstem(tt, S[j-1], -1, params);
	const energy_t e4 = WM2i1j1 +  E_MLstem(tt, S[j-1], S[i+1], params);
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
 * @brief Traceback from W entry.
 * pre: W contains values of row i in interval i..j
 * 
 * @param seq Sequence
 * @param structure Final structure
 * @param W W array
 * @param i row index
 * @param j column index
 */
void trace_W(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta, const std::vector<energy_t> & W, std::vector< energy_t > &WM, std::vector< energy_t > &WM2, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array) {
	if (i+TURN+1>=j) return;
	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,j-1,p_table,up_array);
		return;
	}
	cand_pos_t m=j+1;
	energy_t v=INF;
	energy_t w;
    Dangle dangle =3;
	energy_t vk = INF;
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i;++it ) {
		m = it->first;
		energy_t v_kj = it->fourth >> 2;
		Dangle d = it->fourth & 3;
		w = W[m-1] + v_kj;
		if (W[j] == w) {
		v =it->second;
        dangle = d;
		vk = v_kj;
		break;
		}
	}
	cand_pos_t k=m;
	cand_pos_t l=j;
	pair_type ptype = 0;
    switch(dangle){
		case 0:
			ptype = pair[S[k]][S[l]];
			v = vk - E_ExtLoop(ptype,-1,-1,params);
			break;
        case 1:
			
            k=m+1;
			ptype = pair[S[k]][S[l]];
			v = vk - E_ExtLoop(ptype,S[m],-1,params);
            break;
        case 2:
            l=j-1;
			ptype = pair[S[k]][S[l]];
			v = vk - E_ExtLoop(ptype,-1,S[j],params);
            break;
        case 3:
            if(params->model_details.dangles == 1){
                k=m+1;
                l=j-1;
				ptype = pair[S[k]][S[l]];
				v = vk - E_ExtLoop(ptype,S[m],S[j],params);
            }
            break;

        
    }
	assert(i<=m && m<j);
	assert(v<INF);
	// don't recompute W, since i is not changed
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,m-1,p_table,up_array);
	trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,p_table,up_array);
}

/**
* @brief Traceback from V entry
* 
* @param structure Final Structure
* @param mark_candidates Whether Candidates should be [ ]
* @param i row index
* @param j column index
*/
void trace_V(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array){
	assert( i+TURN+1<=j );
	if (mark_candidates && is_candidate(CL,cand_comp,i,j)) {
		structure[i]='{';
		structure[j]='}';
	} else {
		structure[i]='(';
		structure[j]=')';
	}
	const pair_type ptype_closing = pair[S[i]][S[j]];
	if (exists_trace_arrow_from(ta,i,j)) {
		
		const TraceArrow &arrow = trace_arrow_from(ta,i,j);
		const cand_pos_t k=arrow.k(i,j);
		const cand_pos_t l=arrow.l(i,j);
		assert(i<k);
		assert(l<j);
		
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l, arrow.target_energy(),p_table,up_array);
		return;

	}
	else {

		// try to trace back to a candidate: (still) interior loop case
		cand_pos_t l_min = std::max(i,j-31);
		for ( cand_pos_t l=j-1; l>l_min; l--) {
			for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
				const cand_pos_t k=it->first;
				if(k-i > 31) continue;
				// if (  e == it->second + ILoopE(S,S1,params,ptype_closing,i,j,k,l) ) {
				if (  e == it->second + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params)) ) {
					trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,it->second,p_table,up_array);
					return;
				}
			}
		}
	}
	// is it a hairpin?
	if ( e == HairpinE(seq,S,S1,params,i,j) ) {
		return;
	}
	
	// if we are still here, trace to wm2 (split case);
	// in this case, we know the 'trace arrow'; the next row has to be recomputed
	const std::vector<energy_t> temp = recompute_WM(WM,CL,S,params,n,i+1,j-1,p_table,up_array);
	WM = temp;
	const std::vector<energy_t> temp2 = recompute_WM2(WM,WM2,CL,S,params,n,i+1,j-1,p_table);
	WM2 = temp2;
	
	// Dangle for Multiloop
	cand_pos_t k = i+1;
	cand_pos_t l = j-1;
	if(params->model_details.dangles == 1){
		const std::vector<energy_t> temp3 = recompute_WM(WM,CL,S,params,n,i+2,j-1,p_table,up_array);
		const std::vector<energy_t> temp4 = recompute_WM2(temp3,WM2,CL,S,params,n,i+2,j-1,p_table);
		find_mb_dangle(temp2[j-1],temp4[j-1],temp2[j-2],temp4[j-2],params,S,i,j,k,l,p_table);
		if (k>i+1){
			WM = temp3;
			WM2 = temp4;
		}
	}

	trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,p_table,up_array);
}

/**
* @brief Traceback from WM
*
* @param WM WM array at [i,j]
* @param WM2 WM2 array at [i,j]
* @param i row index
* @param j column index 
* @param e energy in WM(i,j) 
*/
void trace_WM(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params,const short* S, const short* S1, TraceArrows &ta, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n, const bool& mark_candidates, cand_pos_t i, cand_pos_t j, energy_t e,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array){

	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,WM[j-1],p_table,up_array);
		return;
	}
	energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 3;
    cand_pos_t m = j+1;
	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		m = it->first;
		// auto const [v_kj,d] = decode(it->third);
		energy_t v_kj = it->third >> 2;
		Dangle d = it->third & 3;
		if ( e == WM[m-1] + v_kj ) {
            dangle = d;
			vk = v_kj;
			v = it->second;
			// no recomp, same i
		    break;
		} else if ( e == static_cast<energy_t>((m-i)*params->MLbase) + v_kj ) {
            dangle = d;
			vk = v_kj;
			v = it->second;
		    break;
		}
	}
	cand_pos_t k = m;
	cand_pos_t l = j;
	pair_type ptype = 0;
    switch(dangle){
		case 0:
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,-1,params);
			break;
        case 1:
            k=m+1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,S[m],-1,params) - params->MLbase;
            break;
        case 2:
            l=j-1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,S[j],params) - params->MLbase;
            break;
        case 3:
			if(params->model_details.dangles == 1){
				k=m+1;
				l=j-1;
				ptype= pair[S[k]][S[l]];
				v = vk - E_MLstem(ptype,S[m],S[j],params) - 2*params->MLbase;
			}
            break;
    }   

    if ( e == WM[m-1] + vk ) {
		// no recomp, same i
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,m-1,WM[m-1],p_table,up_array);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,p_table,up_array);
		return;
	} else if ( e == static_cast<energy_t>((m-i)*params->MLbase) + vk ) {
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,p_table,up_array);
		return;
	}
	assert(false);
}

/**
* @brief Traceback from WM2
* 
* @param WM WM array at [i,j]
* @param WM2 Wm2 array at [i,j]
* @param i row index
* @param j column index
 */
void trace_WM2(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1, TraceArrows &ta, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n, const bool& mark_candidates,cand_pos_t i, cand_pos_t j,const std::vector<cand_pos_t>& p_table, const std::vector<cand_pos_t>& up_array){

	if (i+2*TURN+3>j) {return;}
	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		
		// same i, no recomputation
		trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,p_table,up_array);
		return;
	}

    cand_pos_t m = j+1;
    energy_t v = INF;
    energy_t vk = INF;
    Dangle dangle = 3;
	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		m = it->first;
		
		energy_t v_kj = it->third >> 2;
		Dangle d = it->third & 3;
		if (e == WM[m-1] + v_kj) {
			vk = v_kj;
            dangle = d;
			v = it->second;
			break;
		}
	}
	cand_pos_t k = m;
	cand_pos_t l = j;
	pair_type ptype = 0;
    switch(dangle){
		case 0:
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,-1,params);
			break;
        case 1:
            k=m+1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,S[m],-1,params)-params->MLbase;
            break;
        case 2:
            l=j-1;
			ptype= pair[S[k]][S[l]];
			v = vk - E_MLstem(ptype,-1,S[j],params)-params->MLbase;
            break;
        case 3:
			if(params->model_details.dangles == 1){
				k=m+1;
				l=j-1;
				ptype= pair[S[k]][S[l]];
				v = vk - E_MLstem(ptype,S[m],S[j],params)-2*params->MLbase;
			}
            break;
    }
    

    if ( e == WM[m-1] + vk ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,m-1,WM[m-1],p_table,up_array);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,v,p_table,up_array);
		return;
	}
	assert(false);
}
/**
* @brief Trace back
* pre: row 1 of matrix W is computed
* @return mfe structure (reference)
*/
const std::string & trace_back(const std::string& seq, const std::vector< cand_list_t >& CL, const SparseMFEFold::Cand_comp& cand_comp, std::string &structure, paramT* params, const short* S, const short* S1,TraceArrows &ta, const std::vector<energy_t>& W, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, const cand_pos_t& n,const std::vector<cand_pos_t>& p_table,const std::vector<cand_pos_t>& up_array, const bool& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,1,n,p_table,up_array);
	structure = structure.substr(1,n);

	return structure;
}

/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
* @param wmij energy at WM(i,j)
* @param wij energy at W(i,j)
*/
void register_candidate(std::vector<cand_list_t> &CL, cand_pos_t const& i, cand_pos_t const& j, energy_t const& e,energy_t const& wmij,energy_t const& wij) {
	assert(i<=j+TURN+1);
	
	CL[j].push_back( cand_entry_t(i,e,wmij,wij) );
}

/**
 * @brief Evaluates whether a pairing can occur based on the restriction
 * 
 * @param i Current i
 * @param j Current j
 * @param p_table Restricted array
 * @param last_j_array Restricted array
 * @param in_pair_array Restricted array
 * @return whether i and j can be non INF 
 */
bool evaluate_restriction(cand_pos_t i, cand_pos_t j, const std::vector<cand_pos_t>& last_j_array, const std::vector<cand_pos_t>& in_pair_array){
	bool evaluate = 1;
	// if(in_pair_array[i]>in_pair_array[j]) evaluate = 0;
	evaluate &= ~(in_pair_array[i]>in_pair_array[j]);
	// if(in_pair_array[i]<in_pair_array[j]) evaluate = 0;
	evaluate &= ~(in_pair_array[i]<in_pair_array[j]);
	// if(in_pair_array[i]==in_pair_array[j]) if(j>last_j_array[i]) evaluate = 0
	evaluate &= ((in_pair_array[i]==in_pair_array[j]) & ~(j>last_j_array[i]));
	
	return evaluate;
}

/**
 * @brief Determines the MFE energy for a given sequence
*/
energy_t fold(const std::string& seq, LocARNA::Matrix<energy_t> &V, const SparseMFEFold::Cand_comp& cand_comp, std::vector<cand_list_t> &CL, const short* S, const short* S1, paramT* params, TraceArrows &ta, std::vector<energy_t> &W, std::vector<energy_t> &WM, std::vector<energy_t> &WM2, std::vector<energy_t> &dmli1, std::vector<energy_t> &dmli2, const cand_pos_t n, const bool garbage_collect,const std::vector<cand_pos_t>& p_table,const std::vector<cand_pos_t>& last_j_array,const std::vector<cand_pos_t>& in_pair_array,const std::vector<cand_pos_t>& up_array) {
	Dangle d = 3;
    if(params->model_details.dangles == 0 || params->model_details.dangles == 1) d = 0;
    
    for (cand_pos_t i=n; i>0; --i) {
		for (cand_pos_t j=i+TURN+1; j<=n; j++ ) {

			bool evaluate = evaluate_restriction(i,j,last_j_array,in_pair_array);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			energy_t wm_split = INF;
			energy_t wm2_split = INF;
			for ( auto it=CL[j].begin();CL[j].end() != it;++it ) {
				const cand_pos_t k=it->first;

				// Decode the energies
				const energy_t v_kj = it->third >> 2;
				const energy_t v_kjw = it->fourth >> 2;

				const bool can_pair = up_array[k-1] >= (k-i);
				
				wm_split = std::min( wm_split, WM[k-1] + v_kj );
				if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
				wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
				w_split = std::min( w_split, W[k-1] + v_kjw );			
			}

			if(p_table[j]<0) w_split = std::min(w_split,W[j-1]);
			if(p_table[j]<0) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			if(p_table[j]<0) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
			
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V


			cand_pos_t i_mod= i%(MAXLOOP+1);

			const pair_type ptype_closing = pair[S[i]][S[j]];
			const bool restricted = p_table[i] == -1 || p_table[j] == -1;

			const bool paired = (p_table[i] == j && p_table[j] == i);
			energy_t v = INF;
			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0 && evaluate && !restricted) { // if i,j form a canonical base pair
				bool canH = !up_array[j-1]<(j-i-1);
				
				energy_t v_h = canH ? HairpinE(seq,S,S1,params,i,j) : INF;
				// info of best interior loop decomposition (if better than hairpin)
				cand_pos_t best_l=0;
				cand_pos_t best_k=0;
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
				cand_pos_t max_k = std::min(j-TURN-2,i+MAXLOOP+1);
				if((p_table[i]<-1 && p_table[j] < -1) || p_table[i] == j) { 
					for ( cand_pos_t k=i+1; k<=max_k; k++) {
						cand_pos_t k_mod= k%(MAXLOOP+1);
						
						energy_t cank = ((up_array[k-1]>=(k-i-1))-1);;
						cand_pos_t min_l=std::max(k+TURN+1, k+j-i- MAXLOOP-2);
						
						for (size_t l=j-1; l>=min_l; --l) {
							assert(k-i+j-l-2<=MAXLOOP);
							energy_t canl = (((up_array[j-1]>=(j-l-1))-1) | cank);
							energy_t v_iloop_kl = INF & canl;
							v_iloop_kl = v_iloop_kl + V(k_mod,l) + E_IntLoop(k-i-1,j-l-1,ptype_closing,rtype[pair[S[k]][S[l]]],S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
							// v_iloop_kl = v_iloop_kl + V(k_mod,l) + ILoopE(S,S1,params,ptype_closing,i,j,k,l);

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
						avoid_trace_arrow(ta);
					} else {
						register_trace_arrow(ta,i,j,best_k,best_l,best_e);
					}
				}
				V(i_mod,j) = v;
			} else {
				V(i_mod,j) = INF;
			} // end if (i,j form a canonical base pair)			

 			cand_pos_t ip1_mod = (i+1)%(MAXLOOP+1);
			energy_t vi1j = V(ip1_mod,j);
			energy_t vij1 = V(i_mod,j-1);
			energy_t vi1j1 = V(ip1_mod,j-1);
				
			// Checking the dangle positions for W
			// if(params->model_details.dangles == 1) d =0;
			const energy_t w_v  = E_ext_Stem(v,vi1j,vij1,vi1j1,S,params,i,j,d,n,p_table);
			// Checking the dangle positions for W
			const energy_t wm_v = E_MLStem(v,vi1j,vij1,vi1j1,S,params,i,j,d,n,p_table);


			cand_pos_t k = i;
            cand_pos_t l = j;
			if(params->model_details.dangles == 1){
                if(d>0){
                    switch(d){
                        case 1:
                           	++k;
							break;
                        case 2:
                            --l;
							break;
                        case 3: 
                            ++k;
                            --l;
							break;
                    }
                    if(exists_trace_arrow_from(ta,k,l) && (wm_v < wm_split || w_v < w_split)) inc_source_ref_count(ta,k,l);	
                }
			}
			
			w  = std::min(w_v, w_split);
			wm = std::min(wm_v, wm_split);
			if ( w_v < w_split || wm_v < wm_split || paired) {
				// Encode the dangles into the energies
				energy_t w_enc = (w_v << 2) | d;
				energy_t wm_enc = (wm_v << 2) | d;
				register_candidate(CL, i, j,V(i_mod,j), wm_enc,w_enc);
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
		for ( cand_list_t &x: CL ) {
			if (x.capacity() > 1.5*x.size()) {
				cand_list_t vec(x.size());
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
void detect_restricted_pairs(const std::string &structure, std::vector<cand_pos_t> &p_table, std::vector<cand_pos_t> &last_j_array, std::vector<cand_pos_t> &in_pair_array){
	cand_pos_t i, j, count = 0, length = structure.length(),last_j=length;
	std::vector<cand_pos_t>  pairs;
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
cand_pos_t num_of_candidates(const std::vector<cand_list_t>& CL_)  {
	cand_pos_t c=0;
	for ( const cand_list_t &x: CL_ ) {
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
cand_pos_t capacity_of_candidates(const std::vector<cand_list_t>& CL_) {
	cand_pos_t c=0;
	for ( const cand_list_t &x: CL_ ) {
		c += x.capacity();
	}
	return c;
}
void seqtoRNA(std::string &sequence){
	bool DNA = false;
    for (char &c : sequence) {
      	if (c == 'T' || c == 't') {
			c = 'U';
			DNA = true;
		}
    }
	noGU = DNA;
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
	cand_pos_t n = seq.length();

	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (n,'.');

	if(restricted != "" && restricted.length() != n ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		exit(0);
	}

	std::string file= "";
	args_info.paramFile_given ? file = parameter_file : file = "";
	if(file!=""){
		vrna_params_load(file.c_str(), VRNA_PARAMETER_FORMAT_DEFAULT);
	}

	

	bool verbose = args_info.verbose_given;

	bool mark_candidates = args_info.mark_candidates_given;

	noGU = args_info.noGU_given;

	seqtoRNA(seq);

	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted);

	if(args_info.dangles_given) sparsemfefold.params_->model_details.dangles = dangle_model;

	// Make replicate mx array in linear space
	std::vector<cand_pos_t> last_j_array;
	std::vector<cand_pos_t> in_pair_array;
	std::vector<cand_pos_t> p_table;
	std::vector<cand_pos_t> up_array;

	last_j_array.resize(n+1,0);
	in_pair_array.resize(n+1,0);
	p_table.resize(n+1,0);
	up_array.resize(n+1,0);
	
	cmdline_parser_free(&args_info);

	std::cout << seq << std::endl;
	
	detect_restricted_pairs(restricted,p_table,last_j_array,in_pair_array);

	cand_pos_t temp = 0;
	for(cand_pos_t i=1;i<=n;++i){
		if(p_table[i]>0) temp = 0;
		up_array[i]= temp;
		++temp;
		
	}
	
	
	energy_t mfe = fold(sparsemfefold.seq_,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.n_,sparsemfefold.garbage_collect_, p_table,last_j_array,in_pair_array,up_array);		
	std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_,p_table,up_array, mark_candidates);
	
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

	std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;
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
