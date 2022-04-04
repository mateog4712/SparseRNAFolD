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

#include <LocARNA/matrix.hh>

#include <limits>

#include <vector>
#include <iterator>

#include <cstring>
#include <string>
#include <cassert>
#include <numeric>

#include "base.hh"
#include "trace_arrow.hh"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
}

#include "cmdline.hh"
#include <omp.h>

#include <fstream>


typedef unsigned short int cand_pos_t;
typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
typedef std::vector< cand_entry_t > cand_list_t;

class SparseMFEFold;

namespace unrolled {
template < typename Iterator, class Functor >
void for_each( Iterator start, Iterator end, Functor f, size_t num_times ) {
    for(Iterator cur = start; cur != end; cur += num_times ) {
        for( int i = 0; i < num_times; ++i ) {
            f( *(cur + i) );
        }
    }
}
}

energy_t ILoopE(auto const& S_,auto const& S1_, auto const& params_, int ptype_closing,size_t i, size_t j, size_t k,  size_t l);
energy_t MbLoopE(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
energy_t Mlstem(auto const& S_, auto const& params_, int ptype_closing,size_t i, size_t j);
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,int* p_table,int *last_j_array, int *in_pair_array,int dangles);
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,int* p_table,int *last_j_array, int *in_pair_array,int dangles);
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,int* p_table,int *last_j_array, int *in_pair_array,int dangles) ;
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,int* p_table,int *last_j_array, int *in_pair_array,int dangles);

bool evaluate_restriction(int i, int j, int *last_j_array, int *in_pair_array);

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
	int dangles_;
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

	
	

	/**
	candidate list for decomposition in W or WM

	@note Avoid separate candidate lists CLW and CLWM for split cases in W and
	WM to save even more space; here, this works after
	reformulating the recursions such that both split-cases recurse to
	V-entries. (compare OCTs)
	*/
	

	// compare candidate list entries by keys (left index i) in descending order
	struct {
	bool operator ()(const cand_entry_t &x, size_t y) const {
		return x.first > y;
	}
	}
	cand_comp;

	


	SparseMFEFold(const std::string &seq, bool garbage_collect, std::string restricted, int dangles)
	: seq_(seq),
	n_(seq.length()),
	params_(scale_parameters()),
	ta_(n_),
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

	restricted_ = restricted;

	dangles_ = dangles;

	
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


/**
* @brief Rotate WM2 arrays to store the previous and previous previous iterations
* 
* @param WM WM array
* @param WM2 WM2 array
* @param dmli1 WM2 from one iteration ago
* @param dmli2 WM2 from two iterations ago
* @param n Length
*/
void rotate_arrays(auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto n){
	

	for (int j = 1; j <= n; j++){
		dmli2[j] = dmli1[j];
		dmli1[j] = WM2[j];
	} 
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
* @param dangles Determines multiloop contribution
* @param p_table Restricted Array
* @return energy_t 
*/
energy_t E_MbLoop(auto const& dmli1, auto const& dmli2, auto const& S, auto const& params, size_t i, size_t j, auto const& dangles, auto const& p_table){

	int e = INF;
	int en = INF;
  	unsigned int tt;
	tt  = pair[S[j]][S[i]];

	/* double dangles */
	switch(dangles){
		case 2:
			if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) {
			e = dmli1[j - 1];

			if (e != INF) {

				int si1 = S[i + 1];
				int sj1 = S[j - 1];

				e += E_MLstem(tt, sj1, si1, params) + params->MLclosing;
				// e += E_MLstem(tt, -1, -1, params) + params->MLclosing;
			}

			}
			break;

		case 1:
			/* new closing pair (i,j) with mb part [i+1,j-1] */
			tt  = pair[S[j]][S[i]];
			if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) {
        		e = dmli1[j - 1];

        		if (e != INF) {

          			e += E_MLstem(tt, -1, -1, params) + params->MLclosing;

        		}
      		}
     		tt  = pair[S[j]][S[i]];
      		/* new closing pair (i,j) with mb part [i+2,j-1] */
      		if (((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) && p_table[i+1] < 0) {
        		en = dmli2[j - 1];

        		if (en != INF) {

          			int si1 =  S[i + 1];

          			en += E_MLstem(tt, -1, si1, params) + params->MLclosing + params->MLbase;
      
        		}
      		}
      		e   = MIN2(e, en);

			/* new closing pair (i,j) with mb part [i+1, j-2] */
			if (((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) && p_table[j-1] < 0) {
				en = dmli1[j - 2];

				if (en != INF) {
					int sj1 = S[j - 1];

					en += E_MLstem(tt, sj1, -1, params) + params->MLclosing + params->MLbase; 
				}
			}
			e   = MIN2(e, en);

			/* new closing pair (i,j) with mb part [i+2.j-2] */
			if (((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) && p_table[i+1] < -1 && p_table[j-2] <-1) {
				e = dmli2[j - 2];

				if (e != INF) {

					int si1 = S[i + 1];
					int sj1 = S[j - 1];

					e += E_MLstem(tt, sj1, si1, params) + params->MLclosing + 2 * params->MLbase;
				}
			}
			e   = MIN2(e, en);
      		break;
		case 3:
			if ((p_table[i] <-1 && p_table[j] <-1) || (p_table[i] == j and p_table[j] == i)) {
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
* @param WM WM array
* @param CL Candidate List
* @param S Sequence Encoding
* @param params Parameters
* @param i Current i
* @param j Current j
* @param n Length
* @param dangles Determines multiloop contribution
* @param p_table Restricted array
* @return energy_t 
*/
energy_t E_MLStem(auto const& vkj, auto const& WM, auto const& CL,auto const& S, auto const& params,size_t i, size_t j, auto const& n, auto const& dangles, auto const& p_table){

	int e = INF;

	int type = pair[S[i]][S[j]];
	
	switch(dangles){

		case 2:

			if ((p_table[i] < 0 && p_table[j] < 0) || (p_table[i] == j && p_table[j] == i)) {
				int en = vkj;
				if (en != INF) {
					if (dangles == 2)
						en += E_MLstem(type, (i == 1) ? S[n] : S[i - 1], S[j + 1], params);
					else
						en += E_MLstem(type, -1, -1, params);

					e = MIN2(e, en);
				}
			}
			break;

		case 1:
			/*
			*  extension with one unpaired nucleotide at the right (3' site)
			*  or full branch of (i,j)
			*/
			if ((p_table[i] < 0 && p_table[j] < 0) || (p_table[i] == j && p_table[j] == i)) {
				int en = vkj;
				if (en != INF) {
					if (dangles == 2)
						en += E_MLstem(type, (i == 1) ? S[n] : S[i - 1], S[j + 1], params);
					else
						en += E_MLstem(type, -1, -1, params);

					e = MIN2(e, en);
				}
			}
			/*
			*
			* i,j-1 should be WM2[j-1]
			*/
			if (p_table[j-1]<-1) {
				int en = WM[j-1];
				if (en != INF) {
					en += params->MLbase;

					e = MIN2(e, en);
				}
			}
			/*
			*  extension with one unpaired nucleotide at 5' site
			*  and all other variants which are needed for odd
			*  dangle models
			*  i+1 should be equal to WM from iteration before
			*/
			if (p_table[i+1]<-1) {
				int en = WM[j];
				if (en != INF) {
					en += params->MLbase;

					e = MIN2(e, en);
				}
			}
			break;
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
auto const recompute_W(auto const &W, auto const& CL, size_t i, size_t max_j, int* p_table) {
	
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
* @param dangles Determines multiloop calculation
* @return auto const 
*/
auto const recompute_WM(auto const& WM, auto const &CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, int* p_table,int dangles) {
	

	assert(i>=1);
	assert(max_j<=n);

	std::vector<energy_t> temp = WM;

	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }
	
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;
		bool paired;
		int mm3 = S[j-1];
		#pragma omp parallel for num_threads(6);
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
			paired = (p_table[k] == j && p_table[j] == k);
			int mm5 = S[k+1];
			const energy_t v_kj = E_MLStem(it->second,WM,CL,S,params,k,j,n,dangles,p_table);
			bool can_pair = true;
			for(int m = i;m<k;++m){
				if(p_table[m]>-1){
					can_pair = false;
					break;
				} 
		
			}
			if(can_pair) wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wm = std::min( wm, temp[k-1]  + v_kj );
			if(paired) break;
		}
		if(p_table[j]<-1) wm = std::min(wm, temp[j-1] + params->MLbase);
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
* @param dangles Determines calculation type
* @return auto const 
*/
auto const recompute_WM2(auto const& WM, auto const& WM2, auto const CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j, int* p_table,int *last_j_array, int *in_pair_array,int dangles) {
	

	assert(i>=1);
	//assert(i+2*TURN+3<=max_j);
	assert(max_j<= n);

	std::vector<energy_t> temp = WM2;

	for ( size_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	#pragma omp parallel for num_threads(6);
	for ( size_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;
		bool paired;
		int mm3 = S[j-1];
		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			
			size_t k = it->first;
			paired = (p_table[k] == j && p_table[j] == k);
			int mm5 = S[k+1];
			energy_t v_kl = E_MLStem(it->second,WM,CL,S,params,k,j,n,dangles,p_table);
			wm2 = std::min( wm2, WM[k-1]  + v_kl );
			if(paired) break;
		}
		if(p_table[j]<-1) wm2 = std::min(wm2, temp[j-1] + params->MLbase);
		// if(evaluate_restriction(i,j,last_j_array,in_pair_array)) wm2=INF;
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
	const cand_list_t &list = CL[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,cand_comp);

	return it!=list.end() && it->first==i;
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
 * @param dangles Determines Multiloop Contribution
 * pre: W contains values of row i in interval i..j
 */
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j,int* p_table,int *last_j_array, int *in_pair_array,int dangles) {
	
	if (i+TURN+1>=j) return;
	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,j-1,p_table,last_j_array,in_pair_array,dangles);
		return;
	}
	
	size_t k=j+1;
	energy_t v=INF;
	int sj1 = (j<n) ? S[j+1] : -1;
	energy_t w;
	for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i;++it ) {
		k = it->first;
		
		int sk1 = (k>1) ? S[k-1] : -1;
		const energy_t v_kj = it->second + vrna_E_ext_stem(pair[S[k]][S[j]],sk1,sj1,params);
		w = W[k-1] + v_kj;
		// std::cout << w << std::endl;
		
		if (W[j] == w) {
		v = it->second;
		break;
		}
	}

	assert(i<=k && k<j);
	assert(v<INF);

	// don't recompute W, since i is not changed
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,k-1,p_table,last_j_array,in_pair_array,dangles);
	trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,v,p_table,last_j_array,in_pair_array,dangles);
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
* @param dangles Determines Multiloop Contribution
* pre: structure is string of size (n+1)
*/
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S,auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e,int* p_table,int *last_j_array, int *in_pair_array,int dangles) {
	assert( i+TURN+1<=j );
	assert( j<=n );

	if (mark_candidates && is_candidate(CL,cand_comp,i,j)) {
		structure[i]='[';
		structure[j]=']';
	} else {
		structure[i]='(';
		structure[j]=')';
	}
	const int ptype_closing = pair[S[i]][S[j]];

	if (exists_trace_arrow_from(ta,i,j)) {
		// trace arrows may exist for interior loop case
		const TraceArrow &arrow = trace_arrow_from(ta,i,j);

		const size_t k=arrow.k(i,j);
		const size_t l=arrow.l(i,j);
		assert(i<k);
		assert(l<j);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l, arrow.target_energy(),p_table,last_j_array,in_pair_array,dangles);
		return;

	} else {

		assert(ptype_closing>0);

		// try to trace back to a candidate: (still) interior loop case
		for ( size_t l=i; l<j; l++) {
		for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
			const size_t k=it->first;
			if (  e == it->second + ILoopE(S,S1,params,ptype_closing,i,j,k,l) ) {
			trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,it->second,p_table,last_j_array,in_pair_array,dangles);
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
	auto const temp = recompute_WM(WM,CL,S,params,n,i+1,j-1,p_table,dangles);
	WM = temp;
	auto const temp2 = recompute_WM2(WM,WM2,CL,S,params,n,i+1,j-1,p_table,last_j_array,in_pair_array,dangles);
	WM2 = temp2;
	
	trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i+1,j-1,p_table,last_j_array,in_pair_array,dangles);
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
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j, energy_t e, int* p_table,int *last_j_array, int *in_pair_array,int dangles) {
	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,WM[j-1],p_table,last_j_array,in_pair_array,dangles);
		return;
	}
	int mm3 = S[j-1];
	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		const size_t k = it->first;
		int mm5 = S[k+1];
		const energy_t v_kj = E_MLStem(it->second,WM,CL,S,params,k,j,n,dangles,p_table);
		if ( e == WM[k-1] + v_kj ) {
		// no recomp, same i
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1],p_table,last_j_array,in_pair_array,dangles);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second,p_table,last_j_array,in_pair_array,dangles);
		return;
		} else if ( e == static_cast<energy_t>((k-i)*params->MLbase) + v_kj ) {
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second,p_table,last_j_array,in_pair_array,dangles);
		return;
		}
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
* @param dangles Determines Multiloop Contribution
* pre: vectors WM and WM2 are recomputed for row i
 */
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j,int* p_table,int *last_j_array, int *in_pair_array,int dangles) {
	if (i+2*TURN+3>j) {return;}

	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		
		// same i, no recomputation
		trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,p_table,last_j_array,in_pair_array,dangles);
		return;
	}
	int mm3 = S[j-1];
	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		size_t k = it->first;
		int mm5 = S[k+1];
		const energy_t v_kj = E_MLStem(it->second,WM,CL,S,params,k,j,n,dangles,p_table);
		if ( e == WM[k-1] + v_kj ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1],p_table,last_j_array,in_pair_array,dangles);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second,p_table,last_j_array,in_pair_array,dangles);
		return;
		}
	}
	assert(false);
}
/**
* @brief Trace back
* pre: row 1 of matrix W is computed
* @return mfe structure (reference)
*/
const std::string & trace_back(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n,int* p_table,int *last_j_array, int *in_pair_array,int dangles,auto const& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,1,n,p_table,last_j_array,in_pair_array,dangles);
	structure = structure.substr(1,n);

	return structure;
}

/* pre: ptype_closing>0 */
energy_t ILoopE(auto const& S, auto const& S1, auto const& params, int ptype_closing,size_t i, size_t j, size_t k,  size_t l)  {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);
	//assert(l<=len); // don't know len here

	// note: enclosed bp type 'turned around' for lib call
	const int ptype_enclosed = rtype[pair[S[k]][S[l]]];

	if (ptype_enclosed==0) return INF;

	return E_IntLoop(k-i-1,j-l-1,ptype_closing,ptype_enclosed,S1[i+1],S1[j-1],S1[k-1],S1[l+1],const_cast<paramT *>(params));
}


/**
* @brief Register a candidate
* @param i start
* @param j end
* @param e energy of candidate "V(i,j)"
*/
void register_candidate(auto &CL, size_t i, size_t j, energy_t e) {
	assert(i<=j+TURN+1);
	CL[j].push_back( cand_entry_t(i, e) );
}

std::pair< energy_t, energy_t > split_cases( auto const& CL, auto const& WM, auto const& S, auto const& params, int i, int j, auto &km1, int n, int *p_table ) {
	energy_t wm_split = INF;
	energy_t wm2_split = INF;
	int mm3 = S[j-1];

	for ( auto const [key,val] : CL) {
		size_t k = key;
		int mm5 = S[k+1];
		bool paired = (p_table[k] == j && p_table[j] == k);
		energy_t v_kj = E_MLStem(val,WM,CL,S,params,k,j,n,dangles,p_table);
		wm_split = std::min( wm_split, WM[k-1] + v_kj );
		bool can_pair = true;
		// checks to see if the unpaired bases till k can happen
		for(int m = i;m<k;++m){
			if(p_table[m]>-1) {
				can_pair = false;
				break;
			}
		}
		if(can_pair) wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );
		wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
		if(wm2_split==WM[k-1] + v_kj) km1 = k-1;
		
		if(paired) return std::make_pair( wm_split, wm2_split );
	}
	return std::make_pair( wm_split, wm2_split );

}
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
bool evaluate_restriction(int i, int j, int *p_table, int *last_j_array, int *in_pair_array, bool multiloop){
	bool evaluate = 1;
	if(in_pair_array[i]>in_pair_array[j]) evaluate = 0;

	if(in_pair_array[i]<in_pair_array[j]) evaluate = 0;

	if(in_pair_array[i]==in_pair_array[j]){
		if(j>last_j_array[i]) evaluate = 0;
	}
	// Resolves the cases where k-1 is the end of a restricted pair but i is less than the beginning of the k-1 pair
	// And where i is the beginning of the restricted pair but k-1 is past the end of the pair 
	if(multiloop){
		if((p_table[j] >0 && i<p_table[j] && j>p_table[j]) || (p_table[i]>0 && j > p_table[i] && i<p_table[i])) evaluate = 1;
	}
	return evaluate;
}

energy_t fold(auto const& seq, auto &V, auto const& cand_comp, auto &CL, auto const& S, auto const& S1, auto const& params, auto &ta, auto &W, auto &WM, auto &WM2, auto &dmli1, auto &dmli2, auto const& n, auto const& dangles, auto const& garbage_collect, int *p_table, int *last_j_array, int *in_pair_array) {
	
	for (size_t i=n; i>0; --i) {
		int si1 = (i>1) ? S[i-1] : -1;
		for ( size_t j=i+TURN+1; j<=n; j++ ) {

			int sj1 = (j<n) ? S[j+1] : -1;
			int mm5 = S[i+1];
			int mm3 = S[j-1];
			bool evaluate = evaluate_restriction(i,j,p_table,last_j_array,in_pair_array,false);
			// ------------------------------
			// W: split case
			bool pairedkj = 0;
			energy_t w_split = INF;
			for ( auto const [key,val] : CL[j] ) {
				size_t k=key;
				int sk1 = (k>1) ? S[k-1] : -1;
				bool unpairedkj = (p_table[k]<-1 && p_table[j]<-1);
				pairedkj = (p_table[k] == j && p_table[j] == k);
				energy_t v_kj = (unpairedkj || pairedkj) ? val + vrna_E_ext_stem(pair[S[k]][S[j]],sk1,sj1,params) : INF;
				if(pairedkj){
					w_split = W[k-1] + v_kj; 
					break;
				}else{
					w_split = std::min( w_split, W[k-1] + v_kj );
				}
			}
			if(p_table[j]<-1) w_split = std::min(w_split,W[j-1]);

			// ------------------------------
			// WM and WM2: split cases
			int km1 = n;
			auto [wm_split, wm2_split] = split_cases( CL[j], WM,S, params,i,j,km1,n,p_table);
			

			if(p_table[j]<-1) wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			if(p_table[j]<-1) wm_split = std::min( wm_split, WM[j-1] + params->MLbase );
			
			// Check to see if wm and wm2 can be split
			bool check = !(evaluate_restriction(i,km1,p_table,last_j_array,in_pair_array,true));
			if(check && km1 != n) wm2_split=wm_split=INF;
			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V



			size_t i_mod=i%(MAXLOOP+1);

			const int ptype_closing = pair[S[i]][S[j]];
			const bool restricted = p_table[i] == -1 || p_table[j] == -1;
			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0 && !restricted && evaluate) { // if i,j form a canonical base pair

				bool canH = true;
				if((p_table[i]>-1 && p_table[i] != j) || (p_table[j]>-1 && p_table[j] != i)) canH = false;
				for(int k=i+1;k<j;k++) if(p_table[k]>-1){canH = false;} // make more efficient later
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
				#pragma omp parallel for num_threads(6);
				for ( size_t k=i+1; k<=max_k; k++) {
					size_t k_mod=k%(MAXLOOP+1);

					size_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;

					for (size_t l=min_l; l<j; l++) {
						bool canI = true;

						
						for(int m = i+1; m<k;m++) if(p_table[m]>-1){canI = false;}
						for(int m = l+1; m<j;m++) if(p_table[m]>-1){canI = false;}
						if((p_table[i]>-1 && p_table[i] != j) || (p_table[j]>-1 && p_table[j] != i) || (p_table[k]>-1 && p_table[k] != l)) canI=false;

						assert(k-i+j-l-2<=MAXLOOP);

						const energy_t v_iloop_kl = canI ? V(k_mod,l) + ILoopE(S,S1,params,ptype_closing,i,j,k,l) : INF;
						if ( v_iloop_kl < v_iloop ) {
							v_iloop = v_iloop_kl;
							best_l=l;
							best_k=k;
							best_e=V(k_mod,l);
						}
					}
				}
				bool unpaired = (p_table[i]<-1 && p_table[j]<-1);
				bool paired = (p_table[i] == j && p_table[j] == i);
				const energy_t v_split = E_MbLoop(dmli1,dmli2,S,params,i,j,dangles,p_table);

				const energy_t v = std::min(v_h,std::min(v_iloop,v_split));

				const energy_t w_v  = (unpaired || paired) ? v + vrna_E_ext_stem(ptype_closing,si1,sj1,params): INF;
				const energy_t wm_v = (unpaired || paired) ? E_MLStem(v,WM,CL,S,params,i,j,n,dangles,p_table): INF;
				
				// update w and wm by v
				if(paired){
					w = w_v;
					wm = wm_v;
				} else if(pairedkj){
					w = w_split;
					wm = wm_split;
				} else{
					w  = std::min(w_v, w_split);
					wm = std::min(wm_v, wm_split);
				}
				
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
				// check whether (i,j) is a candidate; then register
				if ( w_v < w_split || wm_v < wm_split || paired) {
			
					register_candidate(CL, i, j, v );

					// always keep arrows starting from candidates
					inc_source_ref_count(ta,i,j);
				}
				V(i_mod,j) = v;
			} else {
				V(i_mod,j) = INF;
			} // end if (i,j form a canonical base pair)
			W[j]       = w;
			WM[j]      = wm;
			WM2[j]     = wm2_split;
			
			
			

		} // end loop j
		rotate_arrays(WM,WM2,dmli1,dmli2,n);
		// Clean up trace arrows in i+MAXLOOP+1
		if (garbage_collect && i+MAXLOOP+1 <= n) {
			gc_row(ta,i + MAXLOOP + 1 );
		}

		// Reallocate candidate lists in i
		for ( auto &x: CL ) {
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
void detect_restricted_pairs(auto const &structure, int *p_table, int *last_j_array, int *in_pair_array){
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

	std::string restricted;
    args_info.input_structure_given ? restricted = input_structure : restricted = std::string (seq.length(),'.');

	int dangle;
    args_info.dangles_given ? dangle = dangles : dangles = 2;
	if(restricted != "" && restricted.length() != seq.length() ){
		std::cout << "input sequence and structure are not the same size" << std::endl;
		exit(0);
	}

	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;

	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given,restricted, dangles);

	// Make replicate mx array in linear space
	int last_j_array[seq.length()+1] = {0};
	int in_pair_array[seq.length()+1] = {0};
	int p_table[seq.length()+1] = {0};
	detect_restricted_pairs(restricted,p_table,last_j_array,in_pair_array);

	std::cout << seq << std::endl;
	
	energy_t mfe = fold(sparsemfefold.seq_,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_, sparsemfefold.dmli1_, sparsemfefold.dmli2_,sparsemfefold.n_, sparsemfefold.dangles_,sparsemfefold.garbage_collect_, p_table,last_j_array,in_pair_array);	

	std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_,p_table,last_j_array,in_pair_array,dangles, mark_candidates);
	std::ostringstream smfe;
	smfe << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

	std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

	// size_t n=seq.length();

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
