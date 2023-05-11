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
#include <cassert>
#include <numeric>

#include "base.hh"
#include "trace_arrow.hh"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loops/all.h"
}

#include "cmdline.hh"





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


energy_t ILoopE(auto const& S_, auto const& S1_, auto const& params_, int ptype_closing,size_t i, size_t j, size_t k,  size_t l);
void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e );
void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j);
void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e) ;
void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j);
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
	

	bool garbage_collect_;

	LocARNA::Matrix<energy_t> V_; // store V[i..i+MAXLOOP-1][1..n]
	std::vector<energy_t> W_;
	std::vector<energy_t> WM_;
	std::vector<energy_t> WM2_;

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

	


	SparseMFEFold(const std::string &seq, bool garbage_collect)
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

	// init candidate lists
	CL_.resize(n_+1);

	resize(ta_,n_+1);
	}

	

	~SparseMFEFold() {
	free(params_);
	free(S_);
	free(S1_);
	}
};







// ! TRANSLATED: -----------------------------------------------------------------------------------

	energy_t HairpinE(auto const& seq, auto const &S, auto const &S1, auto const& params, size_t i, size_t j) {

	assert(1<=i);
	assert(i<j);
	
	//assert(j<=len); // don't know len here

	const int ptype_closing = pair[S[i]][S[j]];

	if (ptype_closing==0) return INF;

	return E_Hairpin(j-i-1,ptype_closing,S1[i+1],S1[j-1],&seq.c_str()[i-1], const_cast<paramT *>(params));
	}



	/**
	* @brief Recompute row of W
	*
	* @param i row index
	* @param max_j maximum column index
	*/
auto const recompute_W(auto const &W, auto const& CL, size_t i, size_t max_j) {
	//std::cout << "Compute W " <<i<<" "<<max_j<<std::endl;
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
		w = std::min(w,temp[j-1]);

		temp[j] = w;
	}
	return temp;
}



	/**
	* @brief Recompute row of WM
	*
	* @param i row index
	* @param max_j maximum column index
	*/
	auto const recompute_WM(auto const& WM, auto const &CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j) {
	//std::cout << "Compute WM " <<i<<" "<<max_j<<std::endl;

	assert(i>=1);
	assert(max_j<=n);

	std::vector<energy_t> temp = WM;

	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { temp[j]=INF; }

	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
		energy_t wm = INF;

		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>=i ; ++it ) {
			size_t k = it->first;
			const energy_t v_kj = it->second + E_MLstem(pair[S[k]][S[j]],-1,-1,params);
			wm = std::min( wm, static_cast<energy_t>(params->MLbase*(k-i)) + v_kj );
			wm = std::min( wm, temp[k-1]  + v_kj );
		}
		wm = std::min(wm, temp[j-1] + params->MLbase);

		temp[j] = wm;
	}
	return temp;
}

	/**
	* @brief Recompute row of WM2
	*
	* @param i row index
	* @param max_j maximum column index
	*/
	auto const recompute_WM2(auto const& WM, auto const& WM2, auto const CL, auto const& S, auto const &params, auto const& n, size_t i, size_t max_j) {
	//std::cout << "Recompute WM2 " <<i<<" "<<max_j<<std::endl;

	assert(i>=1);
	//assert(i+2*TURN+3<=max_j);
	assert(max_j<= n);

	std::vector<energy_t> temp = WM2;

	for ( size_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { temp[j]=INF; }

	for ( size_t j=i+2*TURN+3; j<=max_j; j++ ) {
		energy_t wm2 = INF;

		for ( auto it = CL[j].begin();CL[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
			size_t k = it->first;
			energy_t v_kl= it->second + E_MLstem(pair[S[k]][S[j]],-1,-1,params);
			wm2 = std::min( wm2, WM[k-1]  + v_kl );
		}
		wm2 = std::min(wm2, temp[j-1] + params->MLbase);

		temp[j] = wm2;
	}
	return temp;
}

/**
	* Test existence of candidate
	*
	* @param i start
	* @param j end
	*
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
	* @param i row index
	* @param j column index
	* pre: W contains values of row i in interval i..j
	*/
	void trace_W(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j) {
	// std::cout << "Trace W "<<i<<" "<<j<<std::endl;
	if (i+TURN+1>=j) return;

	// case j unpaired
	if (W[j] == W[j-1]) {
		trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,j-1);
		return;
	}

	size_t k=j+1;
	energy_t v=INF;

	// determine best split W -> W V
	std::any_of(CL[j].begin(),CL[j].end(),[&]( auto const key_val ) {
		auto const [key,val] = key_val;
		k = key>=i ? key : k;
		
		const energy_t v_kj = val + E_ExtLoop(pair[S[k]][S[j]],-1,-1,params);
		const energy_t w = W[k-1] + v_kj;

		v = W[j] == w && key>=i ? val : v;
		return key<i || W[j] == w;
	});
	// auto it = std::lower_bound(CL[j].begin(),CL[j].end(),W[j],[&]( auto const key_val, auto const comp) {
	// 	auto const [key,val] = key_val;
		
	// 	if(key<i) return false;
	// 	k = key;
	// 	const energy_t v_kj = val + E_ExtLoop(pair[S[k]][S[j]],-1,-1,params);
	// 	const energy_t w = W[k-1] + v_kj;
		
	// 	// v = W[j] == w && key>=i ? val : v;
	// 	return W[j] == comp;
	// });

	// if(*it != *(CL[j].end())){
	// 	auto const [key,val] = *it;
	// 	v=val;
	// 	k=key;
	// }


	assert(i<=k && k<j);
	assert(v<INF);

	// don't recompute W, since i is not changed
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,i,k-1);
	trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,v);
}

	/**
	* @brief Trace from V entry
	*
	* @param i row index
	* @param j column index
	* @param energy energy in V[i,j]
	*
	* pre: structure is string of size (n+1)
	*/
	void trace_V(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates, size_t i, size_t j, energy_t e ) {
	// std::cout << "trace_V "<<i<<" "<<j<<std::endl;

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
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l, arrow.target_energy());
		return;

	} else {

		assert(ptype_closing>0);

		// try to trace back to a candidate: (still) interior loop case
		for ( size_t l=i; l<j; l++) {
		for ( auto it=CL[l].begin(); CL[l].end()!=it && it->first>i; ++it ) {
			const size_t k=it->first;
			if (  e == it->second + ILoopE(S,S1,params,ptype_closing,i,j,k,l) ) {
			trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,l,it->second);
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
	auto const temp = recompute_WM(WM,CL,S,params,n,i+1,j-1);
	WM = temp;
	auto const temp2 = recompute_WM2(WM,WM2,CL,S,params,n,i+1,j-1);
	WM2 = temp2;
	
	trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i+1,j-1);
}
/**
	* @brief Trace from WM
	*
	* @param i row index
	* @param j column index
	* pre: vector WM is recomputed for row i
	*/
	void trace_WM(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j, energy_t e) {
	if (i+TURN+1>j) {return;}

	if ( e == WM[j-1] + params->MLbase ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1,WM[j-1]);
		return;
	}

	for ( auto it=CL[j].begin();CL[j].end() != it && it->first>=i;++it ) {
		const size_t k = it->first;
		const energy_t v_kj = it->second + E_MLstem(pair[S[k]][S[j]],-1,-1,params);
		if ( e == WM[k-1] + v_kj ) {
		// no recomp, same i
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1]);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second);
		return;
		} else if ( e == static_cast<energy_t>((k-i)*params->MLbase) + v_kj ) {
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second);
		return;
		}
	}
	assert(false);
}


/**
	* @brief Trace from WM2
	*
	* @param i row index
	* @param j column index
	* pre: vectors WM and WM2 are recomputed for row i
	* auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& S, auto const& S1, auto const& ta, auto &WM, auto &WM2
	*/ 
	void trace_WM2(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto &WM, auto &WM2, auto const& n, auto const& mark_candidates,size_t i, size_t j) {
	if (i+2*TURN+3>j) {return;}

	const energy_t e = WM2[j];

	// case j unpaired
	if ( e == WM2[j-1] + params->MLbase ) {
		// same i, no recomputation
		trace_WM2(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,j-1);
		return;
	}

	for ( auto it=CL[j].begin();CL[j].end() != it  && it->first>=i+TURN+1;++it ) {
		size_t k = it->first;
		const energy_t v_kj = it->second + E_MLstem(pair[S[k]][S[j]],-1,-1,params);
		if ( e == WM[k-1] + v_kj ) {
		trace_WM(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,i,k-1,WM[k-1]);
		trace_V(seq,CL,cand_comp,structure,params,S,S1,ta,WM,WM2,n,mark_candidates,k,j,it->second);
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
	const std::string & trace_back(auto const& seq, auto const& CL, auto const& cand_comp, auto &structure, auto const& params, auto const& S, auto const& S1, auto &ta, auto const& W, auto &WM, auto &WM2, auto const& n,auto const& mark_candidates=false) {

	structure.resize(n+1,'.');

	/* Traceback */
	trace_W(seq,CL,cand_comp,structure,params,S,S1,ta,W,WM,WM2,n,mark_candidates,1,n);
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

std::pair< energy_t, energy_t > split_cases( auto const& CL, auto const& WM, auto const& S, auto const& params, int i, int j ) {
	energy_t wm_split = INF;
	energy_t wm2_split = INF;

	unrolled::for_each( std::cbegin( CL[j] ), std::cend( CL[j] ), [&]( auto const key_val ) {
		auto const [key,val] = key_val;
		size_t k = key;
		energy_t v_kj = val + E_MLstem(pair[S[k]][S[j]],-1,-1,params);

		wm_split = std::min( wm_split, WM[k-1] + v_kj );
		wm_split = std::min( wm_split,static_cast<energy_t>((k-i)*params->MLbase) + v_kj );

		wm2_split = std::min( wm2_split, WM[k-1] + v_kj );
		});
	return std::make_pair( wm_split, wm2_split );

}

struct w_split_cost {
	size_t j, i;
	short * S;

	const std::vector<energy_t>& WM_;

	paramT * params;
	energy_t operator () ( energy_t prev_best, auto const& key_val ) const {
		size_t k = key_val.first ;

		energy_t v_kj = key_val.second + E_MLstem(pair[S[k]][S[j]],-1,-1,params);
		return std::min( prev_best, 
					std::min(energy_t(v_kj + WM_[k - 1])
							, v_kj + static_cast<energy_t>((k-i)*params->MLbase) ) );
	}
};

struct w_split_cost_2 {
	size_t j, i;
	short * S;
	const std::vector<energy_t>& WM_;
	paramT * params;
	energy_t operator () ( energy_t prev_best, auto const& key_val ) const {
		size_t k = key_val.first ;

		energy_t v_kj = key_val.second + E_MLstem(pair[S[k]][S[j]],-1,-1,params);
		return std::min( prev_best, 
					energy_t(v_kj + WM_[k - 1]) );
	}
};


// typedef unsigned short int cand_pos_t;
// typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
// typedef std::vector< cand_entry_t > cand_list_t;

std::pair< energy_t, energy_t > split_cases_1( auto const& CL, auto const& WM, auto const& S, auto const& params, size_t i, size_t j ) {
	energy_t wm_split = INF;
	energy_t wm2_split = INF;
	// w_split_cost{j,i, S, WM,params};

	wm_split = std::accumulate( 
						CL.begin()
						, CL.end()
						, wm_split
						, w_split_cost{j,i, S, WM,params} );
	wm2_split = std::accumulate( 
						CL.begin()
						, CL.end()
						, wm_split
						, w_split_cost_2{j,i, S,WM,params} );


	return std::make_pair( wm_split, wm2_split );

}



energy_t fold(auto const& seq, auto &V, auto const& cand_comp, auto &CL, auto const& S, auto const& S1, auto const& params, auto &ta, auto &W, auto &WM, auto &WM2, auto const& n, auto const& garbage_collect) {
	for (size_t i=n; i>0; --i) {
		energy_t WM2_ip1_jm1 = INF;
		for ( size_t j=i+TURN+1; j<=n; j++ ) {

			// ------------------------------
			// W: split case
			energy_t w_split = INF;
			unrolled::for_each( std::cbegin( CL[j] ), std::cend( CL[j] ), [&]( auto const key_val ) {
                      auto const [key,val] = key_val;
                      size_t k=key;
                      energy_t v_kj = val + E_ExtLoop(pair[S[k]][S[j]],-1,-1,params);
                      w_split = std::min( w_split, W[k-1] + v_kj ); },1);
			w_split = std::min(w_split,W[j-1]);

			// ------------------------------
			// WM and WM2: split cases
			auto [wm_split, wm2_split] = split_cases_1( CL[j], WM,S, params,i,j);

			wm2_split = std::min( wm2_split, WM2[j-1] + params->MLbase );
			wm_split = std::min( wm_split, WM[j-1] + params->MLbase );

			energy_t w  = w_split; // entry of W w/o contribution of V
			energy_t wm = wm_split; // entry of WM w/o contribution of V


			size_t i_mod=i%(MAXLOOP+1);

			const int ptype_closing = pair[S[i]][S[j]];

			// ----------------------------------------
			// cases with base pair (i,j)
			if(ptype_closing>0) { // if i,j form a canonical base pair

				energy_t v_h = HairpinE(seq,S,S1,params,i,j);

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
				for ( size_t k=i+1; k<=max_k; k++) {
					size_t k_mod=k%(MAXLOOP+1);

					size_t min_l=std::max(k+TURN+1 + MAXLOOP+2, k+j-i) - MAXLOOP-2;

					for (size_t l=min_l; l<j; l++) {

						assert(k-i+j-l-2<=MAXLOOP);

						const energy_t v_iloop_kl = V(k_mod,l) + ILoopE(S,S1,params,ptype_closing,i,j,k,l);

						if ( v_iloop_kl < v_iloop ) {
							v_iloop = v_iloop_kl;
							best_l=l;
							best_k=k;
							best_e=V(k_mod,l);
						}
					}
				}

				const energy_t v_split = WM2_ip1_jm1 + E_MLstem(rtype[ptype_closing],-1,-1,params) + params->MLclosing;
				// this value, conceptually
				// WM2(i+1,j-1), is set in the
				// previous j-iteration before this
				// value in array WM2_[] is overwritten

				const energy_t v = std::min(v_h,std::min(v_iloop,v_split));

				const energy_t w_v  = v + E_ExtLoop(ptype_closing,-1,-1,params);
				const energy_t wm_v = v + E_MLstem(ptype_closing,-1,-1,params);

				// update w and wm by v
				w  = std::min(w_v, w_split);
				wm = std::min(wm_v, wm_split);

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
				if ( w_v < w_split || wm_v < wm_split ) {

					//std::cout << "Reg Cand "<<i<<","<<j<<std::endl;

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

			WM2_ip1_jm1 = WM2[j]; // here, the array WM2_ still
					// contains WM2(i+1,j-1); in
					// the next j-iteration, we
					// need this.
			WM2[j]     = wm2_split;

		} // end loop j

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


size_t num_of_candidates(auto const& CL_)  {
	size_t c=0;
	for ( auto const &x: CL_ ) {
		c += x.size();
	}
	return c;
}

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


	bool verbose;
	verbose = args_info.verbose_given;

	bool mark_candidates;
	mark_candidates = args_info.mark_candidates_given;

	SparseMFEFold sparsemfefold(seq,!args_info.noGC_given);


	cmdline_parser_free(&args_info);

	std::cout << seq << std::endl;
	//std::cout << "Len:\t"<<seq.length()<<std::endl<<std::endl;

	energy_t mfe = fold(sparsemfefold.seq_,sparsemfefold.V_,sparsemfefold.cand_comp,sparsemfefold.CL_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.params_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_,sparsemfefold.garbage_collect_);

	std::string structure = trace_back(sparsemfefold.seq_,sparsemfefold.CL_,sparsemfefold.cand_comp,sparsemfefold.structure_,sparsemfefold.params_,sparsemfefold.S_,sparsemfefold.S1_,sparsemfefold.ta_,sparsemfefold.W_,sparsemfefold.WM_,sparsemfefold.WM2_,sparsemfefold.n_, mark_candidates);

	std::ostringstream smfe;
	smfe << std::setw(6) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

	std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

	size_t n=seq.length();

	float factor=1024;
	
	const std::string unit=" kB";

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
