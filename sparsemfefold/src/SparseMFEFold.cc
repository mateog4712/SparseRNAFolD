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

#include "base.hh"
#include "trace_arrow.hh"

extern "C" {
#include "ViennaRNA/pair_mat.h"
#include "ViennaRNA/loop_energies.h"
}

#include "SparseMFEFold_cmdline.h"





typedef unsigned short int cand_pos_t;
typedef std::pair<cand_pos_t,energy_t> cand_entry_t;
typedef std::vector< cand_entry_t > cand_list_t;

class SparseMFEFold;

int pair_type(SparseMFEFold & fold, size_t i, size_t j);
energy_t ILoopE(SparseMFEFold &fold, int ptype_closing,size_t i, size_t j, size_t k,  size_t l);
void trace_V( SparseMFEFold &fold, size_t i, size_t j, energy_t e );
void trace_W(SparseMFEFold &fold, size_t i, size_t j);
void trace_WM( SparseMFEFold &fold,size_t i, size_t j, energy_t e) ;
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

energy_t HairpinE(SparseMFEFold &fold,size_t i, size_t j) {

	assert(1<=i);
	assert(i<j);
	//assert(j<=len); // don't know len here

	int ptype_closing = pair_type(fold, i,j);

	if (ptype_closing==0) return INF;

	return
	    E_Hairpin(j-i-1,
		      ptype_closing,
		      fold.S1_[i+1],
		      fold.S1_[j-1],
		      &fold.seq_.c_str()[i-1],
		      const_cast<paramT *>(fold.params_));

    }
/**
     * @brief Recompute row of W
     *
     * @param i row index
     * @param max_j maximum column index
     */
    void recompute_W(SparseMFEFold &fold, size_t i, size_t max_j) {
	//std::cout << "Compute W " <<i<<" "<<max_j<<std::endl;

	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { fold.W_[j]=0; }
	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {

	    energy_t w = INF;

	    // note: the loop covers the case W(i,j)=V(i,j),
	    // since this is in the candidate list (TBS)
	    for ( auto it = fold.CL_[j].begin();
		  fold.CL_[j].end()!=it && it->first>=i ; ++it ) {
		w = std::min( w, fold.W_[it->first-1] + it->second );
	    }
	    // case "j unpaired" is not in the CL (anymore)
	    w = std::min(w,fold.W_[j-1]);

	    fold.W_[j] = w;
	}
}
/**
     * @brief Recompute row of WM
     *
     * @param i row index
     * @param max_j maximum column index
     */
    void recompute_WM(SparseMFEFold &fold,size_t i, size_t max_j) {
	//std::cout << "Compute WM " <<i<<" "<<max_j<<std::endl;

	assert(i>=1);
	assert(max_j<=fold.n_);

	for ( size_t j=i-1; j<=std::min(i+TURN,max_j); j++ ) { fold.WM_[j]=INF; }

	for ( size_t j=i+TURN+1; j<=max_j; j++ ) {
	    energy_t wm = INF;

	    for ( auto it = fold.CL_[j].begin();
		  fold.CL_[j].end()!=it && it->first>=i ; ++it ) {
		size_t k = it->first;
		energy_t v_kj =
		    it->second + E_MLstem(pair_type(fold, k,j),-1,-1,fold.params_);
		wm = std::min( wm, static_cast<energy_t>(fold.params_->MLbase*(k-i)) + v_kj );
		wm = std::min( wm, fold.WM_[k-1]  + v_kj );
	    }
	    wm = std::min(wm, fold.WM_[j-1] + fold.params_->MLbase);

	    fold.WM_[j] = wm;
	}
}

/**
     * @brief Recompute row of WM2
     *
     * @param i row index
     * @param max_j maximum column index
     */
    void recompute_WM2(SparseMFEFold &fold, size_t i, size_t max_j) {
	//std::cout << "Recompute WM2 " <<i<<" "<<max_j<<std::endl;

	assert(i>=1);
	//assert(i+2*TURN+3<=max_j);
	assert(max_j<=fold.n_);

	for ( size_t j=i-1; j<=std::min(i+2*TURN+2,max_j); j++ ) { fold.WM2_[j]=INF; }

	for ( size_t j=i+2*TURN+3; j<=max_j; j++ ) {
	    energy_t wm2 = INF;

	    for ( auto it = fold.CL_[j].begin();
		  fold.CL_[j].end()!=it && it->first>i+TURN+1 ; ++it ) {
		size_t k = it->first;
		energy_t v_kl=
		    it->second + E_MLstem(pair_type(fold, k,j),-1,-1,fold.params_);
		wm2 = std::min( wm2, fold.WM_[k-1]  + v_kl );
	    }
	    wm2 = std::min(wm2, fold.WM2_[j-1] + fold.params_->MLbase);

	    fold.WM2_[j] = wm2;
	}
}

/**
     * Test existence of candidate
     *
     * @param i start
     * @param j end
     *
     * @return whether (i,j) is candidate for W/WM splits
     */
    bool is_candidate(SparseMFEFold &fold,size_t i, size_t j) {
	const cand_list_t &list = fold.CL_[j];

	auto it = std::lower_bound(list.begin(),list.end(),i,fold.cand_comp);

	return it!=list.end() && it->first==i;
    }
/**
     * @brief Trace from WM2
     *
     * @param i row index
     * @param j column index
     * pre: vectors WM and WM2 are recomputed for row i
     */
    void trace_WM2(SparseMFEFold &fold,size_t i, size_t j) {
	if (i+2*TURN+3>j) {return;}

	energy_t e = fold.WM2_[j];

	// case j unpaired
	if ( e == fold.WM2_[j-1] + fold.params_->MLbase ) {
	    // same i, no recomputation
	    trace_WM2(fold,i,j-1);
	    return;
	}

	for ( auto it=fold.CL_[j].begin();
	      fold.CL_[j].end() != it  && it->first>=i+TURN+1;
	      ++it ) {
	    size_t k = it->first;
	    energy_t v_kj = it->second + E_MLstem(pair_type(fold,k,j),-1,-1,fold.params_);
	    if ( e == fold.WM_[k-1] + v_kj ) {
		trace_WM(fold,i,k-1,fold.WM_[k-1]);
		trace_V(fold,k,j,it->second);
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
    const std::string & trace_back(SparseMFEFold &fold,bool mark_candidates=false) {
	fold.mark_candidates_=mark_candidates;

	fold.structure_.resize(fold.n_+1,'.');

	/* Traceback */
	trace_W(fold,1,fold.n_);
	fold.structure_ = fold.structure_.substr(1,fold.n_);

	return fold.structure_;
}

/**
     * @brief Trace from WM
     *
     * @param i row index
     * @param j column index
     * pre: vector WM is recomputed for row i
     */
    void trace_WM( SparseMFEFold &fold,size_t i, size_t j, energy_t e) {
	if (i+TURN+1>j) {return;}

	if ( e == fold.WM_[j-1] + fold.params_->MLbase ) {
	    trace_WM(fold,i,j-1,fold.WM_[j-1]);
	    return;
	}

	for ( auto it=fold.CL_[j].begin();
	      fold.CL_[j].end() != it && it->first>=i;
	      ++it ) {
	    size_t k = it->first;
	    energy_t v_kj = it->second
		+ E_MLstem(pair_type(fold,k,j),-1,-1,fold.params_);
	    if ( e == fold.WM_[k-1] + v_kj ) {
		// no recomp, same i
		trace_WM(fold,i,k-1,fold.WM_[k-1]);
		trace_V(fold,k,j,it->second);
		return;
	    } else if ( e == static_cast<energy_t>((k-i)*fold.params_->MLbase) + v_kj ) {
		trace_V(fold,k,j,it->second);
		return;
	    }
	}
	assert(false);
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
    void trace_V( SparseMFEFold &fold, size_t i, size_t j, energy_t e ) {
	// std::cout << "trace_V "<<i<<" "<<j<<std::endl;

	assert( i+TURN+1<=j );
	assert( j<=n_ );

	if (fold.mark_candidates_ && is_candidate(fold,i,j)) {
	    fold.structure_[i]='[';
	    fold.structure_[j]=']';
	} else {
	    fold.structure_[i]='(';
	    fold.structure_[j]=')';
	}

	int ptype_closing = pair_type(fold,i,j);

	if (exists_trace_arrow_from(fold.ta_,i,j)) {
	    // trace arrows may exist for interior loop case
	    const TraceArrow &arrow = trace_arrow_from(fold.ta_,i,j);

	    size_t k=arrow.k(i,j);
	    size_t l=arrow.l(i,j);
	    assert(i<k);
	    assert(l<j);
	    trace_V(fold,k,l, arrow.target_energy());
	    return;

	} else {

	    assert(ptype_closing>0);

	    // try to trace back to a candidate: (still) interior loop case
	    for ( size_t l=i; l<j; l++) {
		for ( auto it=fold.CL_[l].begin(); fold.CL_[l].end()!=it && it->first>i; ++it ) {
		    size_t k=it->first;
		    if (  e == it->second + ILoopE(fold,ptype_closing,i,j,k,l) ) {
			trace_V(fold,k,l,it->second);
			return;
		    }
		}
	    }
	}

	// is this a hairpin?
	if ( e == HairpinE(fold,i,j) ) {
	    return;
	}

	// if we are still here, trace to wm2 (split case);
	// in this case, we know the 'trace arrow'; the next row has to be recomputed
	recompute_WM(fold,i+1,j-1);
	recompute_WM2(fold,i+1,j-1);
	trace_WM2(fold,i+1,j-1);
}
/**
     * @brief Trace from W entry
     *
     * @param i row index
     * @param j column index
     * pre: W contains values of row i in interval i..j
     */
    void trace_W(SparseMFEFold &fold, size_t i, size_t j) {
	// std::cout << "Trace W "<<i<<" "<<j<<std::endl;
	if (i+TURN+1>=j) return;

	// case j unpaired
	if (fold.W_[j] == fold.W_[j-1]) {
	    trace_W(fold,i,j-1);
	    return;
	}

	size_t k=j+1;
	energy_t v=INF;

	// determine best split W -> W V
	for ( auto it = fold.CL_[j].begin();
	      fold.CL_[j].end()!=it && it->first>=i;
	      ++it ) {
	    k = it->first;
	    energy_t v_kj = it->second + E_ExtLoop(pair_type(fold, k,j),-1,-1,fold.params_);
	    energy_t w = fold.W_[k-1] + v_kj;

	    if (fold.W_[j] == w) {
		v = it->second;
		break;
	    }
	}

	assert(i<=k && k<j);
	assert(v<INF);

	// don't recompute W, since i is not changed
	trace_W(fold,i,k-1);
	trace_V(fold,k,j,v);
}

int pair_type(SparseMFEFold &fold, size_t i, size_t j)  {
	return pair[fold.S_[i]][fold.S_[j]];
    }
/* pre: ptype_closing>0 */
energy_t ILoopE(SparseMFEFold &fold, int ptype_closing,size_t i, size_t j, size_t k,  size_t l)  {
	assert(ptype_closing>0);
	assert(1<=i);
	assert(i<k);
	assert(k<l);
	assert(l<j);
	//assert(l<=len); // don't know len here

	// note: enclosed bp type 'turned around' for lib call
	int ptype_enclosed = rtype[pair_type(fold,k,l)];

	if (ptype_enclosed==0) return INF;

	return
	    E_IntLoop(k-i-1,j-l-1,
		      ptype_closing,
		      ptype_enclosed,
		      fold.S1_[i+1],
		      fold.S1_[j-1],
		      fold.S1_[k-1],
		      fold.S1_[l+1],
		      const_cast<paramT *>(fold.params_));
}


/**
     * @brief Register a candidate
     * @param i start
     * @param j end
     * @param e energy of candidate "V(i,j)"
     */
    void register_candidate(SparseMFEFold &fold, size_t i, size_t j, energy_t e) {
		assert(i<=j+TURN+1);
		fold.CL_[j].push_back( cand_entry_t(i, e) );
    }

energy_t fold(SparseMFEFold &fold) {
	for (size_t i=fold.n_; i>0; --i) {
	energy_t WM2_ip1_jm1 = INF;
	for ( size_t j=i+TURN+1; j<=fold.n_; j++ ) {

	// ------------------------------
	// W: split case
	energy_t w_split = INF;
	for ( auto &x : fold.CL_[j] ) {
		size_t k=x.first;
		energy_t v_kj =
		x.second
		+ E_ExtLoop(pair_type(fold, k,j),-1,-1,fold.params_);
		w_split = std::min( w_split, fold.W_[k-1] + v_kj );
	}
	w_split = std::min(w_split,fold.W_[j-1]);

	// ------------------------------
	// WM and WM2: split cases
	energy_t wm_split = INF;
	energy_t wm2_split = INF;
	for ( auto &x : fold.CL_[j] ) {
		size_t k = x.first;
		energy_t v_kj =
		x.second
		+ E_MLstem(pair_type(fold,k,j),-1,-1,fold.params_);

		wm_split = std::min( wm_split, fold.WM_[k-1] + v_kj );
		wm_split = std::min( wm_split,
					static_cast<energy_t>((k-i)*fold.params_->MLbase) + v_kj );

		wm2_split = std::min( wm2_split, fold.WM_[k-1] + v_kj );
	}

	wm2_split = std::min( wm2_split, fold.WM2_[j-1] + fold.params_->MLbase );
	wm_split = std::min( wm_split, fold.WM_[j-1] + fold.params_->MLbase );

	energy_t w  = w_split; // entry of W w/o contribution of V
	energy_t wm = wm_split; // entry of WM w/o contribution of V


	size_t i_mod=i%(MAXLOOP+1);

	int ptype_closing = pair_type(fold,i,j);

	// ----------------------------------------
	// cases with base pair (i,j)
	if(ptype_closing>0) { // if i,j form a canonical base pair

		energy_t v_h = HairpinE(fold,i,j);

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

			energy_t v_iloop_kl =
			fold.V_(k_mod,l)
			+ ILoopE(fold,ptype_closing,i,j,k,l);

			if ( v_iloop_kl < v_iloop ) {
			v_iloop = v_iloop_kl;
			best_l=l;
			best_k=k;
			best_e=fold.V_(k_mod,l);
			}
		}
		}

		energy_t v_split =
		WM2_ip1_jm1 // this value, conceptually
				// WM2(i+1,j-1), is set in the
				// previous j-iteration before this
				// value in array WM2_[] is overwritten
		+ E_MLstem(rtype[ptype_closing],-1,-1,fold.params_)
		+ fold.params_->MLclosing;

		energy_t v = std::min(v_h,std::min(v_iloop,v_split));

		energy_t w_v  = v + E_ExtLoop(ptype_closing,-1,-1,fold.params_);
		energy_t wm_v = v + E_MLstem(ptype_closing,-1,-1,fold.params_);

		// update w and wm by v
		w  = std::min(w_v, w_split);
		wm = std::min(wm_v, wm_split);

		// register required trace arrows from (i,j)
		if ( v_iloop < std::min(v_h,v_split) ) {
		if ( is_candidate(fold,best_k,best_l) ) {
			//std::cout << "Avoid TA "<<best_k<<" "<<best_l<<std::endl;
			avoid_trace_arrow(fold.ta_);
		} else {
			//std::cout<<"Reg TA "<<i<<","<<j<<":"<<best_k<<","<<best_l<<std::endl;
			register_trace_arrow(fold.ta_,i,j,best_k,best_l,best_e);
		}
		}

		// check whether (i,j) is a candidate; then register
		if ( w_v < w_split
			||
			wm_v < wm_split ) {

		//std::cout << "Reg Cand "<<i<<","<<j<<std::endl;

		register_candidate( fold, i, j, v );

		// always keep arrows starting from candidates
		inc_source_ref_count(fold.ta_,i,j);
		}

		fold.V_(i_mod,j) = v;

	} else {
		fold.V_(i_mod,j) = INF;
	} // end if (i,j form a canonical base pair)

	fold.W_[j]       = w;
	fold.WM_[j]      = wm;

	WM2_ip1_jm1 = fold.WM2_[j]; // here, the array WM2_ still
					// contains WM2(i+1,j-1); in
					// the next j-iteration, we
					// need this.
	fold.WM2_[j]     = wm2_split;

	} // end loop j

	// Clean up trace arrows in i+MAXLOOP+1
	if ( fold.garbage_collect_ && i+MAXLOOP+1 <= fold.n_) {
	gc_row(fold.ta_,i + MAXLOOP + 1 );
	}

	// Reallocate candidate lists in i
	for ( auto &x: fold.CL_ ) {
	if (x.capacity() > 1.5*x.size()) {
		cand_list_t vec(x.size());
		copy(x.begin(),x.end(),vec.begin());
		vec.swap(x);
	}
	}

	compactify(fold.ta_);
}

return fold.W_[fold.n_];
}





size_t num_of_candidates(SparseMFEFold &fold)  {
	size_t c=0;
	for ( auto &x: fold.CL_ ) {
		c += x.size();
	}
	return c;
}

TraceArrows & ta(SparseMFEFold &fold){
	return fold.ta_;
    }

size_t capacity_of_candidates(SparseMFEFold &fold) {
	size_t c=0;
	for ( auto &x: fold.CL_ ) {
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

    gengetopt_args_info args_info;

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

    energy_t mfe = fold(sparsemfefold);

    std::string structure = trace_back(sparsemfefold, mark_candidates);

    std::ostringstream smfe;
    smfe << std::setw(6) << std::setiosflags(std::ios::fixed) << std::setprecision(2) << mfe/100.0 ;

    std::cout << structure << " ("<<smfe.str()<<")"<<std::endl;

    size_t n=seq.length();

    float factor=1024;
    const std::string unit=" kB";

    if (verbose) {
		

	std::cout <<std::endl;

	std::cout << "TA cnt:\t"<<sizeT(ta(sparsemfefold))<<std::endl;
	std::cout << "TA max:\t"<<maxT(ta(sparsemfefold))<<std::endl;
	std::cout << "TA av:\t"<<avoidedT(ta(sparsemfefold))<<std::endl;
	std::cout << "TA rm:\t"<<erasedT(ta(sparsemfefold))<<std::endl;
	// std::cout << "TA rm:\t"<<ta(sparsemfefold).erased()<<std::endl;

	std::cout <<std::endl;
	std::cout << "Can num:\t"<<num_of_candidates(sparsemfefold)<<std::endl;
	std::cout << "Can cap:\t"<<capacity_of_candidates(sparsemfefold)<<std::endl;
	std::cout << "TAs num:\t"<<sizeT(ta(sparsemfefold))<<std::endl;
	// std::cout << "TAs cap:\t"<<ta(sparsemfefold).capacity()<<std::endl;
	std::cout << "TAs cap:\t"<<capacityT(ta(sparsemfefold))<<std::endl;
    }

    return 0;
}
