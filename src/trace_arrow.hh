#ifndef TRACE_ARROW_HH
#define TRACE_ARROW_HH

#include "base_types.hh"
#include "simple_map.hh"
#define NDEBUG
#include <cassert>

/**
 * @brief Trace arrow
 *
 * Describes a trace arrow from i,j to k_,l_. The source (X,i,j)
 * is not represented in the data structure. However, each trace arrow
 * is associated with exactly one source.  Source and target matrix
 * types are omitted, since we don't need them here, but in more
 * general scenarios, such information has to be maintained.
 */
class TraceArrow {
    unsigned char k_; //!< offset to target row of arrow
    unsigned char l_; //!< offset to target column of arrow
    energy_t energy_; //!<target energy
    uint count_; //!< counts how many trace arrows point to the source
public:
    /**
     * @brief construct by target coordinates
     * @param i source row
     * @param j source column
     * @param k target row
     * @param l target column
     */
    TraceArrow(cand_pos_tu i,cand_pos_tu j,cand_pos_tu k,cand_pos_tu l,energy_t e)
	: k_(k-i),l_(j-l),energy_(e),count_(0)
    {}

    /**
     * @brief empty c'tor
     */
    TraceArrow() {}

    cand_pos_tu k(cand_pos_tu i,cand_pos_tu j) const {return k_+i;}
    cand_pos_tu l(cand_pos_tu i,cand_pos_tu j) const {return j-l_;}
    energy_t target_energy() const {return energy_;}
    cand_pos_tu source_ref_count() const {return count_;}

    void inc_src() {count_++;}
    void dec_src() {count_--;}

};

/**
 * @brief Collection of trace arrows
 *
 * Stores trace arrows to be accessible by row and col index.  Access
 * by column index is logarithmic. TAs of one row are
 * traversable. Supports garbage collection of TAs. Keeps track of
 * several statistics on TAs.
 */
class TraceArrows {

public:
    typedef SimpleMap< cand_pos_tu, TraceArrow >       trace_arrow_row_map_t;
    typedef std::vector< trace_arrow_row_map_t >  trace_arrow_map_t;
    trace_arrow_map_t trace_arrow_;
    cand_pos_tu n_; //!< sequence length

    cand_pos_tu ta_count_; // count all generated tas
    cand_pos_tu ta_avoid_; // count all avoided tas (since they point to candidates)
    cand_pos_tu ta_erase_; // count all erased tas (in gc)
    cand_pos_tu ta_max_; // keep track of maximum number of tas, existing simultaneously
public:

    /**
     * @brief Construct for sequence of specific length
     * @param n sequence length
     */
    TraceArrows(cand_pos_tu n);

    /**
     * @brief Clear the TA structure to be used again
     * 
     */
    void reset(){
        ta_count_ = 0;
        ta_avoid_ = 0;
        ta_erase_ = 0;
        ta_max_ = 0;
        trace_arrow_.clear();
    }

      
};


/**
* Get target of trace arrow by source (non-const)
*
* @param i source row index
* @param j source column index
*/
TraceArrow & trace_arrow_from(TraceArrows &t, cand_pos_tu i, cand_pos_tu j);

/**
* Check existence of trace arrow by source
*
* @param i source row index
* @param j source column index
* @returns whether trace arrow exists
*/
bool exists_trace_arrow_from(TraceArrows &t,cand_pos_tu i, cand_pos_tu j);



/**
* avoid one trace arrow (for statistics only)
*/
void avoid_trace_arrow(TraceArrows &t);

/**
 * Increment the reference count of the source
 *
 * @param i source row index
 * @param j source column index
 *
 * If no trace arrow from source exists, do nothing
 */
void inc_source_ref_count(TraceArrows &t, cand_pos_tu i, cand_pos_tu j);


/**
 * Register trace arrow
 *
 * @param srctype source matrix type
 * @param i source row
 * @param j source column
 * @param tgttype target matrix type
 * @param k target row
 * @param l target column
 */
void register_trace_arrow(TraceArrows &t,cand_pos_tu i, cand_pos_tu j,cand_pos_tu k, cand_pos_tu l,energy_t e);


void resize(TraceArrows &t,cand_pos_tu n);


/**
     * Garbage collect trace arrow
     *
     * if count = 0 then remove trace arrow; recursively decrement
     * targets and remove if count drops to 0
     */
void gc_trace_arrow(TraceArrows &t, cand_pos_tu i, cand_pos_tu j);

void gc_trace_arrow(TraceArrows &t, TraceArrows &td, cand_pos_tu i, cand_pos_tu j);

void gc_row(TraceArrows &t, cand_pos_tu i );

void gc_row(TraceArrows &t, TraceArrows &td, cand_pos_tu i );


/**
* @brief Compactify heap space
*/
void compactify(TraceArrows &t);



/** @brief Number of trace arrows
* @return number
*/
cand_pos_tu numberT(TraceArrows &t);

cand_pos_tu sizeT(TraceArrows &t);
cand_pos_tu erasedT(TraceArrows &t);
cand_pos_tu avoidedT(TraceArrows &t);
cand_pos_tu maxT(TraceArrows &t);

/** @brief Capacity of trace arrows vectors
* @return capacity
*/
cand_pos_tu capacityT(TraceArrows &t);



#endif // TRACE_ARROW_HH
