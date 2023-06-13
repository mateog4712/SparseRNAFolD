#include "trace_arrow.hh"

TraceArrows::TraceArrows(cand_pos_tu n)
    : n_(n),
      ta_count_(0),
      ta_avoid_(0),
      ta_erase_(0),
      ta_max_(0)
{}


TraceArrow & trace_arrow_from(TraceArrows &t, cand_pos_tu i, cand_pos_tu j) {
return t.trace_arrow_[i].find(j)->second;
}

bool exists_trace_arrow_from(TraceArrows &t,cand_pos_tu i, cand_pos_tu j){
return t.trace_arrow_[i].exists(j);
}


void avoid_trace_arrow(TraceArrows &t){
    t.ta_avoid_++;
}



void register_trace_arrow(TraceArrows &t,cand_pos_tu i, cand_pos_tu j,cand_pos_tu k, cand_pos_tu l,energy_t e) {
// std::cout << "register_trace_arrow "<<i<<" "<<j<<" "<<k<<" "<<l<<std::endl;
t.trace_arrow_[i].push_ascending( j, TraceArrow(i,j,k,l,e) );

inc_source_ref_count(t,k,l);

t.ta_count_++;
t.ta_max_ = std::max(t.ta_max_,t.ta_count_);
}


void inc_source_ref_count(TraceArrows &t, cand_pos_tu i, cand_pos_tu j) {
	// get trace arrow from (i,j) if it exists
	if (! t.trace_arrow_[i].exists(j)) return;

	auto it=t.trace_arrow_[i].find(j);

	TraceArrow &ta=it->second;

	ta.inc_src();
}



/**
 * @brief Resizes the trace arrows list to size n
 * 
 * @param t Trace Arrows list
 * @param n New size
 */
void resize(TraceArrows &t,cand_pos_tu n) {
    t.trace_arrow_.resize(n);
}


void gc_trace_arrow(TraceArrows &t, cand_pos_tu i, cand_pos_tu j) {

    assert( t.trace_arrow_[i].exists(j) );

    auto col = t.trace_arrow_[i].find(j);

    const auto &ta = col->second;
    
    
    if (ta.source_ref_count() == 0) {
        // get trace arrow from the target if the arrow exists
        if (exists_trace_arrow_from(t,ta.k(i,j),ta.l(i,j))) {
            auto &target_ta = trace_arrow_from(t,ta.k(i,j),ta.l(i,j));
            target_ta.dec_src();
            gc_trace_arrow(t,ta.k(i,j),ta.l(i,j));
        }
        
        t.trace_arrow_[i].erase(col);
        t.ta_count_--;
        t.ta_erase_++;
    }
}

void gc_row(TraceArrows &t, cand_pos_tu i ) {
    assert(i<=t.n_);
    // i + TURN + 1
    for (cand_pos_tu j=i+4; j<=t.n_ ; j++) {
	    if (! t.trace_arrow_[i].exists(j)) continue;
	    gc_trace_arrow(t,i,j);
    }
}

void compactify(TraceArrows &t) {
    for ( auto &x: t.trace_arrow_ ) {
	if (x.capacity() > 1.2 * x.size()) {
	    x.reallocate();
	}
    }
}

cand_pos_tu numberT(TraceArrows &t){
    size_t c=0;
    for ( auto &x: t.trace_arrow_ ) {
	c += x.size();
    }
    return c; 
}

cand_pos_tu capacityT(TraceArrows &t){
    size_t c=0;
    for ( auto &x: t.trace_arrow_ ) {
	c += x.capacity();
    }
    return c;
}
cand_pos_tu erasedT(TraceArrows &t){
    return t.ta_erase_;
}
cand_pos_tu sizeT(TraceArrows &t){
    return t.ta_count_;
}
cand_pos_tu avoidedT(TraceArrows &t){
    return t.ta_avoid_;
}
cand_pos_tu maxT(TraceArrows &t){
    return t.ta_max_;
}