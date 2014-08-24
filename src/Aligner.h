#ifndef ALIGNER_H__
#define ALIGNER_H__
#include <fst/fstlib.h>
#include "utf8.h"
#include "DictionaryReader.h"
#include <sys/stat.h>
using namespace fst;

enum MaxType { JOINT_MAX,
	       GP_MAX,
	       PG_MAX };

template<class F>
class Aligner {
 public:
  typedef F FST;
  typedef typename F::Arc Arc;
  typedef typename Arc::StateId StateId;
  typedef typename Arc::Weight Weight;
  typedef typename Arc::Label Label;
  typedef unordered_map<Label, Weight> LWMap;
  typedef unordered_map<Label, LWMap> EMMap;
  typedef unordered_map<Label, int> PenaltyMap;

  Aligner (const DictionaryReader& reader, UTFMap& utfmap, 
	   int gmax, int pmax, bool delg, bool delp, 
	   const string mdelim, const string cdelim,
	   SymbolTable& isyms, SymbolTable& osyms, 
	   const string& far_name, MaxType max_type, int verbose = 0)
    : reader_(reader), utfmap_(utfmap), gmax_(gmax), pmax_(pmax), 
      delg_(delg), delp_(delp), mdelim_(mdelim), cdelim_(cdelim),
      isyms_(isyms), osyms_(osyms),
      far_name_(far_name), max_type_(max_type), verbose_(verbose) { 
    }

  void InitModel () {
    //Book-keeping stuff
    string key_prefix = "";
    string key_suffix = ".fst";
    string key        = "";
    char   keybuf[16];
    int32  generate_keys = 7; //Suitable for up to several million lattices
    string tmpfar_name = far_name_ + ".tmp.far";
    
    FarWriter<Arc>* tmpfar_writer = FarWriter<Arc>::Create (tmpfar_name, FAR_DEFAULT);

    //Write the unweighted, processed alignment lattices to disk
    for (int i = 0; i < reader_.entries.size (); i++) {
      const Entry& entry = reader_.entries [i];
      VectorFst<Arc> fst;
      if (verbose_ >= 2 && i % 10000 == 0)
	cerr << "Processed: " << i << endl;
      EntryToFst (entry, fst);
      if (fst.NumStates () == 0)
	continue;

      if (i == reader_.entries.size () - 1) {
	fst.SetInputSymbols (&isyms_);
	fst.SetOutputSymbols (&osyms_);
      }
      sprintf (keybuf, "%0*d", generate_keys, i + 1);
      key = keybuf;
      tmpfar_writer->Add (key_prefix + key + key_suffix, fst);
    }
    delete tmpfar_writer;

    //Next init the FarReader for EM weight training
    tmpfar_reader_ = FarReader<Arc>::Open (tmpfar_name);
    total_prev_ = Weight::One ();
    //Init the previous model with the current values
    emmap_ = emmap_prev_;

    MaximizationStep ();
  }

  double MaximizationStep () {
    double change = GetChange ();

    if (max_type_ == JOINT_MAX) {
      total_prev_ = total_;
      for (typename EMMap::iterator it = emmap_prev_.begin (); 
	   it != emmap_prev_.end (); ++it) {
	for (typename LWMap::iterator jt = emmap_prev_[it->first].begin ();
	     jt != emmap_prev_[it->first].end (); ++jt) {
	  emmap_[it->first][jt->first] = Divide (jt->second, total_);
	  jt->second = Weight::Zero ();
	}
      }
      total_ = Weight::Zero ();
    } else if (max_type_ == GP_MAX) {
    } else if (max_type_ == PG_MAX) {
    }
    return change;
  }

  void ExpectationStep () {
    if (!tmpfar_reader_)
      return;

    tmpfar_reader_->Reset ();
    vector<Weight> alpha, beta;
    for (int i = 1; !tmpfar_reader_->Done (); tmpfar_reader_->Next (), ++i) {
      if (verbose_ >= 2 && i % 10000 == 0)
	cerr << "E-step processed: " << i << endl;

      string key = tmpfar_reader_->GetKey ();
      VectorFst<Arc> fst (tmpfar_reader_->GetFst ());
      ArcMap (&fst, EMWeightMapper<Arc> (emmap_));
      
      alpha.clear ();
      beta.clear ();
      ShortestDistance (fst, &alpha);
      ShortestDistance (fst, &beta, true);
      
      for (StateIterator<VectorFst<Arc> > siter (fst); !siter.Done ();
	   siter.Next ()) {
	StateId q = siter.Value ();
	for (ArcIterator<VectorFst<Arc> > aiter (fst, q); !aiter.Done ();
	     aiter.Next ()) {
	  const Arc& arc = aiter.Value ();
	  Weight gamma = Divide (Times (Times (alpha [q], arc.weight),
					beta [arc.nextstate]), beta [0]);
	  if (gamma.Value () == gamma.Value ()) {
	    emmap_prev_[arc.ilabel][arc.olabel] =
	      Plus (emmap_prev_[arc.ilabel][arc.olabel], gamma);
	    total_ = Plus (total_, gamma);
	  }
	}
      }
    }
  }

  void PrintOneBestAlignments () {
    tmpfar_reader_->Reset ();
    for (int i = 1; !tmpfar_reader_->Done (); tmpfar_reader_->Next (), ++i) {
      VectorFst<StdArc> ifst;
      Map (tmpfar_reader_->GetFst (), 
	   &ifst,
	   EMWeightPathMapper<Arc,StdArc> (emmap_, ipenalties_, openalties_));
      VectorFst<StdArc> ofst;
      ShortestPath (ifst, &ofst);
      TopSort (&ofst);
      for (StateIterator<VectorFst<StdArc> > siter (ofst); !siter.Done ();
	   siter.Next ()) {
	StateId q = siter.Value ();
	for (ArcIterator<VectorFst<StdArc> > aiter (ofst, q); !aiter.Done ();
	     aiter.Next ()) {
	  const StdArc& arc = aiter.Value ();
	  cout << isyms_.Find (arc.ilabel) << mdelim_ << osyms_.Find (arc.olabel)
	       << ((ofst.Final (arc.nextstate) != Weight::Zero ()) ? "\n" : " ");
	}
      }      
    }
  }

  //Generate and trim an individual alignment FST
  void EntryToFst (const Entry& entry, F& fst) {
    const vector<uint32_t>& graphemes = entry.first;
    const vector<uint32_t>& phonemes  = entry.second;
    StateId istate, ostate, final;
    Label ilabel, olabel;

    final = ((graphemes.size () + 1) * (phonemes.size () + 1) - 1);
    fst.ReserveStates (final);
    for (int i = 0; i <= final; i++)
      fst.AddState ();
    fst.SetStart (0);

    for (int i = 0; i <= graphemes.size (); i++) {
      for (int j = 0; j <= phonemes.size (); j++) {
	istate = i * (phonemes.size () + 1) + j;

	//Insertions
	if (delg_ == true) {
	  for (int m = 1; m <= pmax_; m++) {
	    if (j + m <= phonemes.size ()) {
	      ilabel = 0;
	      olabel = AddSubVecSymbol (phonemes, j, j + m, false);
	      ostate = i * (phonemes.size () + 1) + j + m;
	      fst.AddArc (istate, Arc (ilabel, olabel, Weight::One (), ostate));
	      AddMultigram (ilabel, olabel);
	    }
	  }
	}
	
	//Deletions
	if (delp_ == true) {
	  for (int n = 1; n <= gmax_; n++) {
	    if (i + n <= graphemes.size ()) {
	      ilabel = AddSubVecSymbol (graphemes, i, i + n);
	      olabel = 0;
	      ostate = (i + n) * (phonemes.size () + 1) + j;
	      fst.AddArc (istate, Arc (ilabel, olabel, Weight::One (), ostate));
	      AddMultigram (ilabel, olabel);
	    }
	  }
	}
	
	//Normal arcs
	for (int n = 1; n <= gmax_; n++) {
	  for (int m = 1; m <= pmax_; m++) {
	    if (i + n <= graphemes.size () && j + m <= phonemes.size ()) {
	      if (! (m > 1 && n > 1)) {
		ilabel = AddSubVecSymbol (graphemes, i, i + n);
		olabel = AddSubVecSymbol (phonemes, j, j + m, false);
		ostate = (i + n) * (phonemes.size () + 1) + (j + m);
		fst.AddArc (istate, Arc (ilabel, olabel, 
					 Weight::One (), ostate));
		AddMultigram (ilabel, olabel, Weight::One().Value () * (m + n) );
	      }
	    }
	  }
	}
      }
    }
    
    //Final state
    fst.SetFinal (final, Weight::One ());
    Connect (&fst);
    return;
  }

  //Convenience function for dumping an OpenFst compatible
  // text-format FST with all necessary arcs.  
  //This version will generate unconnected arcs in the case
  // where only one of delg_ or delp_ are set to 'true'.
  void PrintEntryToFst (const Entry& entry) {
    const vector<uint32_t>& graphemes = entry.first;
    const vector<uint32_t>& phonemes  = entry.second;
    int istate, ostate;

    for (int i = 0; i <= graphemes.size (); i++) {
      for (int j = 0; j <= phonemes.size (); j++) {
	istate = i * (phonemes.size () + 1) + j;

	//Insertions
	if (delg_ == true) {
	  for (int m = 1; m <= pmax_; m++) {
	    if (j + m <= phonemes.size ()) {
	      ostate = i * (phonemes.size () + 1) + j + m;
	      cout << istate << " " << ostate << " _ ";
	      PrintSubVec (phonemes, j, j + m);
	      cout << "\n";
	    }
	  }
	}
	//Deletions
	if (delp_ == true) {
	  for (int n = 1; n <= gmax_; n++) {
	    if (i + n <= graphemes.size ()) {
	      ostate = (i + n) * (phonemes.size () + 1) + j;
	      cout << istate << " " << ostate << " ";
	      PrintSubVec (graphemes, i, i + n);
	      cout << " _" << "\n";
	    }
	  }
	}
	//Normal arcs
	for (int n = 1; n <= gmax_; n++) {
	  for (int m = 1; m <= pmax_; m++) {
	    if (i + n <= graphemes.size () && j + m <= phonemes.size ()) {
	      if (! (m > 1 && n > 1)) {
		ostate = (i + n) * (phonemes.size () + 1) + (j + m);
		cout << istate << " " << ostate << " ";
		PrintSubVec (graphemes, i, i + n);
		cout << " ";
		PrintSubVec (phonemes, j, j + m);
		cout << "\n";
	      }
	    }
	  }
	}
      }
    }

    //Final state
    cout << ((graphemes.size () + 1) * (phonemes.size () + 1) - 1) << "\n";
    return;
  }

  //Convenience function for debugging: dump the existing weight map;
  void PrintEMMap () {
    for (typename EMMap::iterator it = emmap_.begin (); it != emmap_.end (); ++it) {
      cout << isyms_.Find (it->first) << endl;
      for (typename LWMap::iterator jt = emmap_[it->first].begin (); 
	   jt != emmap_[it->first].end (); ++jt) 
	cout << "\t" << osyms_.Find (jt->first) << " " << jt->second << endl;
    }
  }

  ~Aligner () {
    delete tmpfar_reader_;
  }

 private:
  void PrintSubVec (const vector<uint32_t>& subseq, int begin, int end) {
    for (int i = begin; i < end; i++) {
      string utfchars;
      for (int j = 0; j < utfmap_[subseq [i]].size (); j++) 
	utf8::append (utfmap_[subseq [i]][j], back_inserter (utfchars));
      
      cout << utfchars << ((i < end - 1) ? "," : "");
    }
  }

  int AddSubVecSymbol (const vector<uint32_t>& subseq, int begin, int end, 
		   bool input = true) {
    string utfchars;
    for (int i = begin; i < end; i++) {
      for (int j = 0; j < utfmap_[subseq [i]].size (); j++) 
	utf8::append (utfmap_[subseq [i]][j], back_inserter (utfchars));
      if (i < end - 1)
	utfchars.append (cdelim_);
    }
    
    if (input) {
      int i = isyms_.AddSymbol (utfchars);
      ipenalties_[i] = end - begin;
      return i;
    } else {
      int i = osyms_.AddSymbol (utfchars);
      openalties_[i] = end - begin;
      return i;
    }
  }

  void AddMultigram (int ilabel, int olabel, Weight weight = Weight::One()) {
    if (olabel == 0 || ilabel == 0)
      weight = 99;

    if (emmap_prev_.find (ilabel) == emmap_prev_.end ()) {
      LWMap lw;
      lw [olabel] = weight;
      emmap_prev_[ilabel] = lw;
    } else if (emmap_prev_[ilabel].find (olabel) == emmap_prev_[ilabel].end ()) {
      emmap_prev_[ilabel][olabel] = Plus (emmap_prev_[ilabel][olabel], weight);
    }
    total_ = Plus (total_, weight);
    //Else do nothing - no need to tally anything during the initialization phase.
  }

  Weight CountEMMap () {
    Weight w = Weight::Zero ();
    for (typename EMMap::iterator it = emmap_.begin (); 
	 it != emmap_.end (); ++it)
      w = Plus (w, -log (emmap_[it->first].size ()));
    return w;
  }
  
  double GetChange () {
    return abs (total_prev_.Value () - total_.Value ());
  }

  //Mapper to update an input alignment lattice with the current 
  // EM arc weights.
  template <class A, class B = A>
  struct EMWeightMapper {
    typedef A FromArc;
    typedef B ToArc;
    typedef typename FromArc::Weight FromWeight;
    typedef typename ToArc::Weight ToWeight;

    EMWeightMapper (EMMap& emmap) : emmap_(emmap) {}

    B operator () (const A& arc) const {
      //Careful with the final weights
      if (arc.nextstate == kNoStateId)
	return ToArc (arc.ilabel, arc.olabel, arc.weight, arc.nextstate);

      ToWeight w = emmap_[arc.ilabel][arc.olabel];
      if (w == ToWeight::Zero () || w != w)
	w = 99;
      return ToArc (arc.ilabel, arc.olabel, w, arc.nextstate);
    }

    MapFinalAction FinalAction () const { return MAP_NO_SUPERFINAL; }

    MapSymbolsAction InputSymbolsAction () const { return MAP_COPY_SYMBOLS; }

    MapSymbolsAction OutputSymbolsAction () const { return MAP_COPY_SYMBOLS; }

    uint64 Properties(uint64 props) const {
      return props & kWeightInvariantProperties;
    }
    
    EMMap& emmap_;      
  };

  //Mapper to update an input alignment lattice with the current 
  // EM arc weights. This also converts to StdArc for use in shortest path.
  template <class A, class B = A>
  struct EMWeightPathMapper {
    typedef A FromArc;
    typedef B ToArc;
    typedef typename FromArc::Weight FromWeight;
    typedef typename ToArc::Weight ToWeight;

    EMWeightPathMapper (EMMap& emmap, PenaltyMap& imap, PenaltyMap& omap)
      : emmap_(emmap), imap_(imap), omap_(omap) {}

    B operator () (const A& arc) const {
      if (arc.nextstate == kNoStateId)
	return B (arc.ilabel, arc.olabel, 
		  convert_weight_(arc.weight), arc.nextstate);

      ToWeight w = convert_weight_(emmap_[arc.ilabel][arc.olabel]);
      if (w == ToWeight::Zero () || w != w)
	w = 99;
      
      w = w.Value () * max (imap_[arc.ilabel], omap_[arc.olabel]);

      return B (arc.ilabel, arc.olabel, w, arc.nextstate);
    }

    MapFinalAction FinalAction () const { return MAP_NO_SUPERFINAL; }

    MapSymbolsAction InputSymbolsAction () const { return MAP_COPY_SYMBOLS; }

    MapSymbolsAction OutputSymbolsAction () const { return MAP_COPY_SYMBOLS; }

    uint64 Properties(uint64 props) const {
      return props & kWeightInvariantProperties;
    }
    
   private:
    EMMap& emmap_;
    PenaltyMap& imap_;
    PenaltyMap& omap_;
    WeightConvert<FromWeight, ToWeight> convert_weight_;
  };

  inline bool file_exists (const string& filename) {
    struct stat buffer;
    return (stat (filename.c_str (), &buffer) == 0);
  }

  const DictionaryReader& reader_;
  UTFMap& utfmap_;
  int gmax_;
  int pmax_;
  bool delg_;
  bool delp_;
  const string mdelim_;
  const string cdelim_;
  SymbolTable& isyms_;
  SymbolTable& osyms_;
  const string& far_name_;
  MaxType max_type_;
  int verbose_;

  EMMap emmap_;
  EMMap emmap_prev_;
  FarReader<Arc>* tmpfar_reader_;
  Weight total_;
  Weight total_prev_;
  PenaltyMap ipenalties_;
  PenaltyMap openalties_;
};
#endif //ALIGNER_H__
