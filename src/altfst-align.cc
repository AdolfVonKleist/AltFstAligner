#include <fst/fstlib.h>
#include <fst/extensions/far/far.h>
#include <fst/extensions/far/stlist.h>
#include <fst/extensions/far/sttable.h>
#include <fst/fst.h>
#include <fst/vector-fst.h>
#include "utf8.h"
#include "DictionaryReader.h"
#include "Aligner.h"
using namespace fst;

DEFINE_string (corpus,   "", "The input pronunciation dictionary to align.");
DEFINE_string (prefix,   "_tmpalign", "The file prefix for temporary .far storage.");
DEFINE_string (gdelim,   "", "Grapheme delimiter.");
DEFINE_string (pdelim,   " ", "Phoneme delimiter.");
DEFINE_string (edelim,   "\t", "Word/Pron delimiter.");
DEFINE_string (cdelim,   "|", "Output chunk delimiter: P|H}f.");
DEFINE_string (mdelim,   "}", "Output graph/phon multigram delimiter: P|H}f.");
DEFINE_string (eps,      "_", "Epsilon symbol used to represent insertions/deletions.");
DEFINE_string (max_type, "joint", "Maximization type: join, pg, [gp].");
DEFINE_double (thresh,   -1, "Minimum change for EM termination.");
DEFINE_bool   (ins,      false, "Permit insertions (graphemic nulls).");
DEFINE_bool   (del,      true, "Permit deletions (phonemic nulls).");
DEFINE_int32  (max_iter, 12, "Maximum number of EM iterations.");
DEFINE_int32  (max_g,    2, "Maximum subsequence length for graphemes.");
DEFINE_int32  (max_p,    2, "Maximum subsequence length for phonemes.");
DEFINE_int32  (verbose,  0, "Verbosity level."); 

int main (int argc, char* argv []) {
  string usage = "altfst-align --corpus=input.dict > aligned.corpus\n\n Usage:";
  set_new_handler (FailedNewHandler);
  SetFlags (usage.c_str (), &argc, &argv, false);
  if (FLAGS_corpus.compare ("") == 0) {
    cerr << argv [0] << ": --corpus parameter must be set!" << endl;
    cerr << " See --help for details." << endl;
    exit (1);
  }

  MaxType max_type = JOINT_MAX;
  if (FLAGS_max_type.compare ("joint") == 0)
    max_type = JOINT_MAX;
  else if (FLAGS_max_type.compare ("pg") == 0)
    max_type = PG_MAX;
  else if (FLAGS_max_type.compare ("gp") == 0) {
    cerr << "--max_type 'gp' not yet supported." << endl;
    cerr << " See --help for details." << endl;
    exit (1);
  } else {
    cerr << "--max_type must be one of 'join', 'gp', 'pg'." << endl;
    cerr << " See --help for details." << endl;
  }

  UTFMap utfmap;
  DictionaryReader reader (FLAGS_gdelim, FLAGS_pdelim, FLAGS_edelim, utfmap);
  reader.ReadTrainingFile (FLAGS_corpus);

  SymbolTable isyms;
  SymbolTable osyms;
  isyms.AddSymbol (FLAGS_eps);
  osyms.AddSymbol (FLAGS_eps);

  Aligner<VectorFst<LogArc> > aligner (reader, utfmap, 
				       FLAGS_max_g, FLAGS_max_p, 
				       FLAGS_ins, FLAGS_del, 
				       FLAGS_mdelim, FLAGS_cdelim,
				       isyms, osyms,
				       FLAGS_prefix, max_type,
				       FLAGS_verbose);
  aligner.InitModel ();
  for (int i = 0; i < FLAGS_max_iter + 1; i++) {
    aligner.ExpectationStep ();
    double change = aligner.MaximizationStep ();

    if (FLAGS_verbose >= 1 && i > 0)
      cerr << "Iter: " << i << " Change: " << change << endl;

    if (FLAGS_thresh > 0  && i > 0 && change < FLAGS_thresh)
      break;
  }

  aligner.PrintOneBestAlignments ();

  return 0;
}
