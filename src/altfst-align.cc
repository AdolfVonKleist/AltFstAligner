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
DEFINE_double (thresh,   0.001, "Minimum change for EM termination.");
DEFINE_bool   (ins,      false, "Permit insertions (graphemic nulls).");
DEFINE_bool   (del,      true, "Permit deletions (phonemic nulls).");
DEFINE_int32 (max_iter, 12, "Maximum number of EM iterations.");
DEFINE_int32  (max_g,    2, "Maximum subsequence length for graphemes.");
DEFINE_int32  (max_p,    2, "Maximum subsequence length for phonemes.");
DEFINE_int32  (verbose,  0, "Verbosity level."); 

int main (int argc, char* argv []) {
  string usage = "altfst-align --corpus=input.dict > aligned.corpus\n\n Usage:";
  set_new_handler (FailedNewHandler);
  SetFlags (usage.c_str (), &argc, &argv, false);

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
				       isyms, osyms,
				       FLAGS_prefix, JOINT_MAX,
				       FLAGS_verbose);
  aligner.InitModel ();
  for (int i = 0; i < FLAGS_max_iter + 1; i++) {
    aligner.ExpectationStep ();
    double change = aligner.MaximizationStep ();

    if (FLAGS_verbose >= 1 && i > 0)
      cerr << "Iter: " << i << " Change: " << change << endl;

    if (i > 0 && change < FLAGS_thresh)
      break;
  }

  aligner.PrintOneBestAlignments ();

  return 1;
}
