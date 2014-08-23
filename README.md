AltFstAligner
=============

An alternative OpenFst-based aligner for grapheme-to-phoneme conversion.

USAGE:
```
#Basic usage
$ cd src/ && make install && cd ..
$ ./altfst-align --corpus=input.dict --verbose=2 > output.corpus

#For numerous options
$ ./altfst-align --help
```
WARNING:
 * Haven't tested on linux yet.
 * Should be UTF8 compatible, but no extensive tests yet.



More or less the same thing that is in Phonetisaurus.
Intended to replace M2MAligner, this version stores the alignment
Fsts offline in a FarArchive.

Turns out this is both way lower memory (25mb vs 300mb...) but 
it is also faster.

Results are a tiny bit different to the original and produce different
accuracies in the G2P model-building and decoding stages:

Original:
```
Words: 12000  Hyps: 12000 Refs: 12000
######################################################################
                          EVALUATION RESULTS
----------------------------------------------------------------------
(T)otal tokens in reference: 75790
(M)atches: 71721  (S)ubstitutions: 3630  (I)nsertions: 362  (D)eletions: 439
% Correct (M/T)           -- %94.63
% Token ER ((S+I+D)/T)    -- %5.85
% Accuracy 1.0-ER         -- %94.15
       --------------------------------------------------------
(S)equences: 12000  (C)orrect sequences: 9070  (E)rror sequences: 2930
% Sequence ER (E/S)       -- %24.42
% Sequence Acc (1.0-E/S)  -- %75.58
######################################################################
```
AltFstAligner:
```
Words: 12000  Hyps: 12000 Refs: 12000
######################################################################
                          EVALUATION RESULTS
----------------------------------------------------------------------
(T)otal tokens in reference: 75774
(M)atches: 71692  (S)ubstitutions: 3630  (I)nsertions: 366  (D)eletions: 452
% Correct (M/T)           -- %94.61
% Token ER ((S+I+D)/T)    -- %5.87
% Accuracy 1.0-ER         -- %94.13
       --------------------------------------------------------
(S)equences: 12000  (C)orrect sequences: 9066  (E)rror sequences: 2934
% Sequence ER (E/S)       -- %24.45
% Sequence Acc (1.0-E/S)  -- %75.55
######################################################################
```
I'm not sure exactly where the difference is coming from yet, but I suspect
that it is not statistically significant. Could be nice for model combination.

