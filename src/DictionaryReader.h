#ifndef DICTIONARY_READER_H__
#define DICTIONARY_READER_H__
#include <fst/fstlib.h>
#include "utf8.h"
using namespace fst;
typedef vector<vector<uint32_t> > UTF8TokenSeq;
typedef pair<vector<uint32_t>, vector<uint32_t> > Entry;
typedef unordered_map<uint32_t, vector<uint32_t> > UTFMap;

struct TokenHasher {
public:
  //This is probably overkill
  uint32_t operator () (const vector<uint32_t>& t) const {
    uint32_t h = 0;
    for (int i = 0; i < t.size (); i++) {
      uint32_t high = h & 0xf8000000;
      h = h << 5;
      h = h ^ (high >> 27);
      h = h ^ t [i];
    }
    return h;
  }
};
typedef unordered_map<vector<uint32_t>, uint32_t, TokenHasher> TokenMap;


class DictionaryReader {
 public:
  DictionaryReader (const string& gdelim, const string& pdelim, 
		    const string& edelim, UTFMap& utfmap) 
    : gdelim_(CharToUnicode (gdelim)), 
    pdelim_(CharToUnicode (pdelim)),
    edelim_(CharToUnicode (edelim)),
    utfmap_(utfmap) { }

  void UTF8StringToUnicode (const string& utf8string, 
			    UTF8TokenSeq& unicode,
			    uint32_t delim) {
    char* str = (char*)utf8string.c_str ();
    char* itr = str;
    char* end = str + strlen (str) + 1;
    
    vector<uint32_t> token;
    do {
      uint32_t code = utf8::next (itr, end);
      if (code == 0)
	continue;

      if (code == delim || delim == 0) {
	unicode.push_back (token);
	token.clear ();
      } else {
	token.push_back (code);
      }
      
    } while (itr < end);

    unicode.push_back (token);
  }

  uint32_t Hasher (const vector<uint32_t>& token) {
    if (token_map_.find (token) == token_map_.end ()) {
      uint32_t h = token_map_.size () + 1;
      token_map_[token] = h;
      return h;
    } else {
      return token_map_[token];
    }
  }

  void TokenizeInput (const string& input) {
    UTF8TokenSeq tokens;
    Entry entry;
    uint32_t h;
    vector<uint32_t> token;
    UTF8StringToUnicode (input, tokens, edelim_);

    if (gdelim_ == 0) {
      for (int i = 0; i < tokens [0].size (); i++) {
	h = Hasher ({tokens [0][i]});
	utfmap_[h] = {tokens [0][i]};
	entry.first.push_back (h);
      }
    } else {
      for (int i = 0; i < tokens [0].size (); i++) {
	if (tokens [0][i] == gdelim_) {
	  h = Hasher (token);
	  utfmap_[h] = token;
	  token.clear ();
	  entry.first.push_back (h);
	} else {
	  token.push_back (tokens [0][i]);
	}
      }
      h = Hasher (token);
      utfmap_[h] = token;
      entry.first.push_back (h);
    }

    if (pdelim_ == 0) {
      for (int i = 0; i < tokens [1].size (); i++) {
	h = Hasher ({tokens [1][i]});
	utfmap_[h] = {tokens [1][i]};
	entry.second.push_back (h);
      }
    } else {
      for (int i = 0; i < tokens [1].size (); i++) {
	if (tokens [1][i] == pdelim_) {
	  h = Hasher (token);
	  utfmap_[h] = token;
	  token.clear ();
	  entry.second.push_back (h);
	} else {
	  token.push_back (tokens [1][i]);
	}
      }
      h = Hasher (token);
      utfmap_[h] = token;
      entry.second.push_back (h);
    }

    entries.push_back (entry);
  }

  void ReadTrainingFile (const string& infile) {
    ifstream ifp (infile.c_str ());
    string line;
      
    if (ifp.is_open ()) {
      while (ifp.good ()) {
	getline (ifp, line);
	if (line.empty ())
	  continue;
	TokenizeInput (line);
      }
    }    
  }

  uint32_t CharToUnicode (const string& delim) {
    char* itr = (char*)delim.c_str ();
    uint32_t delim_code 
      = utf8::next (itr, (char*)delim.c_str () + strlen (delim.c_str ()) + 1);
    return delim_code;
  }

  vector<Entry> entries;

 private:
  uint32_t gdelim_;
  uint32_t pdelim_;
  uint32_t edelim_;
  UTFMap& utfmap_;
  TokenMap token_map_;
};
#endif //DICTIONARY_READER_H__
