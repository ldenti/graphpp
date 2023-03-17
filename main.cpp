#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <string>
#include <vector>
#include <zlib.h>

#include "gcsa/gcsa.h"
#include "gcsa/lcp.h"
#include "kseq.h"

#include "argument_parser.hpp"

KSEQ_INIT(gzFile, gzread)

static const char RCN[128] = {
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   //  0
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 10
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 20
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 30
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 40
    0,   0,   0, 0,   0,   0,   0,   0,   0,   0,   // 50
    0,   0,   0, 0,   0,   'T', 0,   'G', 0,   0,   // 60
    0,   'C', 0, 0,   0,   0,   0,   0,   'N', 0,   // 70
    0,   0,   0, 0,   'A', 0,   0,   0,   0,   0,   // 80
    0,   0,   0, 0,   0,   0,   0,   'T', 0,   'G', // 90
    0,   0,   0, 'G', 0,   0,   0,   0,   0,   0,   // 100
    'N', 0,   0, 0,   0,   0,   'A', 0,   0,   0,   // 110
    0,   0,   0, 0,   0,   0,   0,   0              // 120
};

using namespace std;

int main(int argc, char *argv[]) {
  parse_arguments(argc, argv);

  // Loading the index
  cerr << "Loading the index.." << endl;
  ifstream ifs;
  ifs.open(opt::index_path, std::ifstream::in);
  gcsa::GCSA index;
  index.load(ifs);
  ifs.close();
  ifs.open(opt::index_path + ".lcp", std::ifstream::in);
  gcsa::LCPArray lcp;
  lcp.load(ifs);
  ifs.close();

  gzFile fp = gzopen(opt::fx_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  gcsa::range_type range;
  vector<gcsa::node_type> results;
  int beg, end, p;
  vector<pair<int, int>> specifics;
  while ((l = kseq_read(seq)) >= 0) {
    cerr << "Querying " << seq->name.s << ".." << endl;
    beg = l - 1;
    end = beg;
    while (beg >= 0) {
      // cerr << "Starting from " << beg << "," << end << endl;
      // Backward search until first mismatch
      range = index.charRange(index.alpha.char2comp[seq->seq.s[beg]]);
      --beg;
      while (beg >= 0 && gcsa::Range::length(range) > 0) {
        // cerr << "Backward extending to " << beg << endl;
        range = index.LF(range, index.alpha.char2comp[seq->seq.s[beg]]);
        --beg;
      }
      if (beg < 0)
        // no specific since we matched the entire query
        continue;
      // we must have a specific
      ++beg;
      // cerr << "Mismatch at " << beg << endl;

      // Forward search until first mismatch
      p = beg;
      // cerr << "Starting from " << p << "," << end << endl;
      range = index.charRange(index.alpha.char2comp[RCN[seq->seq.s[p]]]);
      // range = rcindex.charRange(rcindex.alpha.char2comp[seq->seq.s[p]]);
      ++p;
      while (p <= end && gcsa::Range::length(range) > 0) {
        // cerr << "Forward extending to " << p << endl;
        range = index.LF(range, index.alpha.char2comp[RCN[seq->seq.s[p]]]);
        // range = rcindex.LF(range, rcindex.alpha.char2comp[seq->seq.s[p]]);
        ++p;
      }
      --p;
      // cerr << "Mismatch at " << p << endl;
      // cout << "N " << seq->name.s << " " << beg << " " << p << " "
      //      << p - beg + 1 << endl;
      specifics.push_back(make_pair(beg, p));
      beg = p - 1;
      end = beg;
    }

    // assemble specifics
    cerr << "Assembling (round 1).." << endl;
    vector<pair<int, int>> assembled_specifics;
    int i = specifics.size() - 1;
    while (i >= 0) {
      int j;
      for (j = i - 1; j >= 0; --j) {
        if (specifics[j + 1].second < specifics[j].first) {
          // non-overlapping
          assembled_specifics.push_back(make_pair(
              max(0, specifics[i].first - opt::f),
              min(specifics[j + 1].second + opt::f, (int)(seq->seq.l - 1))));
          i = j;
          break;
        }
      }
      if (j < 0) {
        assembled_specifics.push_back(
            make_pair(max(0, specifics[i].first - opt::f),
                      min(specifics[0].second + opt::f, (int)(seq->seq.l - 1))));
        i = j;
      }
    }

    reverse(assembled_specifics.begin(), assembled_specifics.end());
    
    // assemble assembled specifics (2nd round due to flanking)
    cerr << "Assembling (round 2).." << endl;
    vector<pair<uint, uint>> final_assembled_specifics;
    i = assembled_specifics.size() - 1;
    while (i >= 0) {
      int j;
      for (j = i - 1; j >= 0; --j) {
        if (assembled_specifics[j + 1].second < assembled_specifics[j].first) {
          // non-overlapping
          final_assembled_specifics.push_back(make_pair(
              assembled_specifics[i].first, assembled_specifics[j + 1].second));
          i = j;
          break;
        }
      }
      if (j < 0) {
        final_assembled_specifics.push_back(make_pair(
            assembled_specifics[i].first, assembled_specifics[0].second));
        i = j;
      }
    }

    cerr << "Dumping.." << endl;
    for (const auto &spec : final_assembled_specifics) {
      string Q(seq->seq.s, spec.first, spec.second - spec.first + 1);
      cout << ">" << seq->name.s << ":" << spec.first << "-" << spec.second
           << "\n"
           << Q << endl;
    }
  }

  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}
