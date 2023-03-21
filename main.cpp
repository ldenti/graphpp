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

void assemble(const vector<pair<int, int>> &specifics,
              vector<pair<int, int>> &assembled_specifics) {
  int i = specifics.size() - 1;
  while (i >= 0) {
    int j;
    for (j = i - 1; j >= 0; --j) {
      if (specifics[j + 1].second < specifics[j].first) {
        // non-overlapping
        assembled_specifics.push_back(
            make_pair(specifics[i].first, specifics[j + 1].second));
        i = j;
        break;
      }
    }
    if (j < 0) {
      assembled_specifics.push_back(
          make_pair(specifics[i].first, specifics[0].second));
      i = j;
    }
  }
}

void query(gcsa::GCSA index, kseq_t *seq, int flank,
           vector<pair<int, int>> &solutions) {
  cerr << "Querying " << seq->name.s << ".." << endl;
  gcsa::range_type range;
  vector<gcsa::node_type> results;
  vector<pair<int, int>> specifics;
  int beg = seq->seq.l - 1;
  int end = beg;
  int p;
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
    specifics.push_back(
        make_pair(max(0, beg - flank), min(p + flank, (int)(seq->seq.l - 1))));
    beg = p - 1;
    end = beg;
    if (specifics.back().first == 0)
      break;
  }
  cerr << "Assembling " << specifics.size() << " specific strings.." << endl;
  vector<pair<int, int>> assembled_specifics;
  assemble(specifics, assembled_specifics);
  solutions = assembled_specifics;
}

void dump(kseq_t *seq, const vector<pair<int, int>> &solutions) {
  cerr << "Dumping " << solutions.size() << " strings.." << endl;
  for (const auto &spec : solutions) {
    string S(seq->seq.s, spec.first, spec.second - spec.first + 1);
    if (seq->qual.l == 0) {
      cout << ">" << seq->name.s << ":" << spec.first << "-" << spec.second
           << "\n"
           << S << endl;
    } else {
      string Q(seq->qual.s, spec.first, spec.second - spec.first + 1);
      cout << "@" << seq->name.s << ":" << spec.first << "-" << spec.second
           << "\n"
           << S << "\n"
           << "+"
           << "\n"
           << Q << endl;
    }
  }
}

int main(int argc, char *argv[]) {
  parse_arguments(argc, argv);

  // Loading the index
  cerr << "Loading the index.." << endl;
  ifstream ifs;
  ifs.open(opt::index_path, std::ifstream::in);
  gcsa::GCSA index;
  index.load(ifs);
  ifs.close();

  gzFile fp = gzopen(opt::fx_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int l;
  while ((l = kseq_read(seq)) >= 0) {
    vector<pair<int, int>> solutions;
    query(index, seq, opt::f, solutions);
    dump(seq, solutions);
  }

  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}
