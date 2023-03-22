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

void query(const gcsa::GCSA &index, const string &seq, int flank,
           vector<pair<int, int>> &solutions) {
  gcsa::range_type range;
  vector<gcsa::node_type> results;
  vector<pair<int, int>> specifics;
  int seql = seq.size();
  int beg = seql - 1;
  int end = beg;
  int p;
  while (beg >= 0) {
    // cerr << "Starting from " << beg << "," << end << endl;
    // Backward search until first mismatch
    range = index.charRange(index.alpha.char2comp[seq[beg]]);
    --beg;
    while (beg >= 0 && gcsa::Range::length(range) > 0) {
      // cerr << "Backward extending to " << beg << endl;
      range = index.LF(range, index.alpha.char2comp[seq[beg]]);
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
    range = index.charRange(index.alpha.char2comp[RCN[seq[p]]]);
    // range = rcindex.charRange(rcindex.alpha.char2comp[seq->seq.s[p]]);
    ++p;
    while (p <= end && gcsa::Range::length(range) > 0) {
      // cerr << "Forward extending to " << p << endl;
      range = index.LF(range, index.alpha.char2comp[RCN[seq[p]]]);
      // range = rcindex.LF(range, rcindex.alpha.char2comp[seq->seq.s[p]]);
      ++p;
    }
    --p;
    // cerr << "Mismatch at " << p << endl;
    specifics.push_back(
        make_pair(max(0, beg - flank), min(p + flank, seql - 1)));
    beg = p - 1;
    end = beg;
    if (specifics.back().first == 0)
      break;
  }
  vector<pair<int, int>> assembled_specifics;
  assemble(specifics, assembled_specifics);
  solutions = assembled_specifics;
}

void dump(const tuple<string, string, string> &input,
          const vector<pair<int, int>> &output) {
  for (const auto &o : output) {
    string S(get<1>(input), o.first, o.second - o.first + 1);
    if (get<2>(input).size() == 0) {
      cout << ">" << get<0>(input) << ":" << o.first << "-" << o.second << "\n"
           << S << endl;
    } else {
      string Q(get<2>(input), o.first, o.second - o.first + 1);
      cout << "@" << get<0>(input) << ":" << o.first << "-" << o.second << "\n"
           << S << "\n"
           << "+"
           << "\n"
           << Q << endl;
    }
  }
}

bool load_batch(kseq_t *seq, vector<tuple<string, string, string>> &input) {
  int l = 0;
  int i = 0;
  while (i < opt::b && (l = kseq_read(seq)) >= 0) {
    if (seq->qual.l == 0)
      input[i] = make_tuple(string(seq->name.s), string(seq->seq.s), "");
    else
      input[i] = make_tuple(string(seq->name.s), string(seq->seq.s),
                            string(seq->qual.s));
    ++i;
  }
  return l < 0;
}

int main(int argc, char *argv[]) {
  parse_arguments(argc, argv);

  cerr << "Loading the index..";
  ifstream ifs;
  ifs.open(opt::index_path, std::ifstream::in);
  gcsa::GCSA index;
  index.load(ifs);
  ifs.close();
  cerr << " Done." << endl;

  vector<tuple<string, string, string>> input(opt::b);
  vector<vector<pair<int, int>>> output(opt::b);

  gzFile fp = gzopen(opt::fx_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int i = 0;
  int n = 0;
  bool is_last_batch = false;
  while (!is_last_batch) {
    cerr << "Loading   batch " << n << ".."
         << "\r";
    is_last_batch = load_batch(seq, input);
    cerr << "Analyzing batch " << n << ".."
         << "\r";
#pragma omp parallel for num_threads(opt::t)
    for (i = 0; i < input.size(); ++i) {
      query(index, get<1>(input[i]), opt::f, output[i]);
    }
    cerr << "Dumping   batch " << n << ".."
         << "\r";
    for (i = 0; i < input.size(); ++i) {
      dump(input[i], output[i]);
    }
    ++n;
  }
  cerr << "Dumped " << n - 1 << " batches. Done." << endl;

  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}
