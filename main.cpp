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
           vector<string> &solutions) {
  if (seq.size() == 0 || seq.size() > 255)
    // FIXME
    return;
  gcsa::range_type range = index.find(seq);
  // if (gcsa::Range::empty(range))
  //   return;
  vector<gcsa::node_type> results;
  index.locate(range, results);
  cout << seq.size() << " " << results.size() << endl;
  for (const auto &r : results)
    solutions.push_back(gcsa::Node::decode(r));
}

void dump(const tuple<string, string, string> &input,
          const vector<string> &output) {
  for (const auto &o : output)
    cout << get<0>(input) << " " << o << endl;
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
  vector<vector<string>> output(opt::b);

  gzFile fp = gzopen(opt::fx_path.c_str(), "r");
  kseq_t *seq = kseq_init(fp);
  int i = 0;
  int n = 1;
  bool is_last_batch = false;
  while (!is_last_batch) {
    cerr << "Loading   batch " << n << ".."
         << "\r";
    is_last_batch = load_batch(seq, input);
    cerr << "Analyzing batch " << n << ".."
         << "\r";
#pragma omp parallel for num_threads(opt::t)
    for (i = 0; i < input.size(); ++i)
      query(index, get<1>(input[i]), opt::f, output[i]);
    cerr << "Dumping   batch " << n << ".."
         << "\r";
    for (i = 0; i < input.size(); ++i)
      dump(input[i], output[i]);
    ++n;
  }
  cerr << "Dumped " << n - 1 << " batches. Done." << endl;

  kseq_destroy(seq);
  gzclose(fp);

  return 0;
}
