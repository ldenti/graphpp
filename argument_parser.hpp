#ifndef _GRAPHPP_PARSER_HPP_
#define _GRAPHPP_PARSER_HPP_

#include <sstream>

#include <getopt.h>

using namespace std;

static const char *USAGE_MESSAGE =
    "Usage: graphpp [-f FLANK] <graph.gcsa2> <query.fx>\n"
    "\n"
    "      -f, --flank            size of flanking region on both side of a "
    "specific string (default:0)\n"
    "      -h, --help             display this help and exit\n"
    "      -t, --threads          number of threads (default: 1)\n"
    "      -b, --bsize            batch size (default: 100000)\n"
    // "      -v, --verbose                     verbose mode (default: false)\n"
    "\n";

namespace opt {
static int f = 0;
static int b = 100000;
static int t = 1;
static string index_path;
static string fx_path;
} // namespace opt

static const char *shortopts = "f:b:t:h";

static const struct option longopts[] = {
    {"flank", required_argument, NULL, 'f'},
    {"bsize", required_argument, NULL, 'b'},
    {"threads", required_argument, NULL, 't'},
    {"help", no_argument, NULL, 'h'},
    {NULL, 0, NULL, 0}};

void parse_arguments(int argc, char **argv) {
  for (char c;
       (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
    istringstream arg(optarg != NULL ? optarg : "");
    switch (c) {
    case 'f':
      arg >> opt::f;
      break;
    case 'b':
      arg >> opt::b;
      break;
    case 't':
      arg >> opt::t;
      break;
    // case 'v':
    //   opt::verbose = true;
    //   break;
    case '?':
      cerr << USAGE_MESSAGE;
      exit(EXIT_FAILURE);
      break;
    case 'h':
      cout << USAGE_MESSAGE;
      exit(EXIT_SUCCESS);
    }
  }

  if (argc - optind < 2) {
    cerr << USAGE_MESSAGE;
    exit(EXIT_FAILURE);
  }

  opt::index_path = argv[optind++];
  opt::fx_path = argv[optind++];
}

#endif
