#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/modifier.h>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <omp.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <seqan/store.h>
#include <seqan/consensus.h>
#include <seqan/random.h>
#include <set>
#include <queue>
#include "HlaStore.h"
#include "Typing.h"
#include "options.h"

using namespace std;
using namespace seqan;


