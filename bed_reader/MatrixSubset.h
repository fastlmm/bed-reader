#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToDoubleFAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToDoubleCAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToDoubleFAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToDoubleCAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToSingleFAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToSingleCAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToSingleFAAA
#include "MatrixSubsetT.h"


#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToSingleCAAA
#include "MatrixSubsetT.h"

// REALIN float

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToDoubleFAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToDoubleCAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToDoubleFAAA
#include "MatrixSubsetT.h"


#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToDoubleCAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToSingleFAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToSingleCAAA
#include "MatrixSubsetT.h"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToSingleFAAA
#include "MatrixSubsetT.h"


#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToSingleCAAA
#include "MatrixSubsetT.h"




