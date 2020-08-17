#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToDoubleFAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToDoubleCAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToDoubleFAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT double
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToDoubleCAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToSingleFAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleFToSingleCAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToSingleFAAA
#include "MatrixSubsetT.cpp"


#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN double
#define REALOUT float
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## DoubleCToSingleCAAA
#include "MatrixSubsetT.cpp"

// REALIN float

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToDoubleFAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToDoubleCAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToDoubleFAAA
#include "MatrixSubsetT.cpp"


#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT double
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToDoubleCAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#define ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToSingleFAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#define ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleFToSingleCAAA
#include "MatrixSubsetT.cpp"

#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#undef ORDERFIN
#define ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToSingleFAAA
#include "MatrixSubsetT.cpp"


#undef REALIN
#undef REALOUT
#undef SUFFIX
#define REALIN float
#define REALOUT float
#undef ORDERFIN
#undef ORDERFOUT
#define SUFFIX(NAME) NAME ## SingleCToSingleCAAA
#include "MatrixSubsetT.cpp"
