#define REAL double
#define ORDERC
#undef ORDERF
#undef MISSING_VALUE
#define SUFFIX(NAME) NAME ## doubleCAAA
#include "CPlinkBedFileT.cpp"
#undef REAL
#undef SUFFIX

#define REAL float
#define ORDERC
#undef ORDERF
#undef MISSING_VALUE
#define SUFFIX(NAME) NAME ## floatCAAA
#include "CPlinkBedFileT.cpp"
#undef REAL
#undef SUFFIX

#define REAL double
#define ORDERF
#undef ORDERC
#undef MISSING_VALUE
#define SUFFIX(NAME) NAME ## doubleFAAA
#include "CPlinkBedFileT.cpp"
#undef REAL
#undef SUFFIX

#define REAL float
#define ORDERF
#undef ORDERC
#undef MISSING_VALUE
#define SUFFIX(NAME) NAME ## floatFAAA
#include "CPlinkBedFileT.cpp"
#undef REAL
#undef SUFFIX

#define REAL signed char
#define ORDERF
#undef ORDERC
#define MISSING_VALUE -127
#define SUFFIX(NAME) NAME ## int8FAAA
#include "CPlinkBedFileT.cpp"
#undef REAL
#undef SUFFIX

#define REAL signed char
#define ORDERC
#undef ORDERF
#define MISSING_VALUE -127
#define SUFFIX(NAME) NAME ## int8CAAA
#include "CPlinkBedFileT.cpp"
#undef REAL
#undef SUFFIX
