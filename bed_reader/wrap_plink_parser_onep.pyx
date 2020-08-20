import numpy as np 


cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "./CPlinkBedFile.h":

	void _readPlinkBedFilefloatFAAA "readPlinkBedFilefloatFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out, int num_threads)
	void _readPlinkBedFiledoubleFAAA "readPlinkBedFiledoubleFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out, int num_threads)
	void _readPlinkBedFilefloatCAAA "readPlinkBedFilefloatCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out, int num_threads)
	void _readPlinkBedFiledoubleCAAA "readPlinkBedFiledoubleCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out, int num_threads)
	void _readPlinkBedFileint8FAAA "readPlinkBedFileint8FAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, signed char* out, int num_threads)
	void _readPlinkBedFileint8CAAA "readPlinkBedFileint8CAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, vector[size_t] iid_idx_list, vector[int] sid_idx_list, signed char* out, int num_threads)

	void _writePlinkBedFilefloatFAAA "writePlinkBedFilefloatFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, float* inx)
	void _writePlinkBedFiledoubleFAAA "writePlinkBedFiledoubleFAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, double* inx)
	void _writePlinkBedFileint8FAAA "writePlinkBedFileint8FAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, signed char* inx)
	void _writePlinkBedFilefloatCAAA "writePlinkBedFilefloatCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, float* inx)
	void _writePlinkBedFiledoubleCAAA "writePlinkBedFiledoubleCAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, double* inx)
	void _writePlinkBedFileint8CAAA "writePlinkBedFileint8CAAA"(string bed_fn, int input_num_ind, int input_num_snps, bool count_A1, signed char* inx)


	void _ImputeAndZeroMeanSNPsfloatFAAA "ImputeAndZeroMeanSNPsfloatFAAA"( 
		float *SNPs,
		size_t nIndividuals,
		size_t nSNPs,
		const bool betaNotUnitVariance,
		const float betaA,
		const float betaB,
		const bool apply_in_place,
		const bool use_stats,
		float *stats
		)
	void _ImputeAndZeroMeanSNPsdoubleFAAA "ImputeAndZeroMeanSNPsdoubleFAAA"( 
		double *SNPs,
		size_t nIndividuals,
		size_t nSNPs,
		const bool betaNotUnitVariance,
		const double betaA,
		const double betaB,
		const bool apply_in_place,
		const bool use_stats,
		double *stats
		)
	void _ImputeAndZeroMeanSNPsfloatCAAA "ImputeAndZeroMeanSNPsfloatCAAA"( 
		float *SNPs,
		size_t nIndividuals,
		size_t nSNPs,
		bool betaNotUnitVariance,
		float betaA,
		float betaB,
		const bool apply_in_place,
		const bool use_stats,
		float *stats
		)

	void _ImputeAndZeroMeanSNPsdoubleCAAA "ImputeAndZeroMeanSNPsdoubleCAAA"( 
		double *SNPs,
		size_t nIndividuals,
		size_t nSNPs,
		const bool betaNotUnitVariance,
		const double betaA,
		const double betaB,
		const bool apply_in_place,
		const bool use_stats,
		double *stats
		)


def standardizefloatFAAA(np.ndarray[np.float32_t, ndim=2] out, bool betaNotUnitVariance, float betaA, float betaB, bool apply_in_place, bool use_stats, np.ndarray[np.float32_t, ndim=2] stats):
	
	num_ind = out.shape[0]
	num_snps = out.shape[1]

	#http://wiki.cython.org/tutorials/NumpyPointerToC
	_ImputeAndZeroMeanSNPsfloatFAAA(<float*> out.data, num_ind, num_snps, betaNotUnitVariance, betaA, betaB, apply_in_place, use_stats, <float *> stats.data)

	return out, stats



def standardizedoubleFAAA(np.ndarray[np.float64_t, ndim=2] out, bool betaNotUnitVariance, double betaA, double betaB, bool apply_in_place, bool use_stats, np.ndarray[np.float64_t, ndim=2] stats):
	
	num_ind = out.shape[0]
	num_snps = out.shape[1]

	#http://wiki.cython.org/tutorials/NumpyPointerToC
	_ImputeAndZeroMeanSNPsdoubleFAAA(<double*> out.data, num_ind, num_snps, betaNotUnitVariance, betaA, betaB, apply_in_place, use_stats, <double *> stats.data)

	return out, stats



def standardizefloatCAAA(np.ndarray[np.float32_t, ndim=2] out, bool betaNotUnitVariance, float betaA, float betaB, bool apply_in_place, bool use_stats, np.ndarray[np.float32_t, ndim=2] stats):
	
	num_ind = out.shape[0]
	num_snps = out.shape[1]

	#http://wiki.cython.org/tutorials/NumpyPointerToC
	_ImputeAndZeroMeanSNPsfloatCAAA(<float*> out.data, num_ind, num_snps, betaNotUnitVariance, betaA, betaB, apply_in_place, use_stats, <float *> stats.data)

	return out, stats

def standardizedoubleCAAA(np.ndarray[np.float64_t, ndim=2] out, bool betaNotUnitVariance, double betaA, double betaB,  bool apply_in_place, bool use_stats, np.ndarray[np.float64_t, ndim=2] stats):
	
	num_ind = out.shape[0]
	num_snps = out.shape[1]

	#http://wiki.cython.org/tutorials/NumpyPointerToC
	_ImputeAndZeroMeanSNPsdoubleCAAA(<double*> out.data, num_ind, num_snps, betaNotUnitVariance, betaA, betaB, apply_in_place, use_stats, <double *> stats.data)

	return out, stats

#New 
def readPlinkBedFile2floatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=2] out, num_threads):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <float*> out.data, num_threads)
	return out

def readPlinkBedFile2floatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=2] out, num_threads):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <float*> out.data, num_threads)
	return out


def readPlinkBedFile2doubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=2] out, num_threads):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <double*> out.data, num_threads)
	return out

def readPlinkBedFile2doubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=2] out, num_threads):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <double*> out.data, num_threads)
	return out

def readPlinkBedFile2int8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.int8_t, ndim=2] out, num_threads):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFileint8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <signed char*> out.data, num_threads)
	return out

def readPlinkBedFile2int8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iidIdxList, snpIdxList, np.ndarray[np.int8_t, ndim=2] out, num_threads):
	
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFileint8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <signed char*> out.data, num_threads)
	return out

def writePlinkBedFile2floatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float32_t, ndim=2] inx):
	_writePlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <float*> inx.data)

def writePlinkBedFile2floatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float32_t, ndim=2] inx):
	_writePlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <float*> inx.data)

def writePlinkBedFile2doubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float64_t, ndim=2] inx):
	_writePlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <double*> inx.data)

def writePlinkBedFile2doubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.float64_t, ndim=2] inx):
	_writePlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <double*> inx.data)

def writePlinkBedFile2int8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.int8_t, ndim=2] inx):
	_writePlinkBedFileint8FAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <signed char*> inx.data)

def writePlinkBedFile2int8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, np.ndarray[np.int8_t, ndim=2] inx):
	_writePlinkBedFileint8CAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <signed char*> inx.data)

#Old
def readPlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=2] out):
	count_A1 = False
	num_threads = 1
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <float*> out.data, num_threads)
	return out

def readPlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=2] out):
	count_A1 = False
	num_threads = 1
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <float*> out.data, num_threads)
	return out


def readPlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=2] out):
	count_A1 = False
	num_threads = 1
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <double*> out.data, num_threads)
	return out

def readPlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=2] out):
	count_A1 = False
	num_threads = 1
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	#http://wiki.cython.org/tutorials/NumpyPointerToC

	_readPlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, iid_idx_list, sid_idx_list, <double*> out.data, num_threads)
	return out


def writePlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, np.ndarray[np.float32_t, ndim=2] inx):
	count_A1 = False
	_writePlinkBedFilefloatFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <float*> inx.data)

def writePlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, np.ndarray[np.float32_t, ndim=2] inx):
	count_A1 = False
	_writePlinkBedFilefloatCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <float*> inx.data)

def writePlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, np.ndarray[np.float64_t, ndim=2] inx):
	count_A1 = False
	_writePlinkBedFiledoubleFAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <double*> inx.data)

def writePlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, np.ndarray[np.float64_t, ndim=2] inx):
	count_A1 = False
	_writePlinkBedFiledoubleCAAA(bed_fn, input_num_ind, input_num_snps, count_A1, <double*> inx.data)
