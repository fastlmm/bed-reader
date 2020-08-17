import numpy as np

cimport numpy as np
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool

cdef extern from "./MatrixSubset.h":
	

	void _matrixSubsetDoubleFToDoubleFAAA "matrixSubsetDoubleFToDoubleFAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _matrixSubsetDoubleFToDoubleCAAA "matrixSubsetDoubleFToDoubleCAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _matrixSubsetDoubleCToDoubleFAAA "matrixSubsetDoubleCToDoubleFAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _matrixSubsetDoubleCToDoubleCAAA "matrixSubsetDoubleCToDoubleCAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)

	void _matrixSubsetDoubleFToSingleFAAA "matrixSubsetDoubleFToSingleFAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _matrixSubsetDoubleFToSingleCAAA "matrixSubsetDoubleFToSingleCAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _matrixSubsetDoubleCToSingleFAAA "matrixSubsetDoubleCToSingleFAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _matrixSubsetDoubleCToSingleCAAA "matrixSubsetDoubleCToSingleCAAA"(double* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)

	void _matrixSubsetSingleFToDoubleFAAA "matrixSubsetSingleFToDoubleFAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _matrixSubsetSingleFToDoubleCAAA "matrixSubsetSingleFToDoubleCAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _matrixSubsetSingleCToDoubleFAAA "matrixSubsetSingleCToDoubleFAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)
	void _matrixSubsetSingleCToDoubleCAAA "matrixSubsetSingleCToDoubleCAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, double* out)

	void _matrixSubsetSingleFToSingleFAAA "matrixSubsetSingleFToSingleFAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _matrixSubsetSingleFToSingleCAAA "matrixSubsetSingleFToSingleCAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _matrixSubsetSingleCToSingleFAAA "matrixSubsetSingleCToSingleFAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)
	void _matrixSubsetSingleCToSingleCAAA "matrixSubsetSingleCToSingleCAAA"(float* in_, int input_num_ind, int input_num_snps, int input_num_vals, vector[size_t] iid_idx_list, vector[int] sid_idx_list, float* out)


def matrixSubsetDoubleFToDoubleFAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleFToDoubleFAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out
def matrixSubsetDoubleFToDoubleCAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleFToDoubleCAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out
def matrixSubsetDoubleCToDoubleFAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleCToDoubleFAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out
def matrixSubsetDoubleCToDoubleCAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleCToDoubleCAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out

def matrixSubsetDoubleFToSingleFAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleFToSingleFAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
def matrixSubsetDoubleFToSingleCAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleFToSingleCAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
def matrixSubsetDoubleCToSingleFAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleCToSingleFAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
def matrixSubsetDoubleCToSingleCAAA(np.ndarray[np.float64_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetDoubleCToSingleCAAA(<double*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out

def matrixSubsetSingleFToDoubleFAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleFToDoubleFAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out
def matrixSubsetSingleFToDoubleCAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleFToDoubleCAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out
def matrixSubsetSingleCToDoubleFAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleCToDoubleFAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out
def matrixSubsetSingleCToDoubleCAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float64_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleCToDoubleCAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <double*> out.data)
	return out

def matrixSubsetSingleFToSingleFAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleFToSingleFAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
def matrixSubsetSingleFToSingleCAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleFToSingleCAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
def matrixSubsetSingleCToSingleFAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleCToSingleFAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
def matrixSubsetSingleCToSingleCAAA(np.ndarray[np.float32_t, ndim=3] in_, input_num_ind, input_num_snps, input_num_vals, iidIdxList, snpIdxList, np.ndarray[np.float32_t, ndim=3] out):
	cdef vector[size_t] iid_idx_list = iidIdxList
	cdef vector[int] sid_idx_list = snpIdxList
	_matrixSubsetSingleCToSingleCAAA(<float*> in_.data, input_num_ind, input_num_snps, input_num_vals, iid_idx_list, sid_idx_list, <float*> out.data)
	return out
