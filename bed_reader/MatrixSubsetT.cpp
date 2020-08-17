#include "MatrixSubsetT.h"
#include <iostream>
#include <stdio.h>
#include <math.h> 
#include <stdlib.h>

using namespace std;

void SUFFIX(matrixSubset)(REALIN* in_, int in_iid_count, int in_sid_count, int val_shape, std::vector<size_t> iid_index, std::vector<int> sid_index, REALOUT* out)
{
	uint64_t_ out_iid_count = iid_index.size();
	uint64_t_ out_sid_count = sid_index.size();


#ifdef ORDERFIN

	//f-order in "most rapidly changing variable first"
	for (size_t val_index = 0; val_index != (size_t) val_shape; val_index++) {
		for (size_t sid_index_out = 0; sid_index_out != out_sid_count; sid_index_out++) {
			int sid_index_in = sid_index[sid_index_out];
			REALIN* in2 = in_ + in_iid_count * (uint64_t_)sid_index_in + in_iid_count * in_sid_count * (uint64_t_)val_index;

#ifdef ORDERFOUT //fin,fout
			REALOUT* out2 = out + out_iid_count * (uint64_t_)sid_index_out + out_iid_count * out_sid_count * val_index;
#else            //fin,cout
			REALOUT* out2 = out + sid_index_out * val_shape + val_index;

#endif
			for (size_t iid_index_out = 0; iid_index_out != out_iid_count; iid_index_out++) {
				size_t iid_index_in = iid_index[iid_index_out];

#ifdef ORDERFOUT //fin,fout
				out2[iid_index_out] = (REALOUT)in2[iid_index_in];
#else            //fin,cout
				out2[out_sid_count * val_shape * (uint64_t_)iid_index_out] = (REALOUT)in2[iid_index_in];
#endif
			}
		}
	}

#else
	//c-order in "most rapidly changing index is last"
	for (size_t iid_index_out = 0; iid_index_out != out_iid_count; iid_index_out++) {
		size_t iid_index_in = iid_index[iid_index_out];

		REALIN* in2 = in_ + in_sid_count * val_shape * (uint64_t_)iid_index_in;

#ifdef ORDERFOUT //cin,fout
		REALOUT* out2 = out + iid_index_out;
#else            //cin,cout
		REALOUT* out2 = out + out_sid_count * val_shape * (uint64_t_)iid_index_out;
#endif

		for (size_t sid_index_out = 0; sid_index_out != out_sid_count; sid_index_out++) {
			int sid_index_in = sid_index[sid_index_out];
			REALIN* in3 = in2 + val_shape * (uint64_t_)sid_index_in;
			for (size_t val_index = 0; val_index != (size_t) val_shape; val_index++) {
#ifdef ORDERFOUT //cin,fout
				out2[out_iid_count * (uint64_t_)sid_index_out + out_iid_count * out_sid_count * val_index] = (REALOUT)in3[val_index];
#else            //cin,cout
				out2[sid_index_out * val_shape + val_index] = (REALOUT)in3[val_index];
#endif
			}
		}
	}
#endif
}