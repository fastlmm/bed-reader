/*
 * CPlinkBedFile - {PLINK Bed File Access Class}
 *
 *         File Name:   CPlinkBedFile.h
 *     Creation Date:    4 Dec 2011
 *
 *    Module Purpose:   This file defines the CPlinkBedFile class
 *
 *                      A .BED file contains compressed binary genotype values for
 *                         for individuals by SNPs.  In other contexts, we prefer and may
 *                         require the file LayoutMode be LayoutGroupGenotypesBySnp
 *
 *                      The .bed header is three bytes followed immediately by data.
 *
 *                         bedFileMagic1 | bedFileMagic2 | LayoutMode
 *                         [... data ... ]
 *
 *    Change History:   Version 2.00: Reworked to be wrapped in python version by Chris Widmer (chris@shogun-toolbox.org)
 *
 */

#include <vector>
#include <string>
#include <limits>
 //#include <inttypes.h>

using namespace std;
typedef unsigned char BYTE;
typedef unsigned long long uint64_t_;

#if !defined(CPlinkBedFileT_h_consts)
#define CPlinkBedFileT_h_consts

const BYTE bedFileMagic1 = 0x6C;       // 0b01101100 or 'l' (lowercase 'L')
const BYTE bedFileMagic2 = 0x1B;       // 0b00011011 or <esc>

enum LayoutMode
{
	LayoutUnknown = -1
	, LayoutRowMajor = 0                  // all elements of a row are sequential in memory
	, LayoutColumnMajor = 1                  // all elements of a colomn are sequential in memory
	, LayoutGroupGenotypesByIndividual = 0   // all SNP genotypes for a specific individual are seqential in memory
	, LayoutGroupGenotypesBySnp = 1          // all Individual's genotypes for a specific SNP are sequential in memory
};

enum BedGenotype     // integer representation of genotype values in Plink's binary .BED file
{
	bedHomozygousMinor = 0
	, bedMissingGenotype = 1
	, bedHeterozygous = 2
	, bedHomozygousMajor = 3
};

#endif

class SUFFIX(CBedFile)
{
public:
	SUFFIX(CBedFile)();
	SUFFIX(~CBedFile)();

	void Open(const string & filename_, size_t cIndividuals_, size_t cSnps_);                // validate the file matches the extents we expect

	// return the layout mode of the file.  
	//   We work better with files that are LayoutGroupGenotypesBySnp
	LayoutMode  GetLayoutMode();

	// return the compressed length of one line of SNP data (related to cIndividuals)
	// TODO:  MAKE THIS PRIVATE!  Consumers should be using cIndividuals or cSnps
	size_t      CbStride() { return(cbStride); }

	// return the filename associated with this CBedFile
	const string& Filename() { return(filename); }

	// read the data for one SNP (idxSnp) into the BYTE buffer pb
	size_t   ReadLine(BYTE * pb, size_t idxSnp);

	// read the genotype for all the individuals in 'list' at the SNP specified by iSNP
	void     ReadGenotypes(size_t iSnp, bool count_A1, const vector< size_t > & iIndividualList, REAL * pvOutSNP, uint64_t_ startpos, uint64_t_  outputNumSNPs);

private:
	int      NextChar();
	size_t   Read(BYTE * pb, size_t cbToRead);

	static const size_t   cbHeader = 3;         // 
	string   filename;
	FILE* pFile;
	vector< BYTE > rgBytes;
	vector< BedGenotype > rgBedGenotypes;

	LayoutMode  layout;        // 0=RowMajor(all snps per individual together);
							   // 1=ColumnMajor(all individuals per SNP together in memory)
	size_t   cIndividuals;
	size_t   cSnps;
	size_t   cbStride;


#ifdef MISSING_VALUE
	const REAL SUFFIX(unknownOrMissing) = MISSING_VALUE;
#else
	const REAL SUFFIX(unknownOrMissing) = std::numeric_limits<REAL>::quiet_NaN();  // now used by SnpInfo
#endif

	const REAL SUFFIX(homozygousPrimaryAllele) = 0;                // Major Allele
	const REAL SUFFIX(heterozygousAllele) = 1;
	const REAL SUFFIX(homozygousSecondaryAllele) = 2;              // Minor Allele ()

	const REAL SUFFIX(mapBedGenotypeToRealAlleleCountA1)[4] = {
		SUFFIX(homozygousSecondaryAllele),       // look-up 0
		SUFFIX(unknownOrMissing),                // look-up 1
		SUFFIX(heterozygousAllele),              // look-up 2
		SUFFIX(homozygousPrimaryAllele),         // look-up 3
	};

	const REAL SUFFIX(mapBedGenotypeToRealAlleleNoCountA1)[4] = {
		SUFFIX(homozygousPrimaryAllele),         // look-up 0
		SUFFIX(unknownOrMissing),                // look-up 1
		SUFFIX(heterozygousAllele),              // look-up 2
		SUFFIX(homozygousSecondaryAllele),       // look-up 3
	};
};


void SUFFIX(ImputeAndZeroMeanSNPs)(
	REAL* SNPs,
	const size_t nIndividuals,
	const size_t nSNPs,
	const bool betaNotUnitVariance,
	const REAL betaA,
	const REAL betaB,
	const bool apply_in_place,
	const bool use_stats,
	REAL* stats
	);

// to be used by cython wrapper
void SUFFIX(readPlinkBedFile)(std::string bed_fn, int inputNumIndividuals, int inputNumSNPs, bool count_A1, std::vector<size_t> individuals_idx, std::vector<int> snpIdxList, REAL* out, int num_threads);
void SUFFIX(writePlinkBedFile)(std::string bed_fn, int iid_count, int sid_count, bool count_A1, REAL* in);

/*#endif      // CPlinkBedFile_h
*/
