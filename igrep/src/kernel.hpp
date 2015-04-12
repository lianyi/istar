#pragma once
#ifndef IGREP_KERNEL_HPP
#define IGREP_KERNEL_HPP

#define CHARACTER_CARDINALITY 4	/**< One nucleotide is either A, C, G, or T. */
#define MAX_UNSIGNED_INT 0xffffffffUL	/**< The maximum value of an unsigned int. */
#define MAX_UNSIGNED_LONG_LONG 0xffffffffffffffffULL	/**< The maximum value of an unsigned long long. */
#define B 7	/**< Each thread block consists of 2^B (=1<<B) threads. */
#define L 8	/**< Each thread processes 2^L (=1<<L) special codons plus those in the overlapping zone of two consecutive threads. */
// Since each thread block processes 1 << (L + B) special codons, the number of thread blocks will be up to (MAX_SCODON_COUNT + 1 << (L + B) - 1) >> (L + B).
// This program uses 1D CUDA thread organization, so at most 65,536 threads can be specified.
// Therefore, the inequation ((MAX_SCODON_COUNT + (1 << (L + B)) - 1) >> (L + B)) <= 65,536 must hold.
// MAX_SCODON_COUNT = 0.22G ==> L + B >= 12 is required.

/**
 * Transfer necessary parameters to CUDA constant memory.
 * This agrep kernel initialization should be called only once for searching the same corpus.
 * @param[in] scodon_arg The special codon array.
 * @param[in] character_count_arg Actual number of characters.
 * @param[in] match_arg The match array.
 * @param[in] max_match_count_arg Maximum number of matches of one single query.
 */
void initAgrepKernel(const unsigned int *scodon_arg, const unsigned int character_count_arg, const unsigned int *match_arg, const unsigned int max_match_count_arg);

/**
 * Transfer 32-bit mask array and test bit from host to CUDA constant memory.
 * @param[in] mask_array_arg The mask array of a pattern.
 * @param[in] test_bit_arg The test bit.
 */
void transferMaskArray32(const unsigned int *mask_array_arg, const unsigned int test_bit_arg);

/**
 * Transfer 64-bit mask array and test bit from host to CUDA constant memory.
 * @param[in] mask_array_arg The mask array of a pattern.
 * @param[in] test_bit_arg The test bit.
 */
void transferMaskArray64(const unsigned long long *mask_array_arg, const unsigned long long test_bit_arg);

/**
 * Invoke the CUDA implementation of agrep kernel.
 * @param[in] m Pattern length.
 * @param[in] k Edit distance.
 * @param[in] block_count Number of thread blocks.
 */
void invokeAgrepKernel(const unsigned int m, const unsigned int k, const unsigned int block_count);

/**
 * Get the number of matches from CUDA constant memory.
 * @param[out] match_count_arg Number of matches.
 */
void getMatchCount(unsigned int *match_count_arg);

#endif
