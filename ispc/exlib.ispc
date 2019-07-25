/*
 *  Copyright (c) 2003-2010, Mark Borgerding. All rights reserved.
 *  This file is part of KISS FFT - https://github.com/mborgerding/kissfft
 *
 *  SPDX-License-Identifier: BSD-3-Clause
 *  See COPYING file for more information.
 */
inline void aos_to_soa2(const uniform float *uniform SrcPtr, 
                varying float * uniform v0,
                varying float * uniform v1)
{
	// Convert to varying using the shuffle functions

#if TARGET_WIDTH == 4
	static const varying int vPermr = {0, 2, 4, 6};
	static const varying int vPermi = {1, 3, 5, 7};
#elif TARGET_WIDTH == 8
	static const varying int vPermr = {0, 2, 4, 6, 8, 10, 12, 14};
	static const varying int vPermi = {1, 3, 5, 7, 9, 11, 13, 15};
#elif TARGET_WIDTH == 16
	// perf seems whacky - need to verify
	static const varying int vPermr = {0, 2, 4, 6, 8, 10, 12, 14,
    									16, 18, 20, 22, 24, 26, 28, 30};
	static const varying int vPermi = {1, 3, 5, 7, 9, 11, 13, 15,
										17, 19, 21, 23, 25, 27, 29, 31};
#else
	#error "Unknown ISPC_TARGET_* value"
#endif

	// Load in 2x gangs of floats.
	// s1 = 0.r, 0.i, 1.r, 1.i, 2.r, 2.i, 3.r, 3.i
	// s2 = 4.r, 4.i, 5.r, 5.i, 6.r, 6.i, 7.r, 7.i
	const float S1 = SrcPtr[programIndex];
	const float S2 = SrcPtr[programIndex+programCount];

	// Shuffle inputs to get
	//   r = 0.r, 1.r, 2.r, 3.r, 4.r, 5.r, 6.r, 7.r
	//   i = 0.i, 1.i, 2.i, 3.i, 4.i, 5.i, 6.i, 7.i
	*v0 = shuffle(S1, S2, vPermr);
	*v1 = shuffle(S1, S2, vPermi);
}

inline void soa_to_aos2(const varying float &S1,
                 const varying float &S2,
                 uniform float *uniform DstPtr)
{
	// Convert from varying using the shuffle functions

#if TARGET_WIDTH == 4
	static const varying int vPermr = {0, 4, 1, 5};
	static const varying int vPermi = {2, 6, 3, 7};
#elif TARGET_WIDTH == 8
	static const varying int vPermr = {0, 8, 1, 9, 2, 10, 3, 11};
	static const varying int vPermi = {4, 12, 5, 13, 6, 14, 7, 15};
#elif TARGET_WIDTH == 16
	// perf seems whacky - need to verify
	static const varying int vPermr = {0, 16, 1, 17, 2, 18, 3, 19,
    									4, 20, 5, 21, 6, 22, 7, 23};
	static const varying int vPermi = {8, 24, 9, 25, 10, 26, 11, 27,
										12, 28, 13, 29, 14, 30, 15, 31};
#else
	#error "Unknown ISPC_TARGET_* value"
#endif

	// Begin with 2x gangs of floats.
	//   r = 0.r, 1.r, 2.r, 3.r, 4.r, 5.r, 6.r, 7.r
	//   i = 0.i, 1.i, 2.i, 3.i, 4.i, 5.i, 6.i, 7.i

	// Shuffle inputs to get
	// s1 = 0.r, 0.i, 1.r, 1.i, 2.r, 2.i, 3.r, 3.i
	// s2 = 4.r, 4.i, 5.r, 5.i, 6.r, 6.i, 7.r, 7.i
	DstPtr[programIndex] = shuffle(S1, S2, vPermr);
	DstPtr[programIndex+programCount] = shuffle(S1, S2, vPermi);
}
