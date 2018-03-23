/*
	Copyright 2011-2018 Daniel S. Buckstein

	Licensed under the Apache License, Version 2.0 (the "License");
	you may not use this file except in compliance with the License.
	You may obtain a copy of the License at

		http://www.apache.org/licenses/LICENSE-2.0

	Unless required by applicable law or agreed to in writing, software
	distributed under the License is distributed on an "AS IS" BASIS,
	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
	See the License for the specific language governing permissions and
	limitations under the License.
*/

/*
	animal3D SDK: Minimal 3D Animation Framework
	animal3D Math (A3DM) SDK
	By Daniel S. Buckstein

	a3stats.h
	Declarations for miscellaneous statistical functions.
*/

#ifndef __ANIMAL3D_A3DM_STATS_H
#define __ANIMAL3D_A3DM_STATS_H


#include "animal3D/a3/a3config.h"
#include "animal3D/a3/a3types_integer.h"
#include "animal3D/a3/a3types_real.h"


A3_BEGIN_DECL


//-----------------------------------------------------------------------------

// median, mean, variance, standard deviation (with options to return mean)
A3_INLINE a3real a3median(const a3real data[], const a3index n);
A3_INLINE a3real a3medianInt(const a3integer data[], const a3index n);

A3_INLINE a3real a3mean(const a3real data[], const a3index n);
A3_INLINE a3real a3meanInt(const a3integer data[], const a3index n);

A3_INLINE a3real a3variance(const a3real data[], const a3index n, a3real *mean_out);
A3_INLINE a3real a3varianceInt(const a3integer data[], const a3index n, a3real *mean_out);

A3_INLINE a3real a3standardDeviation(const a3real data[], const a3index n, a3real *mean_out);
A3_INLINE a3real a3standardDeviationInt(const a3integer data[], const a3index n, a3real *mean_out);

// factorial, permutations and combinations
A3_INLINE a3bigindex a3factorial(a3bigindex n);
A3_INLINE a3bigindex a3permutations(const a3bigindex n, const a3bigindex k);
A3_INLINE a3bigindex a3combinations(const a3bigindex n, const a3bigindex k);


//-----------------------------------------------------------------------------


A3_END_DECL


#ifdef A3_OPEN_SOURCE
#include "_inl/a3stats_impl.inl"
#endif	// A3_OPEN_SOURCE

#endif	// !__ANIMAL3D_A3DM_STATS_H