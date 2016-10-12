/*
 *   Copyright (c) 2004-2005 Massachusetts Institute of Technology.
 *   All Rights Reserved.
 *
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *   Authors: Alexandr Andoni (andoni@mit.edu), Piotr Indyk (indyk@mit.edu)
 *	 redefine zhaoshuai 2016.10.04
*/

#ifndef __RANDOM_INCLUDED__
#define __RANDOM_INCLUDED__

#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <windows.h>
#include <assert.h>

#ifndef PI
#define PI (3.14159)
#endif

/*
#define IntT int
#define LongUns64T long long unsigned
#define Uns32T unsigned
#define Int32T int
#define BooleanT int
#define TRUE 1
#define FALSE 0
*/

#ifdef REAL_LONG_DOUBLE
#define RealT long double
#define SQRT sqrtl
#define ABS fabsl
#define LOG logl
#define COS cosl
#define FLOOR_INT32(x) ((Int32T)(floorl(x)))
#define CEIL(x) ((int)(ceill(x)))
#define POW(x, y) (powl(x, y))
#define FPRINTF_REAL(file, value) {fprintf(file, "%0.3Lf", value);}
#define FSCANF_REAL(file, value) {fscanf(file, "%Lf", value);}
#endif

#ifdef REAL_DOUBLE
#define RealT double
#define SQRT sqrt
#define ABS fabs
#define LOG log
#define COS cos
#define FLOOR_INT32(x) ((Int32T)(floor(x)))
#define CEIL(x) ((int)(ceil(x)))
#define POW(x, y) (pow(x, y))
#define FPRINTF_REAL(file, value) {fprintf(file, "%0.3lf", value);}
#define FSCANF_REAL(file, value) {fscanf(file, "%lf", value);}
#define EXP exp
#define ERF erf
#define ERFC erfc
#endif

#ifdef REAL_FLOAT
#define RealT float
#define SQRT sqrtf
#define ABS fabsf
#define LOG logf
#define COS cosf
#define FLOOR_INT32(x) ((Int32T)(floorf(x)))
#define CEIL(x) ((int)(ceilf(x)))
#define POW(x, y) (powf(x, y))
#define FPRINTF_REAL(file, value) {fprintf(file, "%0.3f", value);}
#define FSCANF_REAL(file, value) {fscanf(file, "%f", value);}
#define EXP expf
#define ERF erf
#define ERFC erfc
#endif

// The state vector for generation of random numbers.
//char rngState[256];
//#pragma comment( lib,"winmm.lib" )

namespace flann
{
	/*随机产生p稳定分布数据的数据发生器*/
	//template<typename TYPE_T>
	class RandomGen
	{
		public:

			/*构造函数，初始化随机数种子*/
			RandomGen()
			{
				//srand((unsigned)timeGetTime());
				srand((unsigned int)time(NULL));
			}
			/**/
			~RandomGen(){};

			/*
			**	Generate a random real distributed uniformly in [rangeStart,
			**	rangeEnd]. Input must satisfy: rangeStart <= rangeEnd. The
			**	granularity of generated random reals is given by RAND_MAX.
			*/
			float gen_uniform_random(float range_start, float range_end)
			{
				assert(range_start <= range_end);
				float r;
				r = (range_start + ((range_end - range_start) * (float)rand() / (float)RAND_MAX));
				assert(r >= range_start && r <= range_end);
				return(r);
			}

			/*
			**	Generate a random real from normal distribution N(0,1).
			**  Use Box-Muller transform to generate a point from normal
			**	distribution.
			**
			**	Box-Muller方法
			*/
			float gen_gaussian_random()
			{
				float x1 = 0, x2 = 0;

				while (x1 == 0)											//cannot take log of 0
				{
					x1 = gen_uniform_random(0.0, 1.0);
				}						
				x2 = gen_uniform_random(0.0,1.0);
				float z;
				z = sqrtf(-2.0 * logf(x1)) * cosf(2.0 * PI * x2);
				return(z);
			}

			/*
			**	Generate a random real from Cauchy distribution N(0,1).
			*/
			float gen_cauchy_random()
			{
				float x, y;
				x = gen_gaussian_random();
				y = gen_gaussian_random();
				if (fabsf(y) < 0.000001)
				{
					y = 0.000001;
				}
				return ( x / y);
			}

			/*
			** Generate a random 32-bits unsigned (Uns32T) in the range
			** [rangeStart, rangeEnd]. Inputs must satisfy: rangeStart <=
			** rangeEnd.
			*/
			unsigned int gen_random_uns32(unsigned int  rangeStart, unsigned int  rangeEnd)
			{
				assert(rangeStart <= rangeEnd);
				unsigned int r;
				if (RAND_MAX >= rangeEnd - rangeStart) 
				{
					r = rangeStart + (unsigned int)((rangeEnd - rangeStart + 1.0) * rand() / (RAND_MAX + 1.0));
				}
				else 
				{
					r = rangeStart + (unsigned int)( 
						(rangeEnd - rangeStart + 1.0)* 
						( (unsigned)rand() * ( (unsigned)RAND_MAX + 1) + (unsigned)rand() ) 
						/( (unsigned)RAND_MAX* ((unsigned)RAND_MAX + 1) + (unsigned)RAND_MAX + 1.0)  );
				}
				assert(r >= rangeStart && r <= rangeEnd);
				return (r);
			}


	};
}

#endif //
