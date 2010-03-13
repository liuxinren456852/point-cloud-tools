#include "stdafx.h"
#include "define.h"

// No padding
static_assert(sizeof(float3) == 3*sizeof(float), "wrong size");
static_assert(sizeof(double3) == 3*sizeof(double), "wrong size");

// No end-of-struct padding
typedef Point3<float> test_vec2f[2];
typedef Point3<double> test_vec2d[2];
static_assert(sizeof(test_vec2f) == 2*3*sizeof(float), "wrong size");
static_assert(sizeof(test_vec2d) == 2*3*sizeof(double), "wrong size");

// Natural order preserved 
static_assert(offsetof(float3,x)==0, "wrong offset");
static_assert(offsetof(float3,y)==sizeof(float), "wrong offset");
static_assert(offsetof(float3,z)==2*sizeof(float), "wrong offset");

static_assert(offsetof(double3,x)==0, "wrong offset");
static_assert(offsetof(double3,y)==sizeof(double), "wrong offset");
static_assert(offsetof(double3,z)==2*sizeof(double), "wrong offset");
