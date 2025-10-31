// simdutf_wrapper.cpp - Separate compilation unit for simdutf implementation
// Disable advanced SIMD features that might cause issues
#define SIMDUTF_IMPLEMENTATION_ICELAKE 0
#define SIMDUTF_IMPLEMENTATION_HASWELL 0  
// Rick: disable as I'm not sure how to properly enable SSE etc with R packages
#define SIMDUTF_IMPLEMENTATION_WESTMERE 0
//#define SIMDUTF_IMPLEMENTATION_WESTMERE 1  // Only use basic SSE2/SSSE3
#define SIMDUTF_IMPLEMENTATION_ARM64 0
#define SIMDUTF_IMPLEMENTATION_PPC64 0
#define SIMDUTF_IMPLEMENTATION_RVV 0
#define SIMDUTF_IMPLEMENTATION_FALLBACK 1

#include "simdutf/simdutf.cpp"
