#ifndef _GPU_COMMON_H_
#define _GPU_COMMON_H_

// #define DEBUG

#define ERR_EQ(X,Y) do { if ((X) == (Y)) { \
    fprintf(stderr,"Error in %s at %s:%d\n",__func__,__FILE__,__LINE__); \
    exit(-1);}} while(0)

#define ERR_NE(X,Y) do { if ((X) != (Y)) { \
    fprintf(stderr,"Error in %s at %s:%d\n",__func__,__FILE__,__LINE__); \
    exit(-1);}} while(0)

#define CUDA_CALL(X) ERR_NE((X),cudaSuccess)
#define CUDA_ERROR() \
    do {								\
	cudaError_t cerr = cudaGetLastError();				\
	if (cerr != cudaSuccess) {					\
	    const char *ptr = cudaGetErrorString( cerr );			\
	    fprintf(stderr,"Error in %s at %s:%d, %s\n",__func__,__FILE__,__LINE__,ptr); \
	    exit(-1);							\
	}								\
    } while (0);

#define DIVCEIL(X,Y)  ((X+Y-1)/Y)

#define CALIDX2(Y,NY,Z,NZ)  ((Y)+((NY)*(Z)))
#define CALIDX3(X,NX,Y,NY,Z,NZ)  CALIDX2(X,NX, CALIDX2(Y,NY,Z,NZ),(NY)*(NZ))
#define CALIDX4(W,NW,X,NX,Y,NY,Z,NZ)  CALIDX3(W,NW,X,NX, CALIDX2(Y,NY,Z,NZ),(NY)*(NZ))

#define MAX(X,Y) ((X)>(Y) ? (X) : (Y))
#define MIN(X,Y) ((X)<(Y) ? (X) : (Y))
#define ABS(X,Y) ((X)>(Y) ? (X)-(Y) : (Y)-(X))
#define NINT(X)  ((X)>=0 ? (int)((X)+0.5) : (int)((X)-0.5))

#endif // _GPU_COMMON_H_
