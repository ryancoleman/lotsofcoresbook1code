// This code supplements the white paper
//    "Multithreaded Transposition of Square Matrices
//     with Common Code for 
//     Intel Xeon Processors and Intel Xeon Phi Coprocessors"
// available at the following URL:
//     http://research.colfaxinternational.com/post/2013/08/12/Trans-7110.aspx
// You are free to use, modify and distribute this code as long as you acknowledge
// the above mentioned publication.
// (c) Colfax International, 2013

#ifndef __INCLUDED_TRANSPOSE_H__
#define __INCLUDED_TRANSPOSE_H__

// Allow compile with single or double precision
// by specifying the compiler flag -DSINGLE or -DDOUBLE
#ifdef SINGLE
#define FTYPE float
#elif defined DOUBLE
#define FTYPE double
#endif

void Transpose(FTYPE* const A, const int n);

#endif
