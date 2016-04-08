//------------------------------------------------------------------------------
//
//  PROGRAM: Matrix Multipliplication blocked
//
//  PURPOSE: This is a simple matrix multiplication program
//
//                C  = A * B
//
//           A and B are set to constant matrices so we
//           can make a quick test of the multiplication.
//
//           A work-group updates a block of the product matrix
//           from blocks of A and B copied into local memory.
//
//  USAGE:   The matrices are constant matrices, square and the order is
//           set as a constant, ORDER (see mult.h).
//
//  HISTORY: Written by Tim Mattson using an algorithm suggested  
//           by Tom Deakin and Simon McIntosh-Smith.  Tom also
//           provided lots of debugging help, November 2013 
//
//  LICENSE: This work is licensed under the Creative Commons
//           Attribution 4.0 International License.
//           To view a copy of this license, visit
//           http://creativecommons.org/licenses/by/4.0/
//           or send a letter to:
//              Creative Commons,
//              444 Castro Street, Suite 900,
//              Mountain View, California, 94041, USA.
//
//------------------------------------------------------------------------------

#include "matmul.hpp"
#include "matmul_lib.hpp"

// pick up device type from compiler command line or from the default type
#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_DEFAULT
#endif

int main(void)
{

    int N;   // A[N][N], B[N][N], C[N][N]
    int sz;  // number of elements in each matrix
    float tmp;

    N = ORDER;

    sz = N * N;

    std::vector<float> h_A(sz); // Matrix A on the host
    std::vector<float> h_B(sz); // Matrix B on the host
    std::vector<float> h_C(sz); // Matrix C on the host

    cl::Buffer d_A;    // matrix A on the device
    cl::Buffer d_B;    // matrix B on the device
    cl::Buffer d_C;    // matrix C on the device

    initmat(N, N, N, h_A, h_B, h_C);

    printf("\n===== Sequential, matrix mult (dot prod), order %d on CPU ======\n",ORDER);
 
    zero_mat(N, N, h_C);

    util::Timer timer;

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            tmp = 0.0f;
            for (int k = 0; k < N; k++) {
                tmp += h_A[i*N+k] * h_B[k*N+j];
            }
            h_C[i*N+j] = tmp;
        }
    }
              
    double rtime = static_cast<double>(timer.getTimeMilliseconds()) / 1000.0;

    results(N, N, N, h_C, rtime);

    printf("\n===== Parallel matrix mult (dot prod), order %d on device ======\n",ORDER);

    switch (DEVICE) {
      case CL_DEVICE_TYPE_DEFAULT: printf("DEVICE=DEFAULT\n"); break;
      case CL_DEVICE_TYPE_CPU:     printf("DEVICE=CPU\n"); break;
      case CL_DEVICE_TYPE_GPU:     printf("DEVICE=GPU\n"); break;
      default:                     printf("DEVICE=%d\n", DEVICE); break;
    }
 
    zero_mat(N, N, h_C);
    try
    {
   
       cl::Context context(DEVICE);

       // Load in kernel source, creating a program object for the context.
       // Build program explicitly so I can catch errors and display
       // compiler error messages (should any be generated)

       cl::Program program(context, util::loadProgram("matmul_kernel.cl"));
       try
       {
           program.build();
       }
       catch (cl::Error error)
       {
          // If it was a build error then show the error
          if (error.err() == CL_BUILD_PROGRAM_FAILURE)
           {
               std::vector<cl::Device> devices;
               devices = context.getInfo<CL_CONTEXT_DEVICES>();
               std::string built = program.getBuildInfo<CL_PROGRAM_BUILD_LOG>(devices[0]);
               std::cerr << built << "\n";
           }
           throw error;
       }


        // Get the command queue
        cl::CommandQueue queue(context);


        // Create the kernel functor
 
        auto mmul = cl::make_kernel<int, cl::Buffer, cl::Buffer, cl::Buffer, 
                                    cl::LocalSpaceArg, cl::LocalSpaceArg>   
                                    (program, "mmul");

        util::Timer timer;


        d_A   = cl::Buffer(context, begin(h_A), end(h_A), true);
        d_B   = cl::Buffer(context, begin(h_B), end(h_B), true);
        d_C   = cl::Buffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * sz);

        // Work-group computes a block of C.  This size is also set
        // in a #define inside the kernel function.  Note this blocksize
        // must evenly divide the matrix order
        int blocksize = 16;  

        cl::LocalSpaceArg A_block = cl::Local(sizeof(float) * blocksize*blocksize);
        cl::LocalSpaceArg B_block = cl::Local(sizeof(float) * blocksize*blocksize);
 
        mmul(
            cl::EnqueueArgs(
            queue,
            cl::NDRange(N,N),
            cl::NDRange(blocksize,blocksize)),
            N, 
            d_A,
            d_B,
            d_C,
            A_block,
            B_block);

        cl::copy(queue, d_C, begin(h_C), end(h_C));

        double rtime = static_cast<double>(timer.getTimeMilliseconds()) / 1000.0;

        results(N, N, N, h_C, rtime);
          
    }
    catch (cl::Error err) {
        std::cout << "Exception\n";
        std::cerr 
            << "ERROR: "
            << err.what()
            << std::endl;
 
    }

}
