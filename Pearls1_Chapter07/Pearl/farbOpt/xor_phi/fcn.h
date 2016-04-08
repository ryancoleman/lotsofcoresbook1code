#define N_INPUT 2
#define N_OUTPUT 1
#define N_H1 1
#define EXAMPLE_SIZE  (N_INPUT+N_OUTPUT)

// Number of Parameters
// (N_OUTPUT + N_H1) // Neuron Offsets
//  + (N_INPUT*N_H1) // connections from I to H1
//  + (N_H1*N_OUTPUT) // connections from H1 to O
//  + (N_INPUT*N_OUTPUT); // connections from I to O
#define N_PARAM  ( (N_OUTPUT + N_H1) + (N_INPUT*N_H1) + (N_H1*N_OUTPUT )+ (N_INPUT*N_OUTPUT) ) 

//The final 2 ops are for the sum in the objective function
#define FLOP_ESTIMATE ( 2*(N_PARAM - N_OUTPUT-N_H1) + N_OUTPUT + 2 + N_H1*G_ESTIMATE )

inline float myFunc(const int vIndex, const float * restrict p, 
		    const float * restrict ex, const int nExamples, float * pred)
{
  register float h1;
  register float o;
  float in[2];
  
  in[0] = ex[IN(0,nExamples,vIndex)];
  in[1] = ex[IN(1,nExamples,vIndex)];

  h1 = p[0];
  o = p[1];
  h1 += in[0] * p[2];
  h1 += in[1] * p[3];
  h1 = G(h1);
  o += in[0] * p[4];
  o += in[1] * p[5];
  o += h1 * p[6];
#ifdef DO_PRED
  pred[0] = o;
#endif
  return o - ex[OUT(0,nExamples,vIndex)];
}
