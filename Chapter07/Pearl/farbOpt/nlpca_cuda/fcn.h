#define N_INPUT (2)
#define N_H1 (10)
#define N_H2 (1)
#define N_H3 (10)
#define N_OUTPUT (0)
#define EXAMPLE_SIZE (2)
inline float myFunc(const int vIndex, const float  * restrict P, 
                    const float * restrict ex, const int nExamples,
                    float * restrict pred)
{
   float in[2];
   in[0] = ex[IN(0, nExamples, vIndex)];
   in[1] = ex[IN(1, nExamples, vIndex)];
   register float h1_0 = P[0];
   register float h1_1 = P[1];
   register float h1_2 = P[2];
   register float h1_3 = P[3];
   register float h1_4 = P[4];
   register float h1_5 = P[5];
   register float h1_6 = P[6];
   register float h1_7 = P[7];
   register float h1_8 = P[8];
   register float h1_9 = P[9];
   h1_0 += in[0] * P[10];
   h1_1 += in[0] * P[11];
   h1_2 += in[0] * P[12];
   h1_3 += in[0] * P[13];
   h1_4 += in[0] * P[14];
   h1_5 += in[0] * P[15];
   h1_6 += in[0] * P[16];
   h1_7 += in[0] * P[17];
   h1_8 += in[0] * P[18];
   h1_9 += in[0] * P[19];
   h1_0 += in[1] * P[20];
   h1_1 += in[1] * P[21];
   h1_2 += in[1] * P[22];
   h1_3 += in[1] * P[23];
   h1_4 += in[1] * P[24];
   h1_5 += in[1] * P[25];
   h1_6 += in[1] * P[26];
   h1_7 += in[1] * P[27];
   h1_8 += in[1] * P[28];
   h1_9 += in[1] * P[29];
   h1_0 = G(h1_0);
   h1_1 = G(h1_1);
   h1_2 = G(h1_2);
   h1_3 = G(h1_3);
   h1_4 = G(h1_4);
   h1_5 = G(h1_5);
   h1_6 = G(h1_6);
   h1_7 = G(h1_7);
   h1_8 = G(h1_8);
   h1_9 = G(h1_9);
   register float h2_0 = P[30];
   h2_0 += h1_0 * P[31];
   h2_0 += h1_1 * P[32];
   h2_0 += h1_2 * P[33];
   h2_0 += h1_3 * P[34];
   h2_0 += h1_4 * P[35];
   h2_0 += h1_5 * P[36];
   h2_0 += h1_6 * P[37];
   h2_0 += h1_7 * P[38];
   h2_0 += h1_8 * P[39];
   h2_0 += h1_9 * P[40];
   register float h3_0 = P[41];
   register float h3_1 = P[42];
   register float h3_2 = P[43];
   register float h3_3 = P[44];
   register float h3_4 = P[45];
   register float h3_5 = P[46];
   register float h3_6 = P[47];
   register float h3_7 = P[48];
   register float h3_8 = P[49];
   register float h3_9 = P[50];
   h3_0 += h2_0 * P[51];
   h3_1 += h2_0 * P[52];
   h3_2 += h2_0 * P[53];
   h3_3 += h2_0 * P[54];
   h3_4 += h2_0 * P[55];
   h3_5 += h2_0 * P[56];
   h3_6 += h2_0 * P[57];
   h3_7 += h2_0 * P[58];
   h3_8 += h2_0 * P[59];
   h3_9 += h2_0 * P[60];
   h3_0 = G(h3_0);
   h3_1 = G(h3_1);
   h3_2 = G(h3_2);
   h3_3 = G(h3_3);
   h3_4 = G(h3_4);
   h3_5 = G(h3_5);
   h3_6 = G(h3_6);
   h3_7 = G(h3_7);
   h3_8 = G(h3_8);
   h3_9 = G(h3_9);
   register float o,sum = 0.f;
   o = P[61];
   o += h3_0 * P[62];
   o += h3_1 * P[63];
   o += h3_2 * P[64];
   o += h3_3 * P[65];
   o += h3_4 * P[66];
   o += h3_5 * P[67];
   o += h3_6 * P[68];
   o += h3_7 * P[69];
   o += h3_8 * P[70];
   o += h3_9 * P[71];
#ifdef DO_PRED
   pred[0] = o;
#endif
   o -= in[0];
   sum += o*o;
   o = P[72];
   o += h3_0 * P[73];
   o += h3_1 * P[74];
   o += h3_2 * P[75];
   o += h3_3 * P[76];
   o += h3_4 * P[77];
   o += h3_5 * P[78];
   o += h3_6 * P[79];
   o += h3_7 * P[80];
   o += h3_8 * P[81];
   o += h3_9 * P[82];
#ifdef DO_PRED
   pred[1] = o;
#endif
   o -= in[1];
   sum += o*o;
   return(sum);
}

#define N_PARAM (83)
#define FLOP_ESTIMATE (128 + 20 * G_ESTIMATE)
