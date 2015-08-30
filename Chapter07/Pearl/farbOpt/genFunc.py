#!/usr/bin/python
# Rob Farber
nInput=2
nH1=10
nH2=1
nH3=10

print "#define N_INPUT (%d)" % (nInput) 
print "#define N_H1 (%d)" % (nH1) 
print "#define N_H2 (%d)" % (nH2) 
print "#define N_H3 (%d)" % (nH3) 
print "#define N_OUTPUT (%d)" % (0) 
print "#define EXAMPLE_SIZE (%d)" % (nInput) 

index=0
flopEstimate=0;
gCalls=0;

#print start of inline function
print "inline float myFunc(const int vIndex, const float  * restrict P, " 
print "                    const float * restrict ex, const int nExamples,"
print "                    float * restrict pred)"
print "{"
print "   float in[%d];" % (nInput)

for i in range(0,nInput):
    print "   in[%d] = ex[IN(%d, nExamples, vIndex)];" % (i,i)

for i in range(0,nH1):
   print "   register float h1_%d = P[%d];" % (i,index)
   index += 1

#input to h1
for i in range(0,nInput):
    for j in range(0,nH1):
        print "   h1_%d += in[%d] * P[%d];" % (j,i,index)
        index += 1
        flopEstimate += 2

for j in range(0,nH1):
    print "   h1_%d = G(h1_%d);" % (j,j)
    gCalls += 1
   
for i in range(0,nH2):
   print "   register float h2_%d = P[%d];" % (i,index)
   index += 1
   
for i in range(0,nH1):
    for j in range(0,nH2):
        print "   h2_%d += h1_%d * P[%d];" % (j,i,index)
        index += 1
        flopEstimate += 2

for i in range(0,nH3):
   print "   register float h3_%d = P[%d];" % (i,index)
   index += 1
   
for i in range(0,nH2):
    for j in range(0,nH3):
        print "   h3_%d += h2_%d * P[%d];" % (j,i,index)
        index += 1
        flopEstimate += 2

for j in range(0,nH3):
    print "   h3_%d = G(h3_%d);" % (j,j)
    gCalls += 1
   
print "   register float o,sum = 0.f;"

for i in range(0,nInput):
    print "   o = P[%d];" % (index)
    index += 1
    for j in range(0,nH3):
        print "   o += h3_%d * P[%d];" % (j,index)
        index += 1
        flopEstimate += 2

    print "#ifdef DO_PRED"
    print "   pred[%d] = o;" %(i)
    print "#endif"
    print "   o -= in[%d];" % (i)
    print "   sum += o*o;"
    flopEstimate += 3

print "   return(sum);"
print "}"
print 
print "#define N_PARAM (%d)" % (index) 
flopEstimate += 2 # to account for the square and global sum
print "#define FLOP_ESTIMATE (%d + %d * G_ESTIMATE)" % (flopEstimate, gCalls) 
        
