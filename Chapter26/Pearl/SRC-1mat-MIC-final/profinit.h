#ifndef PROFINIT_H
#define PROFINIT_H 1

/*
#ifdef ADD_
#define profinit profinit_
#define profstart profstart_
#define profend profend_
#define profstat profstat_
#endif
*/

#ifdef __cplusplus
extern "C" {
#endif

void profinit_();
void profstart_(char *name);
void profend_(char *name);
void profstat_();


void profinit();
void profstart(char *name);
void profend(char *name);
void profstat();

#ifdef __cplusplus
};
#endif


#endif
