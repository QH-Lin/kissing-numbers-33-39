/* Independent check of the equatorial lines in dimension 38, written from scratch (no shared code with
   d38_verify.py): every line is in the Leech lattice defined by the given Golay basis (integer model with minimum
   squared norm 32), has squared norm 48, the lines are pairwise distinct with |v.w| <= 24, and |v.u| <= 24 for all
   196560 minimal vectors u (generated here from the Golay code). Together with the base configuration this gives
   cos <= 1/2 for every pair involving a new point (see README).
   build: cc -O3 -o d38_check_equator verify/d38_check_equator.c
   usage: ./d38_check_equator certificates/d38/base-data/golay_basis.txt certificates/d38/equator_lines_6855.txt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
static uint32_t code[4096];
static int member(const int*v){ int e=((v[0]%2)+2)%2; uint32_t w=0; long s=0; for(int i=0;i<24;i++){ if((((v[i]%2)+2)%2)!=e) return 0; int t=(v[i]-e)/2; if(((t%2)+2)%2) w|=1u<<i; s+=v[i]; }
  int ok=0; for(int k=0;k<4096;k++) if(code[k]==w) ok=1; return ok && (((s%8)+8)%8)==4*e; }
int main(int argc,char**argv){
  FILE*f=fopen(argv[1],"r"); char line[256]; uint32_t B[12]; int nb=0;
  while(fgets(line,sizeof line,f)){ uint32_t m=0; int k=0; for(char*p=line;*p;p++) if(*p=='0'||*p=='1'){ if(*p=='1') m|=1u<<k; k++; } if(k==24) B[nb++]=m; } fclose(f);
  for(int s=0;s<4096;s++){ uint32_t c=0; for(int b=0;b<12;b++) if(s>>b&1) c^=B[b]; code[s]=c; }
  static int V[200000][24]; int n=0; f=fopen(argv[2],"r");
  while(fgets(line,sizeof line,f)){ int k=0; char*p=line; while(k<24){ char*e; long x=strtol(p,&e,10); if(e==p) break; V[n][k++]=x; p=e; } if(k==24) n++; } fclose(f);
  long bad=0; for(int i=0;i<n;i++){ int s=0; for(int k=0;k<24;k++) s+=V[i][k]*V[i][k]; if(s!=48||!member(V[i])) bad++; }
  printf("%d lines; membership/norm failures: %ld\n",n,bad);
  int mx=0; long pbad=0; for(int i=0;i<n;i++) for(int j=i+1;j<n;j++){ int s=0; for(int k=0;k<24;k++) s+=V[i][k]*V[j][k]; if(s<0) s=-s; if(s>mx) mx=s; if(s>24) pbad++; }
  printf("pairwise max |v.w| = %d (limit 24, 48 would mean a repeated line), violations %ld\n",mx,pbad);
  /* minimal shell (norm 32): (+-2^8 on octads, even # of minus), (-+3, +-1^23) odd type, (+-4^2) */
  static int S[200000][24]; int ns=0; int v[24];
  for(int k=0;k<4096;k++){ uint32_t c=code[k]; if(__builtin_popcount(c)!=8) continue; int pos[8],np=0; for(int i=0;i<24;i++) if(c>>i&1) pos[np++]=i;
    for(int sg=0;sg<256;sg++){ memset(v,0,sizeof v); for(int i=0;i<8;i++) v[pos[i]]=(sg>>i&1)?-2:2; if(member(v)) memcpy(S[ns++],v,sizeof v); } }
  for(int a=0;a<24;a++) for(int b=a+1;b<24;b++) for(int s=0;s<4;s++){ memset(v,0,sizeof v); v[a]=(s&1)?-4:4; v[b]=(s&2)?-4:4; if(member(v)) memcpy(S[ns++],v,sizeof v); }
  for(int k=0;k<4096;k++){ uint32_t c=code[k]; for(int a=0;a<24;a++){ for(int i=0;i<24;i++){ int inc=c>>i&1; v[i]= (i==a) ? (inc?3:-3) : (inc?-1:1); } if(member(v)) memcpy(S[ns++],v,sizeof v); } }
  printf("minimal shell vectors generated: %d (expected 196560)\n",ns);
  int mx2=0; long cbad=0; for(int i=0;i<n;i++) for(int j=0;j<ns;j++){ int s=0; for(int k=0;k<24;k++) s+=V[i][k]*S[j][k]; if(s<0) s=-s; if(s>mx2) mx2=s; if(s>24) cbad++; }
  printf("new-vs-minimal-shell max |v.u| = %d (limit 24), violations %ld\n",mx2,cbad);
  int fail = bad||pbad||cbad||ns!=196560; printf(fail ? "CHECK FAILED\n" : "ALL CHECKS PASSED\n"); return fail; }
