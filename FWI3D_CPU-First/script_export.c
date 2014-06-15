#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define NMAX 1
#define TTIME 101
#define kyori 50
void check_em1d();
void export_spectrum();

int main()
{
//    check_em1d();
    export_spectrum();

}

/////////////////////////////////////
void check_em1d()
{
  int i,k;
  char str[256];
  FILE *if1;
  double dtime,amp,amp2;
// read e-field
  for(i=0;i<NMAX;i++){
      sprintf(str,"./data/rex_000_%03d.dat",i);
      if1 = fopen(str,"r");
      for(k=0;k<TTIME;k++){
          fscanf(if1, "%lf   %lf   %lf\n",&dtime, &amp, &amp2);
      }
      printf("%05d  %05d   %10.3lf   %12.5lf    %12.4lf\n",kyori*i, TTIME, dtime, amp,log10(fabs(amp)));
      fclose(if1);
  }

}
/////////////////////////////////////
void export_spectrum()
{
  int i,k;
  char str[256];
  FILE *if1;
  double dtime,amp,amp2, amp3, amp4;
// read e-field
  for(i=0;i<NMAX;i++){
      if1 = fopen("./data/spectrum_g.dat","r");
      for(k=0;k<TTIME;k++){
          fscanf(if1, "%lf   %lf   %lf   %lf    %lf\n",&dtime, &amp, &amp2, &amp3, &amp4);
      }
      printf("%05d  %10.3e   %12.5e    %12.4e    %12.4e\n", \
              kyori*i, dtime, amp, amp2, amp3);
      fclose(if1);
  }

}
