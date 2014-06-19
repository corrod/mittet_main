#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#define POSIX 450.0
#define POSIY 450.0
#define POSIZ 360.0
void export_model();
void export_model_h();

int main()
{
    export_model_h();
}
/////////////////////////////////////
void export_model()
{
  int i,j,k;
  int NN=31;
  char str[256];
  FILE *if1;
  double dtime,amp,amp2, amp3, amp4;
  if1 = fopen("./old_seg/sig_015XYZ.dat","r");

  for(k=0;k<NN;k++){
    for(j=0;j<NN;j++){
      for(i=0;i<NN;i++){
          fscanf(if1,"%lf   %lf   %lf  %lf\n", &amp, &amp2, &amp3, &amp4);
          if(amp== POSIX && amp2==POSIY){
              printf("%lf   %lf\n", amp3, amp4);
          }

      }
    }
  }

}
/////////////////////////////////////
void export_model_h()
{
  int i,j,k;
  int NN=31;
  char str[256];
  FILE *if1;
  double dtime,amp,amp2, amp3, amp4;
  if1 = fopen("./old_seg/sig_015XYZ.dat","r");

  for(k=0;k<NN;k++){
    for(j=0;j<NN;j++){
      for(i=0;i<NN;i++){
          fscanf(if1,"%lf   %lf   %lf  %lf\n", &amp, &amp2, &amp3, &amp4);
          if(amp2== POSIY && amp3==POSIZ){
              printf("%lf   %lf\n", amp, amp4);
          }

      }
    }
  }

}
