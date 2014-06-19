#include<stdlib.h>
#include<stdio.h>
#include<math.h>

int main(){
    double i;
    double tmp;

    for(i=-100;i<100;i=i+0.1){
        tmp = sin(M_PI*i)/M_PI/i;
        printf("%lf\n",tmp);
    }
}

