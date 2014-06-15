#include <stdio.h>
#include <math.h>

#define max 10000                       //Nの限界値
#define pi 3.1415926535         //円周率

int main()
{
        int k,n,N;
        double tmp[max+1],f[max+1],ReF,ImF,om, dt, dtmp[max+1];

        FILE *function;
        FILE *fourier;
//フーリエ変換したいデータをあらかじめファイル function.data に保存しておく
        //function=fopen("./data/anal1.dat","r");
        function=fopen("./data/waveformX.dat","r");
        //function=fopen("./data/forward_x_000_000.dat","r");
        //function=fopen("./data/inverse_g.dat","r");
        //function=fopen("./data/rex_000_000.dat","r");
        fourier=fopen("./data/fourier.data","w");

//データの読み込み（ただし実数の成分のみ）
        for(N=0;N<max;N++) {
          tmp[N]=f[N] = 0.f;
          fscanf(function,"%lf   %lf\n",&tmp[N], &f[N]);
          //fscanf(function,"%lf   %lf   %lf\n",&tmp[N], &f[N], &dtmp[N]);
        }

        N=10000;
        if(tmp[0] != 0.f){
            dt = tmp[0];
        }else{
            dt = tmp[1];
        }
        printf("%d\n",N);
        printf("dt   =  %e\n",dt);
//実数部分と虚数部分に分けてフーリエ変換
        for(n=0;n<N;n++) {
                ReF=ImF=0.0;
                om = 1.f/N/dt;
                for(k=0;k<N;k++) {
                        ReF+=f[k]*cos(2*pi*k*n/N);
                        ImF+=f[k]*sin(2*pi*k*n/N);
                }
                fprintf(fourier,"%f %e %e %e\n",n*om,ReF,ImF,ReF*ReF+ImF*ImF);             //最後にスペクトルの強さを追加
        }

        fclose(function);
        fclose(fourier);
        
        return 0;
}
