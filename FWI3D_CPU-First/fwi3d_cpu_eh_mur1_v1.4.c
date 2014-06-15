#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "fwi3d_cpu_alloex1_v1.2.h"
#include "fwi3d_cpu_eh_mur1_v1.2.h"
//////////////////////////////////////////
void init_eh_mur_3d()
{
    int i;

//  for xy plane
    for(i=0;i<4*ix*iy;i++){
        ex1_xy[i] = 0.f;
        ey1_xy[i] = 0.f;
        ex2_xy[i] = 0.f;
        ey2_xy[i] = 0.f;
    }
//  for yz plane
    for(i=0;i<4*iy*iz;i++){
        ey1_yz[i] = 0.f;
        ez1_yz[i] = 0.f;
        ey2_yz[i] = 0.f;
        ez2_yz[i] = 0.f;
    }
//  for zx plane
    for(i=0;i<4*ix*iz;i++){
        ez1_zx[i] = 0.f;
        ex1_zx[i] = 0.f;
        ez2_zx[i] = 0.f;
        ex2_zx[i] = 0.f;
    }
}
//////////////////////////////////////////
void e_mur1(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    e_mur1_xy1(EX, EY, EZ, HX, HY, HZ);
    e_mur1_xy2(EX, EY, EZ, HX, HY, HZ);

    e_mur1_yz1(EX, EY, EZ, HX, HY, HZ);
    e_mur1_yz2(EX, EY, EZ, HX, HY, HZ);

    e_mur1_zx1(EX, EY, EZ, HX, HY, HZ);
    e_mur1_zx2(EX, EY, EZ, HX, HY, HZ);
}
//////////////////////////////////////////
void e_mur1_xy1(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double ptz = (cmax*dt-dz)/(cmax*dt+dz);

    for(i=0;i<ix-1;i++){
        for(j=1;j<iy-1;j++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            EX[ijk] = ex1_xy[ijk + ix*iy] + ptz*(EX[ijk + ix*iy] - ex1_xy[ijk]);
        }
    }

    for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            EY[ijk] = ey1_xy[ijk + ix*iy] + ptz*(EY[ijk + ix*iy] - ey1_xy[ijk]);
        }
    }

    for(i=0;i<ix-1;i++){
        for(j=1;j<iy-1;j++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            ex1_xy[ijk + ix*iy] = EX[ijk + ix*iy];
            ex1_xy[ijk        ] = EX[ijk        ];
        }
    }

    for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            ey1_xy[ijk + ix*iy] = EY[ijk + ix*iy];
            ey1_xy[ijk        ] = EY[ijk        ];
        }
    }
}
//////////////////////////////////////////
void e_mur1_xy2(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double ptz = (cmax*dt-dz)/(cmax*dt+dz);

    for(i=0;i<ix-1;i++){
        for(j=1;j<iy-1;j++){
            k=iz-1;
            ijk = k*ix*iy + j*ix + i;
            EX[ijk] = ex1_xy[j*ix +i + 2*ix*iy] + ptz*(EX[ijk - ix*iy] - ex1_xy[j*ix+i + 3*ix*iy]);
        }
    }

    for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
            k=iz-1;
            ijk = k*ix*iy + j*ix + i;
            EY[ijk] = ey1_xy[j*ix +i + 2*ix*iy] + ptz*(EY[ijk - ix*iy] - ey1_xy[j*ix +i + 3*ix*iy]);
        }
    }

    for(i=0;i<ix-1;i++){
        for(j=1;j<iy-1;j++){
            k=iz-1;
            ijk = k*ix*iy + j*ix + i;
            ex1_xy[j*ix +i + 2*ix*iy] = EX[ijk - ix*iy];
            ex1_xy[j*ix +i + 3*ix*iy] = EX[ijk        ];
        }
    }

    for(j=0;j<iy-1;j++){
        for(i=1;i<ix-1;i++){
            k=iz-1;
            ijk = k*ix*iy + j*ix + i;
            ey1_xy[j*ix +i + 2*ix*iy] = EY[ijk - ix*iy];
            ey1_xy[j*ix +i + 3*ix*iy] = EY[ijk        ];
        }
    }
}
//////////////////////////////////////////
void e_mur1_yz1(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double ptx = (cmax*dt-dx)/(cmax*dt+dx);

    for(j=0;j<iy-1;j++){
        for(k=1;k<iz-1;k++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            EY[ijk] = ey1_yz[iy*iz + j*iz +k] + ptx*(EY[ijk + 1] - ey1_yz[j*iz +k]);
        }
    }

    for(k=0;k<iz-1;k++){
        for(j=1;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            EZ[ijk] = ez1_yz[iy*iz + j*iz +k] + ptx*(EZ[ijk + 1] - ez1_yz[j*iz +k]);
        }
    }

    for(j=0;j<iy-1;j++){
        for(k=1;k<iz-1;k++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            ey1_yz[iy*iz + j*iz +k] = EY[ijk + 1];
            ey1_yz[j*iz + k       ] = EY[ijk    ];
        }
    }

    for(k=0;k<iz-1;k++){
        for(j=1;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            ez1_yz[iy*iz + j*iz +k] = EZ[ijk + 1];
            ez1_yz[j*iz + k       ] = EZ[ijk    ];
        }
    }
}
//////////////////////////////////////////
void e_mur1_yz2(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double ptx = (cmax*dt-dx)/(cmax*dt+dx);

    for(j=0;j<iy-1;j++){
        for(k=1;k<iz-1;k++){
            i=ix-1;
            ijk = k*ix*iy + j*ix + i;
            EY[ijk] = ey1_yz[2*iy*iz + j*iz +k] + ptx*(EY[ijk - 1] - ey1_yz[3*iy*iz + j*iz +k]);
        }
    }

    for(k=0;k<iz-1;k++){
        for(j=1;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            EZ[ijk] = ez1_yz[iy*iz + j*iz +k] + ptx*(EZ[ijk + 1] - ez1_yz[3*iy*iz + j*iz +k]);
        }
    }

    for(j=0;j<iy-1;j++){
        for(k=1;k<iz-1;k++){
            i=ix-1;
            ijk = k*ix*iy + j*ix + i;
            ey1_yz[2*iy*iz + j*iz +k] = EY[ijk - 1];
            ey1_yz[3*iy*iz + j*iz +k] = EY[ijk    ];
        }
    }

    for(k=0;k<iz-1;k++){
        for(j=1;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            ez1_yz[2*iy*iz + j*iz +k] = EZ[ijk - 1];
            ez1_yz[3*iy*iz + j*iz +k] = EZ[ijk    ];
        }
    }
}
//////////////////////////////////////////
void e_mur1_zx1(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double pty = (cmax*dt-dx)/(cmax*dt+dx);

    for(i=0;i<ix-1;i++){
        for(k=1;k<iz-1;k++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            EX[ijk] = ex1_zx[ix*iz + i*iz +k] + pty*(EX[ijk + ix] - ex1_zx[i*iz +k]);
        }
    }

    for(k=0;k<iz-1;k++){
        for(i=1;i<ix-1;i++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            EZ[ijk] = ez1_zx[ix*iz + i*iz +k] + pty*(EZ[ijk + ix] - ez1_zx[i*iz +k]);
        }
    }

    for(i=0;i<ix-1;i++){
        for(k=1;k<iz-1;k++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            ex1_zx[ix*iz + i*iz +k] = EX[ijk + ix];
            ex1_zx[i*iz + k       ] = EX[ijk     ];
        }
    }

    for(k=0;k<iz-1;k++){
        for(i=1;i<ix-1;i++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            ez1_zx[ix*iz + i*iz +k] = EZ[ijk + ix];
            ez1_zx[i*iz + k       ] = EZ[ijk     ];
        }
    }
}
//////////////////////////////////////////
void e_mur1_zx2(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double pty = (cmax*dt-dx)/(cmax*dt+dx);

    for(i=0;i<ix-1;i++){
        for(k=1;k<iz-1;k++){
            j=iy-1;
            ijk = k*ix*iy + j*ix + i;
            EX[ijk] = ex1_zx[2*ix*iz + i*iz +k] + pty*(EX[ijk - ix] - ex1_zx[3*ix*iz + i*iz +k]);
        }
    }

    for(k=0;k<iz-1;k++){
        for(i=1;i<ix-1;i++){
            j=iy-1;
            ijk = k*ix*iy + j*ix + i;
            EZ[ijk] = ez1_zx[2*ix*iz + i*iz +k] + pty*(EZ[ijk - ix] - ez1_zx[3*ix*iz + i*iz +k]);
        }
    }

    for(i=0;i<ix-1;i++){
        for(k=1;k<iz-1;k++){
            j=iy-1;
            ijk = k*ix*iy + j*ix + i;
            ex1_zx[2*ix*iz + i*iz +k] = EX[ijk - ix];
            ex1_zx[3*ix+iz + i*iz +k] = EX[ijk     ];
        }
    }

    for(k=0;k<iz-1;k++){
        for(i=1;i<ix-1;i++){
            j=iy-1;
            ijk = k*ix*iy + j*ix + i;
            ez1_zx[2*ix*iz + i*iz +k] = EZ[ijk - ix];
            ez1_zx[3*ix*iz + i*iz +k] = EZ[ijk     ];
        }
    }
}
//////////////////////////////////////////
//      Mur の2次吸収境界条件 (p.65)
//
//           ----------------
//       z  /|     ④      /|
//        /  |    (6)    /  |
//       ----------------   |
//       |(2)|y         |(5)|
//       |   -----------|----
//       |  /    (1)    |  /
//       |/   ③        |/
//       ---------------- x
//
//      ○: 表の面, (*): 隠れた面
//      1: z=0, x-y plane
//      2: x=0, y-z plane
//      3: y=0, x-z plane
//      4: z=z, x-y plane
//      5: x=x, y-z plane
//      6: y=y, x-z plane
//
//////////////////////////////////////////
void e_mur2(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    e_mur2_xy(EX, EY, EZ, HX, HY, HZ);
    e_mur2_yz(EX, EY, EZ, HX, HY, HZ);
    e_mur2_zx(EX, EY, EZ, HX, HY, HZ);
}
//////////////////////////////////////////
//  x-y平面に対して ( z軸方向に伝播 )   //
//////////////////////////////////////////
void e_mur2_xy(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double czd   = (cmax*dt-dz)/(cmax*dt+dz);
    double czu   = (cmax*dt-dz)/(cmax*dt+dz);
    double czz   = (2.f*dz) / (cmax*dt+dz);
    double czfxd = (dz*pow((cmax*dt),2.f)) / (2.f*(dx*dx)*(cmax*dt+dz));
    double czfyd = (dz*pow((cmax*dt),2.f)) / (2.f*(dy*dy)*(cmax*dt+dz));

//  ************** 壁 k=1,nzに対して *************
//  -------------- Exに対して -------------
//  1次吸収境界条件
    for(j=1;j<iy-1;j++){
        i=0;
        k=0;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exz1[j*ix+i+ ix*iy] + czd*(EX[ijk + ix*iy] - exz1[j*ix+i]);
        k=iz-1;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exz1[j*ix+i + 2*ix*iy] + czu*(EX[ijk - ix*iy] - exz1[j*ix+i + 3*ix*iy]);

        i=ix-2;
        k=0;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exz1[j*ix+i+ ix*iy] + czd*(EX[ijk + ix*iy] - exz1[j*ix+i]);
        k=iz-1;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exz1[j*ix+i + 2*ix*iy] + czu*(EX[ijk - ix*iy] - exz1[j*ix+i + 3*ix*iy]);
    }

    for(i=0;i<ix-2;i++){
         j=1;
         k=0;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exz1[j*ix+i + ix*iy] + czd*(EX[ijk + ix*iy] - exz1[j*ix+i]);
         k=iz-1;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exz1[j*ix+i + 2*ix*iy] + czu*(EX[ijk - ix*iy] - exz1[j*ix+i + 3*ix*iy]);

         j=iy-2;
         k=0;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exz1[j*ix+i + ix*iy] + czd*(EX[ijk + ix*iy] - exz1[j*ix+i]);
         k=iz-1;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exz1[j*ix+i + 2*ix*iy] + czu*(EX[ijk - ix*iy] - exz1[j*ix+i + 3*ix*iy]);
    }
//  2次吸収境界条件
    for(i=1;i<ix-2;i++){
      for(j=2;j<iy-2;j++){
        k=0;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = -exz2[j*ix+i + ix*iy] + czd*(EX[ijk+ix*iy] + exz2[j*ix+i]) \
             + czz*(exz1[j*ix+i]   + exz1[j*ix+i + ix*iy]) \
             + czfxd* (exz1[j*ix+i+1] - 2.f*exz1[j*ix+i] + exz1[j*ix+i-1]  \
                  + exz1[j*ix+i+1 +ix*iy] - 2.f*exz1[j*ix+i+ix*iy] + exz1[j*ix+i-1 +ix*iy]) \
             + czfyd* (exz1[(j+1)*ix+i] - 2.f*exz1[j*ix+i] + exz1[(j-1)*ix+i]  \
                  + exz1[(j+1)*ix+i +ix*iy] - 2.f*exz1[j*ix+i+ix*iy] + exz1[(j-1)*ix+i +ix*iy]);
        k=iz-1;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = -exz2[j*ix+i + 2*ix*iy] + czd*(EX[ijk-ix*iy] + exz2[j*ix+i+3*ix*iy]) \
             + czz*(exz1[j*ix+i+3*ix*iy]   + exz1[j*ix+i + 2*ix*iy]) \
             + czfxd* (exz1[j*ix+i+1+3*ix*iy] - 2.f*exz1[j*ix+i+3*ix*iy] + exz1[j*ix+i-1+3*ix*iy]  \
                  + exz1[j*ix+i+1 +2*ix*iy] - 2.f*exz1[j*ix+i+2*ix*iy] + exz1[j*ix+i-1 +2*ix*iy]) \
             + czfyd* (exz1[(j+1)*ix+i+3*ix*iy] - 2.f*exz1[j*ix+i+3*ix*iy] + exz1[(j-1)*ix+i+3*ix*iy]  \
                  + exz1[(j+1)*ix+i +2*ix*iy] - 2.f*exz1[j*ix+i+2*ix*iy] + exz1[(j-1)*ix+i +2*ix*iy]);
      }
    }

    for(i=0;i<ix-1;i++){
        for(j=1;j<iy-1;j++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            exz2[j*ix+i         ] = exz1[j*ix+i         ];
            exz2[j*ix+i +1*ix*iy] = exz1[j*ix+i +1*ix*iy];
            exz2[j*ix+i +2*ix*iy] = exz1[j*ix+i +2*ix*iy];
            exz2[j*ix+i +3*ix*iy] = exz1[j*ix+i +3*ix*iy];
        }
    }
    for(i=0;i<ix-1;i++){
        for(j=1;j<iy-1;j++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            exz1[j*ix+i         ] = EX[j*ix+i         ];
            exz1[j*ix+i +1*ix*iy] = EX[j*ix+i +1*ix*iy];
            exz1[j*ix+i +2*ix*iy] = EX[j*ix+i +(iz-2)*ix*iy];
            exz1[j*ix+i +3*ix*iy] = EX[j*ix+i +(iz-1)*ix*iy];
        }
    }
//  -------------- Eyに対して -------------
//  1次吸収境界条件
    for(j=0;j<iy-1;j++){
        i=1;
        k=0;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyz1[j*ix+i+ ix*iy] + czd*(EY[ijk + ix*iy] - eyz1[j*ix+i]);
        k=iz-1;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyz1[j*ix+i + 2*ix*iy] + czu*(EY[ijk - ix*iy] - eyz1[j*ix+i + 3*ix*iy]);

        i=ix-2;
        k=0;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyz1[j*ix+i+ ix*iy] + czd*(EY[ijk + ix*iy] - eyz1[j*ix+i]);
        k=iz-1;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyz1[j*ix+i + 2*ix*iy] + czu*(EY[ijk - ix*iy] - eyz1[j*ix+i + 3*ix*iy]);
    }

    for(i=2;i<ix-2;i++){
         j=0;
         k=0;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyz1[j*ix+i + ix*iy] + czd*(EY[ijk + ix*iy] - eyz1[j*ix+i]);
         k=iz-1;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyz1[j*ix+i + 2*ix*iy] + czu*(EY[ijk - ix*iy] - eyz1[j*ix+i + 3*ix*iy]);

         j=iy-2;
         k=0;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyz1[j*ix+i + ix*iy] + czd*(EY[ijk + ix*iy] - eyz1[j*ix+i]);
         k=iz-1;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyz1[j*ix+i + 2*ix*iy] + czu*(EY[ijk - ix*iy] - eyz1[j*ix+i + 3*ix*iy]);
    }
//  2次吸収境界条件
    for(i=2;i<ix-2;i++){
      for(j=1;j<iy-2;j++){
        k=0;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = -eyz2[j*ix+i + ix*iy] + czd*(EY[ijk+ix*iy] + eyz2[j*ix+i]) \
             + czz*(eyz1[j*ix+i]   + eyz1[j*ix+i + ix*iy]) \
             + czfxd* (eyz1[j*ix+i+1] - 2.f*eyz1[j*ix+i] + eyz1[j*ix+i-1]  \
                  + eyz1[j*ix+i+1 +ix*iy] - 2.f*eyz1[j*ix+i+ix*iy] + eyz1[j*ix+i-1 +ix*iy]) \
             + czfyd* (eyz1[(j+1)*ix+i] - 2.f*eyz1[j*ix+i] + eyz1[(j-1)*ix+i]  \
                  + eyz1[(j+1)*ix+i +ix*iy] - 2.f*eyz1[j*ix+i+ix*iy] + eyz1[(j-1)*ix+i +ix*iy]);
        k=iz-1;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = -eyz2[j*ix+i + 2*ix*iy] + czd*(EY[ijk-ix*iy] + eyz2[j*ix+i+3*ix*iy]) \
             + czz*(eyz1[j*ix+i+3*ix*iy]   + eyz1[j*ix+i + 2*ix*iy]) \
             + czfxd* (eyz1[j*ix+i+1 +3*ix*iy] - 2.f*eyz1[j*ix+i +3*ix*iy] + eyz1[j*ix+i-1 +3*ix*iy]  \
                     + eyz1[j*ix+i+1 +2*ix*iy] - 2.f*eyz1[j*ix+i +2*ix*iy] + eyz1[j*ix+i-1 +2*ix*iy]) \
             + czfyd* (eyz1[(j+1)*ix+i +3*ix*iy] - 2.f*eyz1[j*ix+i +3*ix*iy] + eyz1[(j-1)*ix+i +3*ix*iy]  \
                     + eyz1[(j+1)*ix+i +2*ix*iy] - 2.f*eyz1[j*ix+i +2*ix*iy] + eyz1[(j-1)*ix+i +2*ix*iy]);
      }
    }

    for(i=1;i<ix-1;i++){
        for(j=0;j<iy-1;j++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            eyz2[j*ix+i         ] = eyz1[j*ix+i         ];
            eyz2[j*ix+i +1*ix*iy] = eyz1[j*ix+i +1*ix*iy];
            eyz2[j*ix+i +2*ix*iy] = eyz1[j*ix+i +2*ix*iy];
            eyz2[j*ix+i +3*ix*iy] = eyz1[j*ix+i +3*ix*iy];
        }
    }
    for(i=1;i<ix-1;i++){
        for(j=0;j<iy-1;j++){
            k=0;
            ijk = k*ix*iy + j*ix + i;
            eyz1[j*ix+i         ] = EY[j*ix+i         ];
            eyz1[j*ix+i +1*ix*iy] = EY[j*ix+i +1*ix*iy];
            eyz1[j*ix+i +2*ix*iy] = EY[j*ix+i +(iz-2)*ix*iy];
            eyz1[j*ix+i +3*ix*iy] = EY[j*ix+i +(iz-1)*ix*iy];
        }
    }
}
//////////////////////////////////////////
//  y-z平面に対して ( x軸方向に伝播 )   //
//////////////////////////////////////////
void e_mur2_yz(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double cxd   = (cmax*dt-dx)/(cmax*dt+dx);
    double cxu   = (cmax*dt-dx)/(cmax*dt+dx);
    double cxx   = (2.f*dx) / (cmax*dt+dx);
    double cxfyd = (dx*pow((cmax*dt),2.f)) / (2.f*(dy*dy)*(cmax*dt+dx));
    double cxfzd = (dx*pow((cmax*dt),2.f)) / (2.f*(dz*dz)*(cmax*dt+dx));

//  ************** 壁 i=1,nxに対して *************
//  -------------- Eyに対して -------------
//  1次吸収境界条件
    for(k=1;k<iz-1;k++){
        j=0;
        i=0;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyx1[k*iy+j+ iy*iz] + cxd*(EY[ijk + 1] - eyx1[k*iy+j]);
        i=ix-1;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyx1[k*iy+j + 2*iy*iz] + cxu*(EY[ijk - 1] - eyx1[k*iy+j + 3*iy*iz]);

        j=iy-2;
        i=0;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyx1[k*iy+j+ iy*iz] + cxd*(EY[ijk + 1] - eyx1[k*iy+j]);
        i=ix-1;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = eyx1[k*iy+j + 2*iy*iz] + cxu*(EY[ijk - 1] - eyx1[k*iy+j + 3*iy*iz]);
    }

    for(j=1;j<iy-2;j++){
         k=1;
         i=0;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyx1[k*iy+j + iy*iz] + cxd*(EY[ijk + 1] - eyx1[k*iy+j]);
         i=ix-1;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyx1[k*iy+j + 2*iy*iz] + cxu*(EY[ijk - 1] - eyx1[k*iy+j + 3*iy*iz]);

         k=iz-2;
         i=0;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyx1[k*iy+j + iy*iz] + cxd*(EY[ijk + 1] - eyx1[k*iy+j]);
         i=ix-1;
         ijk = k*ix*iy + j*ix + i;
         EY[ijk] = eyx1[k*iy+j + 2*iy*iz] + cxu*(EY[ijk - 1] - eyx1[k*iy+j + 3*iy*iz]);
    }
//  2次吸収境界条件
    for(k=2;k<iz-2;k++){
      for(j=1;j<iy-2;j++){
        i=0;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = -eyx2[k*ix+j + ix*iy] + cxd*(EY[ijk+1] + eyx2[k*iy+j]) \
             + cxx*(eyx1[k*iy+j]   + eyx1[k*iy+j + iy*iz]) \
             + cxfyd* (eyx1[k*iy+j+1       ]   - 2.f*eyx1[k*iy+j      ] + eyx1[k*iy+j-1       ]  \
                     + eyx1[k*iy+j+1 +iy*iz]   - 2.f*eyx1[k*iy+j+iy*iz] + eyx1[k*iy+j-1 +iy*iz]) \
             + cxfzd* (eyx1[(k+1)*iy+j     ]   - 2.f*eyx1[k*iy+j      ] + eyx1[(k-1)*iy+j     ]  \
                     + eyx1[(k+1)*iy+j +iy*iz] - 2.f*eyx1[k*iy+j+iy*iz] + eyx1[(k-1)*iy+j +iy*iz]);
        i=ix-1;
        ijk = k*ix*iy + j*ix + i;
        EY[ijk] = -eyx2[k*iy+j + 2*iy*iz] + cxd*(EY[ijk-1] + eyx2[k*iy+j+3*iy*iz]) \
             + cxx*(eyx1[k*iy+j+3*iy*iz]   + eyx1[k*iy+j + 2*iy*iz]) \
             + cxfyd* (eyx1[k*iy+j+1 +3*iy*iz]   - 2.f*eyx1[k*iy+j+3*iy*iz] + eyx1[k*iy+j-1 +3*iy*iz]  \
                     + eyx1[k*iy+j+1 +2*iy*iz]   - 2.f*eyx1[k*iy+j+2*iy*iz] + eyx1[k*iy+j-1 +2*iy*iz]) \
             + cxfzd* (eyx1[(k+1)*iy+j +3*iy*iz] - 2.f*eyx1[k*iy+j+3*iy*iz] + eyx1[(k-1)*iy+j +3*iy*iz]  \
                     + eyx1[(k+1)*iy+j +2*iy*iz] - 2.f*eyx1[k*iy+j+2*iy*iz] + eyx1[(k-1)*iy+j +2*iy*iz]);
      }
    }

    for(k=1;k<iz-1;k++){
        for(j=0;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            eyx2[k*iy+j         ] = eyx1[k*iy+j         ];
            eyx2[k*iy+j +1*iy*iz] = eyx1[k*iy+j +1*iy*iz];
            eyx2[k*iy+j +2*iy*iz] = eyx1[k*iy+j +2*iy*iz];
            eyx2[k*iy+j +3*iy*iz] = eyx1[k*iy+j +3*iy*iz];
        }
    }
    for(k=1;k<iz-1;k++){
        for(j=0;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            eyx1[k*iy+j         ] = EY[k*iy*ix + j*ix      ];
            eyx1[k*iy+j +1*iy*iz] = EY[k*iy*ix + j*ix +1   ];
            eyx1[k*iy+j +2*iy*iz] = EY[k*iy*ix + j*ix +ix-2];
            eyx1[k*iy+j +3*iy*iz] = EY[k*iy*ix + j*ix +ix-1];
        }
    }
//  -------------- Ezに対して -------------
//  1次吸収境界条件
    for(k=0;k<iz-1;k++){
        j=1;
        i=0;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezx1[k*iy+j+ iy*iz] + cxd*(EZ[ijk + 1] - ezx1[k*iy+j]);
        i=ix-1;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezx1[k*iy+j + 2*iy*iz] + cxu*(EZ[ijk - 1] - ezx1[k*iy+j + 3*iy*iz]);

        j=iy-2;
        i=0;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezx1[k*iy+j+    iy*iz] + cxd*(EZ[ijk + 1] - ezx1[k*iy+j]);
        i=ix-1;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezx1[k*iy+j + 2*iy*iz] + cxu*(EZ[ijk - 1] - ezx1[k*iy+j + 3*iy*iz]);
    }

    for(j=2;j<iy-2;j++){
         k=0;
         i=0;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezx1[k*iy+j +   iy*iz] + cxd*(EZ[ijk + 1] - ezx1[k*iy+j]);
         i=ix-1;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezx1[k*iy+j + 2*iy*iz] + cxu*(EZ[ijk - 1] - ezx1[k*iy+j + 3*iy*iz]);

         k=iz-2;
         i=0;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezx1[k*iy+j +   iy*iz] + cxd*(EZ[ijk + 1] - ezx1[k*iy+j]);
         i=ix-1;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezx1[k*iy+j + 2*iy*iz] + cxu*(EZ[ijk - 1] - ezx1[k*iy+j + 3*iy*iz]);
    }
//  2次吸収境界条件
    for(k=1;k<iz-2;k++){
      for(j=2;j<iy-2;j++){
        i=0;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = -ezx2[k*iy+j + iy*iz] + cxd*(EZ[ijk+1] + ezx2[k*iy+j]) \
             + cxx*(ezx1[k*iy+j]   + ezx1[k*iy+j + iy*iz]) \
             + cxfyd* (ezx1[k*iy+j+1         ] - 2.f*ezx1[k*iy+j      ] + ezx1[k*iy+j-1]  \
                     + ezx1[k*iy+j+1 +iy*iz  ] - 2.f*ezx1[k*iy+j+iy*iz] + ezx1[k*iy+j-1 +iy*iz]) \
             + cxfzd* (ezx1[(k+1)*iy+j       ] - 2.f*ezx1[k*iy+j      ] + ezx1[(k-1)*iy+j]  \
                     + ezx1[(k+1)*iy+j +iy*iz] - 2.f*ezx1[k*iy+j+iy*iz] + ezx1[(k-1)*iy+j +iy*iz]);
        i=ix-1;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = -ezx2[k*iy+j + 2*iy*iz] + cxd*(EZ[ijk-1] + ezx2[k*iy+j+3*iy*iz]) \
             + cxx*(ezx1[k*iy+j+3*iy*iz]   + ezx1[k*iy+j + 2*iy*iz]) \
             + cxfyd* (ezx1[k*iy+j+1   +3*iy*iz] - 2.f*ezx1[k*iy+j +3*iy*iz] + ezx1[k*iy+j-1   +3*iy*iz]  \
                     + ezx1[k*iy+j+1   +2*iy*iz] - 2.f*ezx1[k*iy+j +2*iy*iz] + ezx1[k*iy+j-1   +2*iy*iz]) \
             + cxfzd* (ezx1[(k+1)*iy+j +3*iy*iz] - 2.f*ezx1[k*iy+j +3*iy*iz] + ezx1[(k-1)*iy+j +3*iy*iz]  \
                     + ezx1[(k+1)*iy+j +2*iy*iz] - 2.f*ezx1[k*iy+j +2*iy*iz] + ezx1[(k-1)*iy+j +2*iy*iz]);
      }
    }

    for(k=0;k<iz-1;k++){
        for(j=1;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            ezx2[k*iy+j         ] = ezx1[k*iy+j         ];
            ezx2[k*iy+j +1*iy*iz] = ezx1[k*iy+j +1*iy*iz];
            ezx2[k*iy+j +2*iy*iz] = ezx1[k*iy+j +2*iy*iz];
            ezx2[k*iy+j +3*iy*iz] = ezx1[k*iy+j +3*iy*iz];
        }
    }
    for(k=0;k<iz-1;k++){
        for(j=1;j<iy-1;j++){
            i=0;
            ijk = k*ix*iy + j*ix + i;
            ezx1[k*iy+j         ] = EZ[k*iy*ix + j*ix        ];
            ezx1[k*iy+j +1*iy*iz] = EZ[k*iy*ix + j*ix +1     ];
            ezx1[k*iy+j +2*iy*iz] = EZ[k*iy*ix + j*ix +(ix-2)];
            ezx1[k*iy+j +3*iy*iz] = EZ[k*iy*ix + j*ix +(ix-1)];
        }
    }
}
//////////////////////////////////////////
//  x-z平面に対して ( y軸方向に伝播 )   //
//////////////////////////////////////////
void e_mur2_zx(double *EX, double *EY, double *EZ, double *HX, double *HY, double *HZ)
{
    int i,j,k,ijk;
    double cyd   = (cmax*dt-dy)/(cmax*dt+dy);
    double cyu   = (cmax*dt-dy)/(cmax*dt+dy);
    double cyy   = (2.f*dy) / (cmax*dt+dx);
    double cyfxd = (dy*pow((cmax*dt),2.f)) / (2.f*(dx*dx)*(cmax*dt+dy));
    double cyfzd = (dy*pow((cmax*dt),2.f)) / (2.f*(dz*dz)*(cmax*dt+dy));

//  ************** 壁 j=1,nyに対して *************
//  -------------- Exに対して -------------
//  1次吸収境界条件
    for(k=1;k<iz-1;k++){
        i=0;
        j=0;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exy1[k*ix+i   + ix*iz] + cyd*(EX[ijk + ix] - exy1[k*ix+i]);
        j=iy-1;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exy1[k*ix+i + 2*ix*iz] + cyu*(EX[ijk - ix] - exy1[k*ix+i + 3*ix*iz]);

        i=ix-2;
        j=0;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exy1[k*ix+i   + ix*iz] + cyd*(EX[ijk + ix] - exy1[k*ix+i]);
        j=iy-1;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = exy1[k*ix+i + 2*ix*iz] + cyu*(EX[ijk - ix] - exy1[k*ix+i + 3*ix*iz]);
    }

    for(i=1;i<ix-2;i++){
         k=1;
         j=0;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exy1[k*ix+i +   ix*iz] + cyd*(EX[ijk + ix] - exy1[k*ix+i]);
         j=iy-1;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exy1[k*ix+i + 2*ix*iz] + cyu*(EX[ijk - ix] - exy1[k*ix+i + 3*ix*iz]);

         k=iz-2;
         j=0;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exy1[k*ix+i +   ix*iz] + cyd*(EX[ijk + ix] - exy1[k*ix+i]);
         j=iy-1;
         ijk = k*ix*iy + j*ix + i;
         EX[ijk] = exy1[k*ix+i + 2*ix*iz] + cyu*(EX[ijk - ix] - exy1[k*ix+i + 3*ix*iz]);
    }
//  2次吸収境界条件
    for(k=2;k<iz-2;k++){
      for(i=1;i<ix-2;i++){
        j=0;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = -exy2[k*ix+i + ix*iz] + cyd*(EX[ijk+ix] + exy2[k*ix+i]) \
             + cyy*(exy1[k*ix+i]   + exy1[k*ix+i + ix*iz]) \
             + cyfxd* (exy1[k*ix+i+1         ] - 2.f*exy1[k*ix+i      ] + exy1[k*ix+i-1       ]  \
                     + exy1[k*ix+i+1 +  ix*iz] - 2.f*exy1[k*ix+i+ix*iz] + exy1[k*ix+i-1 +ix*iz]) \
             + cyfzd* (exy1[(k+1)*ix+i       ] - 2.f*exy1[k*ix+i      ] + exy1[(k-1)*ix+i     ]  \
                     + exy1[(k+1)*ix+i +ix*iz] - 2.f*exy1[k*ix+i+ix*iz] + exy1[(k-1)*ix+i +ix*iz]);
        j=iy-1;
        ijk = k*ix*iy + j*ix + i;
        EX[ijk] = -exy2[k*ix+i + 2*ix*iz] + cyd*(EX[ijk-ix] + exy2[k*ix+i+3*ix*iz]) \
             + cyy*(exy1[k*ix+i+3*ix*iz]   + exy1[k*ix+i + 2*ix*iz]) \
             + cyfxd* (exy1[k*ix+i+1   +3*ix*iz] - 2.f*exy1[k*ix+i+3*ix*iz] + exy1[k  *ix+i-1 +3*ix*iz]  \
                     + exy1[k*ix+i+1   +2*ix*iz] - 2.f*exy1[k*ix+i+2*ix*iz] + exy1[k  *ix+i-1 +2*ix*iz]) \
             + cyfzd* (exy1[(k+1)*ix+i +3*ix*iz] - 2.f*exy1[k*ix+i+3*ix*iz] + exy1[(k-1)*ix+i +3*ix*iz]  \
                     + exy1[(k+1)*ix+i +2*ix*iz] - 2.f*exy1[k*ix+i+2*ix*iz] + exy1[(k-1)*ix+i +2*ix*iz]);
      }
    }

    for(k=1;k<iz-1;k++){
        for(i=0;i<ix-1;i++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            exy2[k*ix+i         ] = exy1[k*ix+i         ];
            exy2[k*ix+i +1*ix*iz] = exy1[k*ix+i +1*ix*iz];
            exy2[k*ix+i +2*ix*iz] = exy1[k*ix+i +2*ix*iz];
            exy2[k*ix+i +3*ix*iz] = exy1[k*ix+i +3*ix*iz];
        }
    }
    for(k=1;k<iz-1;k++){
        for(i=0;i<ix-1;i++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            exy1[k*ix+i         ] = EX[k*ix*iy + i           ];
            exy1[k*ix+i +1*ix*iz] = EX[k*ix*iy + i +1*ix     ];
            exy1[k*ix+i +2*ix*iz] = EX[k*ix*iy + i +(iy-2)*ix];
            exy1[k*ix+i +3*ix*iz] = EX[k*ix*iy + i +(iy-1)*ix];
        }
    }
//  -------------- Ezに対して -------------
//  1次吸収境界条件
    for(k=0;k<iz-1;k++){
        i=1;
        j=0;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezy1[k*ix+i+    ix*iz] + cyd*(EZ[ijk + ix] - ezy1[k*ix+i]);
        j=iy-1;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezy1[k*ix+i + 2*ix*iz] + cyu*(EZ[ijk - ix] - ezy1[k*ix+i + 3*ix*iz]);

        i=ix-2;
        j=0;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezy1[k*ix+i+    ix*iz] + cyd*(EZ[ijk + ix] - ezy1[k*ix+i]);
        j=iy-1;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = ezy1[k*ix+i + 2*ix*iz] + cyu*(EZ[ijk - ix] - ezy1[k*ix+i + 3*ix*iz]);
    }

    for(i=2;i<ix-2;i++){
         k=0;
         j=0;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezy1[k*ix+i +   ix*iz] + cyd*(EZ[ijk + ix] - ezy1[k*ix+i]);
         j=iy-1;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezy1[k*ix+i + 2*ix*iz] + cyu*(EZ[ijk - ix] - ezy1[k*ix+i + 3*ix*iz]);

         k=iz-2;
         j=0;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezy1[k*ix+i +   ix*iz] + cyd*(EZ[ijk + ix] - ezy1[k*ix+i]);
         j=iy-1;
         ijk = k*ix*iy + j*ix + i;
         EZ[ijk] = ezy1[k*ix+i + 2*ix*iz] + cyu*(EZ[ijk - ix] - ezy1[k*ix+i + 3*ix*iz]);
    }
//  2次吸収境界条件
    for(k=1;k<iz-2;k++){
      for(i=2;i<ix-2;i++){
        j=0;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = -ezy2[k*ix+i + ix*iz] + cyd*(EZ[ijk+ix] + ezy2[k*ix+i]) \
             + cyy*(ezy1[k*ix+i]   + ezy1[k*ix+i + ix*iz]) \
             + cyfxd* (ezy1[k*ix+i+1         ] - 2.f*ezy1[k*ix+i      ] + ezy1[k*ix+j-1]  \
                     + ezy1[k*ix+i+1 +ix*iz  ] - 2.f*ezy1[k*ix+i+ix*iz] + ezy1[k*ix+j-1 +ix*iz]) \
             + cyfzd* (ezy1[(k+1)*ix+i       ] - 2.f*ezy1[k*ix+i      ] + ezy1[(k-1)*ix+i]  \
                     + ezy1[(k+1)*ix+i +ix*iz] - 2.f*ezy1[k*ix+i+ix*iz] + ezy1[(k-1)*ix+i +ix*iz]);
        j=iy-1;
        ijk = k*ix*iy + j*ix + i;
        EZ[ijk] = -ezy2[k*ix+i + 2*ix*iz] + cyd*(EZ[ijk-ix] + ezy2[k*ix+i+3*ix*iz]) \
             + cyy*(ezy1[k*ix+i+ 3*ix*iz] + ezy1[k*ix+i + 2*ix*iz]) \
             + cyfxd* (ezy1[k*ix+i+1   +3*ix*iz] - 2.f*ezy1[k*ix+i +3*ix*iz] + ezy1[k*ix+i-1   +3*ix*iz]  \
                     + ezy1[k*ix+i+1   +2*ix*iz] - 2.f*ezy1[k*ix+i +2*ix*iz] + ezy1[k*ix+i-1   +2*ix*iz]) \
             + cyfzd* (ezy1[(k+1)*ix+i +3*ix*iz] - 2.f*ezy1[k*ix+i +3*ix*iz] + ezy1[(k-1)*ix+i +3*ix*iz]  \
                     + ezy1[(k+1)*ix+i +2*ix*iz] - 2.f*ezy1[k*ix+i +2*ix*iz] + ezy1[(k-1)*ix+i +2*ix*iz]);
      }
    }

    for(k=0;k<iz-1;k++){
        for(i=1;i<ix-1;i++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            ezy2[k*ix+i         ] = ezy1[k*ix+i         ];
            ezy2[k*ix+i +1*ix*iz] = ezy1[k*ix+i +1*ix*iz];
            ezy2[k*ix+i +2*ix*iz] = ezy1[k*ix+i +2*ix*iz];
            ezy2[k*ix+i +3*ix*iz] = ezy1[k*ix+i +3*ix*iz];
        }
    }
    for(k=0;k<iz-1;k++){
        for(i=1;i<ix-1;i++){
            j=0;
            ijk = k*ix*iy + j*ix + i;
            ezy1[k*ix+i         ] = EZ[k*ix*iy+i         ];
            ezy1[k*ix+i +1*ix*iz] = EZ[k*ix*iy+i +1*ix];
            ezy1[k*ix+i +2*ix*iz] = EZ[k*ix*iy+i +(iy-2)*ix];
            ezy1[k*ix+i +3*ix*iz] = EZ[k*ix*iy+i +(iy-1)*ix];
        }
    }
}
