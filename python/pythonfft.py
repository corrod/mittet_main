#-*- coding:utf-8 -*-
import scipy as sp
import pylab as pl

#サンプルレートは1000とする。
sample_rate=1000.00

#０から０．６の間に0.001間隔で値を入れる。
t=sp.r_[0:0.6:1/sample_rate]

#サンプルの数を確認
N=len(t)

#サンプルとなる波を作成
s=sp.sin(2*sp.pi*50*t) + sp.sin(2*sp.pi*70*t+sp.pi/4)

#サンプル波を確認
pl.plot(t,s)
pl.show()
#サンプル波をFFT
S=sp.fft(s)

#周波数ビンを作成
f=sample_rate*sp.r_[0:(N/2)]/N

#周波数ビンの数を確認
n=len(f)

#FFTの結果を可視化
pl.plot(f,abs(S[0:n])/N)
pl.show()
