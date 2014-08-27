インダクタンスL
http://keisan.casio.jp/exec/user/1225887110
http://homepage2.nifty.com/kaoru~i/coil.htm

L = 長岡係数A * 真空の透磁率μ0 * 比透磁率μr * コイルの半径a * 巻き数N * コイルの長さb
    A * (4PI * 10d-7 * μ0) * (PI * a**2) * N**2 / b


k=1/sqrt(1+(L/(2*R))^2)
k2=k^2
長岡係数 Kn=(4/(3*pi*sqrt(1-k2))*(((1-k2)/k2)*ellipticK(k)-((1-2*k2)/k2)*ellipticE(k)-k)



磁場の強度Hと磁束密度B

真空中では
B = μ0H

媒質中では
B = μ0H + M
M = χH

特に線形媒質中では
B = μH


誘導起電力V
V = - NdΦ/dt
Φ = ∫BdS

自己誘導V
V = - LdI/dt

相互誘導
