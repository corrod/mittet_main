cd C:\WORK\‹Æ–±\Îì“‡”d–dH‹Æ\’T¸ŽŽŒ±\ŒŸ“¢‚R\waveforms\”gŒ`•`‰æ

cordinate = csvread('coodinate.txt');

data = csvread('receiver-h-200.dat');

cordinate1 = transpose(cordinate);

data(:,269) = [];

for i = 1:3276
    data(i, :) = data(i, :) * i^1;
end

z = 0:0.001:3.275;
wigb(data(:,:),1,cordinate1,z); axis([0 3360 0 3.2]);
xlabel('Distance [mm]');
ylabel('Travel Time [s]');
Title('Source -> 001');
