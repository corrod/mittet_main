cd ./

cordinate = dlmread('./pycordinate.dat');


data = dlmread('./combine_pattern2_12diff.dat');

cordinate1 = transpose(cordinate);

data(:,95) = [];

for i = 1:4351
    data(i, :) = data(i, :) * i^1;
end

z = 0:1.25846858562e-6:0.00547559681601
wigb(data(:,:),1,cordinate1,z); axis([19 104 0 0.00547559681601]);
xlabel('Distance [mm]');
ylabel('Travel Time [s]');