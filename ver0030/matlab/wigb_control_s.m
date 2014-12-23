cd ./

cordinate = dlmread('./out_residual/pycordinate.dat');

%data = csvread('./out_residual/test.dat'); 

data = dlmread('./combine_residual_abso.dat');

cordinate1 = transpose(cordinate);

data(:,85) = [];

for i = 1:4351
    data(i, :) = data(i, :) * i^1;
end

%z = 0:1.25846858562e-6:0.00547559681601;
z = 0:1.25846858562e-6:0.00547559681601
wigb(data(:,:),1,cordinate1,z); axis([0 85 0 0.00547559681601]);
xlabel('Distance [mm]');
ylabel('Travel Time [s]');
%Title('residual_differentialmat');
