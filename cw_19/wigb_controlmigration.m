cd ./ver0030r_a

cordinate = dlmread('./pycordinate.dat');


data = dlmread('./migration.dat');

cordinate1 = transpose(cordinate);

data(:,95) = [];

for i = 1:19
    data(i, :) = data(i, :) * i^1;
end

% z = 0 : 1.25846858562e-6 : 2.391090312678e-05
z = 0 : 1.25846858562e-6 : 2.51693717124e-05
% wigb(data(:,:),1,cordinate1,z); axis([19 114 0 2.391090312678e-05]);
wigb(data(:,:),1,cordinate1,z); axis([19 114 0 2.51693717124e-05]);
xlabel('Distance [mm]');
ylabel('Travel Time [s]');