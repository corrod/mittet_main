#!/usr/bin/ruby

#####################################################
# Ruby script for making vertical section in FDTD.  #
# Written by Naoto Imamrua (2012/08/12)             #
#####################################################

dir_base = File.split("data2")

# Set the coordinate of plane
yplane = 250.0
zplane = 250.0

600.times{|i| #iterative start
char_in1 ="time_plane_%05d.dat"%[i*10].to_s
char_in2 ="time_slice_%05d.dat"%[i*10].to_s
char_out="time_ex_%05d.dat"%[i*10].to_s

linenum = -1
# Set the file information
fooE1 = File.open(File.join(dir_base,char_in1),"w")
fooE2 = File.open(File.join(dir_base,char_in2),"w")

# Starting extract

File.open(File.join(dir_base,char_out)){|f|
  f.each_line do |line|
    lnx = line.slice(0,12).to_f
    lny = line.slice(16,27).to_f
    lnz = line.slice(30,42).to_f

    if lnx == yplane then
      fooE1.puts line
    end
    if lny == yplane \
    && lnz == zplane then
      fooE2.puts line
    end
  end
}
fooE1.close
fooE2.close

}  # iterative end
