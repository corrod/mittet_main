#!/usr/bin/ruby

#####################################################
# Ruby script for making vertical section in FDTD.  #
# Written by Naoto Imamrua (2012/10/08)             #
#####################################################

dir_base = File.split("data")

# Set the coordinate of plane

char_in1 ="conductivity.dat"
char_out="section_cond.dat"

linenum = -1
# Set the file information
fooE1 = File.open(File.join(dir_base,char_out),"w")

# Starting extract
ax=[]
ay=[]
az=[]
File.open(File.join(dir_base,char_in1)){|f|
  f.each_line do |line|
    ax << line.slice(0,12).to_f
    ay << line.slice(16,27).to_f
    az << line.slice(30,42).to_f
  end
}
max_x = ax.max
max_y = ay.max
max_z = az.max
p max_x
p max_y
p max_z
xplane = max_x/2
yplane = max_y/2
zplane = max_z/2
fooE1.close

linenum = -1
# Set the file information
fooE1 = File.open(File.join(dir_base,char_out),"w")

# Starting extract

File.open(File.join(dir_base,char_in1)){|f|
  f.each_line do |line|
    lnx = line.slice(0,12).to_f
    lny = line.slice(16,27).to_f
    lnz = line.slice(30,42).to_f

    if lny == yplane then
      fooE1.puts line
    end
  end
}
fooE1.close
