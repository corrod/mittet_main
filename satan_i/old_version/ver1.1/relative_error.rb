#!/usr/bin/ruby

##############################################################
# Ruby script for calculating relative error in time domain. #
# Written by Naoto Imamura (2012/10/04)                      #
##############################################################

dir_base = File.split("data")

char1 ="rex_000_000.dat"
char2 ="anal2.dat"
char3 ="relative_error.dat"

linenum = 0
refmax  = 0
array1  =[]
# Set the file information
fooE1 = File.open(File.join(dir_base,char1),"r")
fooE2 = File.open(File.join(dir_base,char2),"r")
fooE3 = File.open(File.join(dir_base,char3),"w")

# Starting READ
File.open(File.join(dir_base,char2)){|f|
  f.each_line do |line|
    ltime = line.slice(0,12).to_f
    lnum  = line.slice(14,29).to_f
    if lnum.abs > refmax.abs
      refmax = lnum
    end
    array1 << lnum
    linenum=linenum+1
end
}
fooE2.close

linenum=0
File.open(File.join(dir_base,char1)){|f|
  f.each_line do |line|
    ltime = line.slice(0,10).to_f
    lrec  = line.slice(12,27).to_f
    rel   = (lrec - array1[linenum]).abs / refmax
    fooE3.printf("%f    %f\n", ltime, rel);
    linenum=linenum+1
end
}

fooE1.close
fooE3.close
