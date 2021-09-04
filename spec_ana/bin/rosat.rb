#!/usr/bin/env ruby

ra1=ARGV[0].to_s
ra2=ARGV[1].to_s
ra3=ARGV[2].to_s
dec1=ARGV[3].to_s
dec2=ARGV[4].to_s
dec3=ARGV[5].to_s
ra=ra1+" "+ra2+" "+ra3
dec=dec1+" "+dec2+" "+dec3
radius=2.0
inner=1.0

_ra = ra.gsub(' ','+').gsub('.','%2E')
_dec = dec.gsub('+','%2B').gsub(' ','+').gsub('.','%2E')

head = 'http://heasarc.gsfc.nasa.gov/cgi-bin/Tools/xraybg/'

require "open-uri"

xml = ""
open(head+"xraybg.pl?Entry=#{_ra}%2C+#{_dec}&CoordSys=J2000&NR=GRB%2FSIMBAD%2BSesame%2FNED&radius=#{radius}&region=annulus&inner_radius=#{inner}&spectrum=create&scaling=hist") do |f|
  xml = f.read
end


a=xml.index('xraybgfits.pl')
b=xml.index('">FITS spectrum')

fitsfile=xml[a...b]

piname="rass_orig.pi"
str='wget '+head+fitsfile+' -O '+piname
system(str)

a=xml.index('/cgi-bin/W3Browse/startfv.pl')
b=xml.index('">ROSAT PSPCc response matrix')

fitsfile=xml[a...b]

respname="pspcc.rsp"
str='wget '+'http://heasarc.gsfc.nasa.gov'+fitsfile+' -O '+respname
system(str)
