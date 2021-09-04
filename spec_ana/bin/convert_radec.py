import sys
import math
from astropy import units as u
from astropy.coordinates import SkyCoord

if len(sys.argv) != 3 : 
   print "INPUT ERROR in python";
   sys.exit();


ra0=float(sys.argv[1])
dec0=float(sys.argv[2])

c=SkyCoord(ra=ra0*u.degree, dec=dec0*u.degree, frame='icrs')

ra="%d %d %lf" % (c.ra.hms[0],c.ra.hms[1],c.ra.hms[2])
dec="%d %d %lf" % (c.dec.dms[0],abs(c.dec.dms[1]),abs(c.dec.dms[2]))
print ra,dec


