# attempt to determine the type of a file from its name

# $Id: filetype.awk 19707 2010-10-29 17:59:36Z d3y133 $

/\.fh$/			{print "Fortran-header"; next;}

/\.F$/ || /\.f$/	{print "Fortran"; next;}

/\.h$/			{print "C-header"; next;}

/\.c$/			{print "C"; next;}

/[Mm]akefile$/		{print "Makefile"; next;}

			{print "Unknown"; next;}
