
#
# Nasty little awk script to insert a CVS id comment line
# before the second executable statement
#

# $Id: fortranid.awk 19707 2010-10-29 17:59:36Z d3y133 $

BEGIN {
   FIRST = 0;
   DONE  = 0;
}

DONE == 1  {print; next;}

/^[ \t][ \t][ \t][ \t][ \t][ \t]/ {
		if (FIRST) {
			printf("C$Id: fortranid.awk 19707 2010-10-29 17:59:36Z d3y133 $\n");
			DONE = 1;
		} else {
			FIRST = 1;
		}
		print;
		next;
	}

		{print;}
