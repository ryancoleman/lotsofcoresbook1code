Archive file(s) in this directory dated in 2014 will correspond to 
version(s) used in the book.

This code assumes that MPI is working an that /opt/intel is mounted 
via a network file-system like NFS.

To run with multiple cards on a standard Centos 7 distribution:
	1) Disable the firewall (or open ports for MPI in the firewall)

	   a) See http://research.colfaxinternational.com/file.axd?file=2013%2F10%2FColfax_Heterogeneous_Clustering_Xeon_Phi.pdf

	2) Disable SELinux or set it to permissive

	3) Enable IPv4 forwarding in the virtual network that connects the coprocessors

	  a) http://www.centos.org/docs/5/html/Virtual_Server_Administration/s1-lvs-forwarding-VSA.html

To run with one card (and potentially avoid the issues with SELinux and 
port forarding), use *_oneCard.zip



