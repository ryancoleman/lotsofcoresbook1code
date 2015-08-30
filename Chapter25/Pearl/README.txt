This file contains the information about building and running of Asian
Options Pricing demo with dynamic load balancing ("Boss-workers" model)
from Intel Xeon Phi GEMS book:
Chapter 25: Heterogeneous MPI application optimization with ITAC

1. Code compilation:
$> make

2. Code run:
In order to run the code MPI libraries and scripts should be shared between
the host and coprocessors. One way to do this, is to configure NFS share of
/opt/intel/impi folder or even /opt/intel folder on coprocessors. But this
operation needs to be done with the root privileges:

$> sudo micctrl --addnfs=/opt/intel --dir=/opt/intel --server=localhost
$> sudo service mpss restart

Alternatively, required files can be copied directly onto coprocessors, but we
will let the reader to set up the environment on her/his own.

For multiple devices firewall and SELinux should be turned off, and IP
packets forwarding enabled. Run the following commands with the root
privileges:

$> su # root access
$> echo 0 > /selinux/enforce
$> service iptables stop
$> /sbin/sysctl -w net.ipv4.ip_forward=1

or for the last step modify line from /etc/sysctl.conf to:
net.ipv4.ip_forward = 1

After the preparations, just run the following command:

$> make run

or 

$> ./run-options.sh cpu
$> ./run-options.sh mic
$> ./run-options.sh both

This will produce output similar to the following one:

#         Worker   Share Performance  Effic.
       c001-n003   16.0%    8.00e+06   93.5%
  c001-n003-mic0   21.0%    1.04e+07   94.1%
  c001-n003-mic1   21.0%    1.03e+07   95.1%
  c001-n003-mic2   21.0%    1.03e+07   95.0%
  c001-n003-mic3   21.0%    1.05e+07   93.8%
# Calculation   3 of  10 took 2.245 seconds
# Net performance: 4.67e+07 paths/second

Efficiency is measured by using ratio between times used on useful
calculations, and time with included MPI communication with the BOSS rank.
