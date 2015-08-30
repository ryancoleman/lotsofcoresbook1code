Copyright (c) The University of Tennessee.  All rights reserved.

$COPYRIGHT$

Additional copyrights may follow

$HEADER$

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer listed in this license
in the documentation and/or other materials provided with the distribution.

- Neither the name of the copyright holders nor the names of its contributors
  may be used to endorse or promote products derived from this software
without specific prior written permission.

The copyright holders provide no reassurances that the source code provided
does not infringe any patent, copyright, or any other intellectual property
rights of third parties.  The copyright holders disclaim any liability to any
recipient for claims brought against recipient by any third party for
infringement of that parties intellectual property rights.


$NOTES$
To compile the test driver, please execute:
  make clean ; make pdlltdriver3.exe ; make pdlltdriver4.exe
The executable pdlltdriver3.exe represents the version running on MIC device. The executable pdlltdriver4.exe represents the version only using CPU.
Users should be able to follow the example in "run_llt_mic.sh" to run the test code.
Users need to make sure the MKL library can provide the scalapack
functionalities, which are needed for compilation and linking of the code. If the user does not have a MKL with scalapack support. The user can install the scalapack separately. For the user convinience, a tarball of scalapack-2.0.0.tgz is also included in the package.
