# Copyright (2012)2 (03-2014)3 Intel Corporation All Rights Reserved. 
# The source code contained or described herein and all documents related 
# to the source code ("Material") are owned by Intel Corporation or its suppliers 
# or licensors. Title to the Material remains with Intel Corporation or its 
# suppliers and licensors. The Material contains trade secrets and proprietary 
# and confidential information of Intel or its suppliers and licensors. The 
# Material is protected by worldwide copyright and trade secret laws and treaty 
# provisions. No part of the Material may be used, copied, reproduced, modified, 
# published, uploaded, posted, transmitted, distributed, or disclosed in any way 
# without Intel’s prior express written permission.
#
# No license under any patent, copyright, trade secret or other intellectual property 
# right is granted to or conferred upon you by disclosure or delivery of the 
# Materials, either expressly, by implication, inducement, estoppel or otherwise. 
# Any license under such intellectual property rights must be express 
# and approved by Intel in writing.




##################################################
##################################################
#########                              ###########
######### Master Makefile for ISO-3DFD ###########
#########                              ###########
##################################################
######### cedric.andreolli@intel.com   ###########
######### philippe.thierry@intel.com   ###########
##################################################
##################################################

include config/make.config

notarget:
	@echo ""
	@echo "####################################################################################"
	@echo "                             ISO-3DFD MASTER MAKEFILE"
	@echo "####################################################################################"
	@echo "For information about compilation and execution, refer to the file README.txt. This"
	@echo "file will provide information about each versions of the code."
	@echo ""
	@echo "This Makefie can't be run independently. A configuration file must be accessible at"
	@echo "the location config/Make.config. Feel free to edit it if you want." 
	@echo ""
	@echo ""
	@echo "To build the last version, call 'make last' which will call the Makefile with the following options :"
	@echo "make build model=$(CONFIG_MODEL) simd=$(simd) version=$(CONFIG_VERSION)"
	@echo ""
	@echo ""
	@echo "The default direcoty is $(CONFIG_ROOT_DIR)"
	@echo ""
	@echo ""
	@echo ""
	@echo "The default behavior can be re-defined in the file Makefile.config."
	@echo ""
	@echo "To get a decription of the option, run 'make help'."
	@echo ""
	@echo ""
	@echo ""

build:
	@echo ""
	@echo ""
	@echo "make last is first calling make clean"
	make clean
	@echo "Now that make clean has been called, calling make build model=\"$(CONFIG_MODEL)\" simd=\"$(simd)\" version=\"$(CONFIG_VERSION)\""
	cd $(CONFIG_S_DIR);make all model="$(CONFIG_MODEL)" simd="$(simd)" version="$(CONFIG_VERSION)" 
	cd $(CONFIG_ROOT_DIR)

last:
	make build

clean:
	@echo ""
	@echo ""
	@echo "Calling make clean model=\"$(CONFIG_MODEL)\" simd=\"$(simd)\" version=\"$(CONFIG_VERSION)\""
	cd $(CONFIG_S_DIR);make clean;$(CONFIG_RM) $(CONFIG_EXEC);cd $(CONFIG_ROOT_DIR)


cleanall:
	@echo ""
	@echo ""
	@echo "Cleaning all the project..."
	$(CONFIG_RM) $(CONFIG_BIN_DIR)/*

