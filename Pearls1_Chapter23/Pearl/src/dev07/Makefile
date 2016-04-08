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


## Main config file
include ../../config/make.config


## Objects and sources
CSRCS	= $(wildcard *.cc)
COBJS	= $(CSRCS:.cc=.o)


default:
	@echo "Calling make build model=$(CONFIG_MODEL) simd=$(CONFIG_SIMD) version=$(CONFIG_VERSION) root=$(CONFIG_ROOT_DIR)"
	@echo "CSRCS=$(CSRCS)"
	@echo "COBJS=$(COBJS)"

all:$(COBJS)
	$(CONFIG_CC) $(CONFIG_LDFLAGS) $(COBJS) $(CONFIG_LIBS) -o $(CONFIG_EXEC)  

clean:
	$(CONFIG_RM) $(COBJS) $(CONFIG_EXEC) *.o



