
DEV_CC_DIR ?= /usr/linux-?1om-4.7/bin
DEV_DIR ?= /opt/intel/mic
COI_DIR ?= $(DEV_DIR)/coi
COI_INC ?= -I $(COI_DIR)/include

HOST_FLAGS ?= $(COI_INC) -lcoi_host -Wl,--enable-new-dtags
DEBUG_FLAGS = -g -O0 -D_DEBUG
RELEASE_FLAGS = -DNDEBUG -O3 
DEBUG_DIR = debug
RELEASE_DIR = release

# for GCC compiler
HOST_CC ?= g++
DEV_CC ?= $(DEV_CC_DIR)/x86_64-?1om-linux-g++
DEV_FLAGS ?= $(COI_INC) -lcoi_device -rdynamic -Wl,--enable-new-dtags

# for ICC compiler
#HOST_CC ?= icc 
#DEV_CC ?= icc
#DEV_FLAGS ?= $(COI_INC) -mmic -lcoi_device -rdynamic -Wl,--enable-new-dtags

all: $(DEBUG_DIR)/sink_mic $(RELEASE_DIR)/sink_mic \
     $(DEBUG_DIR)/source_stream $(RELEASE_DIR)/source_stream \
     $(DEBUG_DIR)/source_normal $(RELEASE_DIR)/source_normal \
     $(DEBUG_DIR)/source_pinned $(RELEASE_DIR)/source_pinned 

$(DEBUG_DIR)/source_stream: source_stream.cpp
	mkdir -p $(DEBUG_DIR)
	$(HOST_CC) -L$(COI_DIR)/host-linux-debug/lib -Wl,-rpath=$(COI_DIR)/host-linux-debug/lib $(HOST_FLAGS) $(DEBUG_FLAGS) -o $@ $<

$(RELEASE_DIR)/source_stream: source_stream.cpp
	mkdir -p $(RELEASE_DIR)
	$(HOST_CC) -L$(COI_DIR)/host-linux-release/lib -Wl,-rpath=$(COI_DIR)/host-linux-release/lib $(HOST_FLAGS) $(RELEASE_FLAGS) -o $@ $<

$(DEBUG_DIR)/source_normal: source_normal.cpp
	mkdir -p $(DEBUG_DIR)
	$(HOST_CC) -L$(COI_DIR)/host-linux-debug/lib -Wl,-rpath=$(COI_DIR)/host-linux-debug/lib $(HOST_FLAGS) $(DEBUG_FLAGS) -o $@ $<

$(RELEASE_DIR)/source_normal: source_normal.cpp
	mkdir -p $(RELEASE_DIR)
	$(HOST_CC) -L$(COI_DIR)/host-linux-release/lib -Wl,-rpath=$(COI_DIR)/host-linux-release/lib $(HOST_FLAGS) $(RELEASE_FLAGS) -o $@ $<

$(DEBUG_DIR)/source_pinned: source_pinned.cpp
	mkdir -p $(DEBUG_DIR)
	$(HOST_CC) -L$(COI_DIR)/host-linux-debug/lib -Wl,-rpath=$(COI_DIR)/host-linux-debug/lib $(HOST_FLAGS) $(DEBUG_FLAGS) -o $@ $<

$(RELEASE_DIR)/source_pinned: source_pinned.cpp
	mkdir -p $(RELEASE_DIR)
	$(HOST_CC) -L$(COI_DIR)/host-linux-release/lib -Wl,-rpath=$(COI_DIR)/host-linux-release/lib $(HOST_FLAGS) $(RELEASE_FLAGS) -o $@ $<

$(DEBUG_DIR)/sink_mic: sink.cpp
	mkdir -p $(DEBUG_DIR)
	$(DEV_CC) -L$(COI_DIR)/device-linux-debug/lib $(DEV_FLAGS) $(DEBUG_FLAGS) -o $@ $<

$(RELEASE_DIR)/sink_mic: sink.cpp
	mkdir -p $(RELEASE_DIR)
	$(DEV_CC) -L$(COI_DIR)/device-linux-release/lib $(DEV_FLAGS) $(RELEASE_FLAGS) -o $@ $<

clean:
	rm -rf debug
	rm -rf release
