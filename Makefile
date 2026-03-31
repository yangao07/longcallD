.DEFAULT_GOAL := all

GCC_CHECK := $(shell gcc --version | head -n 1 | grep -i "clang")

# Check if the OS is macOS or linux
UNAME_S := $(shell uname -s)

EXTRA_CFLAGS ?=
EXTRA_LDFLAGS ?=
EXTRA_LIBS ?=
PROFILER_LIBDIR ?=
GMON ?= gmon.out

ifneq ($(gprof),)
ifneq ($(pg),)
$(error use either gprof=1 or pg=1, not both)
endif
pg := $(gprof)
endif

ifeq ($(GCC_CHECK),) # gcc
	CXXFLAGS = -std=c++11
else # clang
	CXXFLAGS = -std=c++11 -stdlib=libc++
endif

# add -fno-tree-vectorize to avoid certain vectorization errors in O3 optimization
# right now, we are using -O3 for the best performance, and no vectorization errors were found
EXTRA_FLAGS = -Wall -Wno-misleading-indentation -Wno-unused-function #-Wno-unused-variable -Wno-alloc-size-larger-than

# Define the version number
VERSION=0.0.10
# Get the Git commit hash
GIT_COMMIT := $(shell git rev-parse --short HEAD 2> /dev/null)
ifneq ($(GIT_COMMIT),)
	LONGCALLD_VERSION = $(VERSION)-$(GIT_COMMIT)
else
	LONGCALLD_VERSION = $(VERSION)
endif

HTSLIB_DIR  = ./htslib
HTSLIB      = $(HTSLIB_DIR)/libhts.a

EDLIB_DIR     = ./edlib
EDLIB_INC_DIR = ./edlib/include
EDLIB         = $(EDLIB_DIR)/src/edlib.o

ABPOA_DIR   = ./abPOA
ABPOA_LIB   = ./abPOA/lib/libabpoa.a
ABPOA_INC_DIR = $(ABPOA_DIR)/include

WFA2_DIR    = ./WFA2-lib
WFA2_LIB    = $(WFA2_DIR)/lib/libwfa.a

LIB_PATH    =
LINK_LDFLAGS =
LINK_LIBS    = $(EXTRA_LIBS)
LIB         = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) $(LINK_LDFLAGS) -lm -lz -lpthread -llzma -lbz2 -lcurl $(LINK_LIBS)
INCLUDE     = -I $(HTSLIB_DIR) -I $(EDLIB_INC_DIR) -I $(ABPOA_INC_DIR) -I $(WFA2_DIR)

# Try linking against libdeflate
ifeq ($(shell echo "int main() {return 0;}" | ${CC} -x c - $(LIB_PATH) -ldeflate >/dev/null 2>&1 && echo "yes"),yes)
	LIB += -ldeflate
endif

ifeq ($(UNAME_S),Linux) # Linux
	LIB += -lcrypto
	ifneq ($(portable),)
		LIB += -static-libgcc -static-libstdc++
		ifneq ($(opt_lib),)
			LIB = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) $(LINK_LDFLAGS) -static-libgcc -static-libstdc++ -L${opt_lib} -lm -lz -lpthread -llzma -lbz2 -lcurl -lssl -lcrypto -lssh2 -ldeflate -lzstd $(LINK_LIBS)
		else
			ifneq ($(OPT_LIB),)
				LIB = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) $(LINK_LDFLAGS) -static-libgcc -static-libstdc++ -L${OPT_LIB} -lm -lz -lpthread -llzma -lbz2 -lcurl -lssl -lcrypto -lssh2 -ldeflate -lzstd $(LINK_LIBS)
			endif
		endif
	endif
endif

# for debug
ifneq ($(debug),)
	EXTRA_FLAGS  += -D __DEBUG__
endif

ABPOA_GDB_LIB = ./abPOA/lib/libabpoa_sse41.a
ABPOA_NOR_LIB = ./abPOA/lib/libabpoa.a
WFA_GDB_LIB   = ./WFA2-lib/lib/libwfa_gdb.a
WFA_NOR_LIB   = ./WFA2-lib/lib/libwfa.a

BUILD_MODE_NAME = release
OPT_FLAGS = -O3
THIRD_PARTY_OPT_FLAGS = -O3
PROFILE_LDFLAGS =
PROFILE_LIBS =

ifneq ($(gdb),)
ifneq ($(profile),)
$(error gdb and profile cannot be enabled together)
endif
endif

ifneq ($(profile),)
ifneq ($(pg),)
$(error pg and profile cannot be enabled together)
endif
endif

# for gdb
ifneq ($(gdb),)
	BUILD_MODE_NAME = gdb
	OPT_FLAGS = -O0 -g
	THIRD_PARTY_OPT_FLAGS = -O0 -g
	ABPOA_LIB = $(ABPOA_GDB_LIB)
	WFA2_LIB  = $(WFA_GDB_LIB)
else
	ABPOA_LIB = $(ABPOA_NOR_LIB)
	WFA2_LIB  = $(WFA_NOR_LIB)
endif

ifneq ($(profile),)
	BUILD_MODE_NAME = profile
	OPT_FLAGS = -O3 -g3 -fno-omit-frame-pointer -DNDEBUG
	THIRD_PARTY_OPT_FLAGS = -O3 -g3 -fno-omit-frame-pointer -DNDEBUG
	ifeq ($(UNAME_S),Darwin)
		ifeq ($(PROFILER_LIBDIR),)
			PROFILER_PREFIX := $(shell brew --prefix gperftools 2>/dev/null)
			ifneq ($(PROFILER_PREFIX),)
				PROFILER_LIBDIR := $(PROFILER_PREFIX)/lib
			endif
		endif
	endif
	ifneq ($(PROFILER_LIBDIR),)
		PROFILE_LDFLAGS += -L$(PROFILER_LIBDIR) -Wl,-rpath,$(PROFILER_LIBDIR)
	endif
	ifeq ($(UNAME_S),Linux)
		PROFILE_LIBS += -Wl,--no-as-needed -lprofiler -Wl,--as-needed
	else
		PROFILE_LIBS += -lprofiler
	endif
endif

HTSLIB_BUILD_CFLAGS = -Wall -fvisibility=hidden $(THIRD_PARTY_OPT_FLAGS) $(EXTRA_CFLAGS)
HTSLIB_BUILD_LDFLAGS = $(EXTRA_LDFLAGS) $(PROFILE_LDFLAGS)
WFA2_BUILD_CC_FLAGS = -Wall -fPIE $(THIRD_PARTY_OPT_FLAGS) $(EXTRA_CFLAGS)
CFLAGS = $(OPT_FLAGS) $(EXTRA_FLAGS) $(EXTRA_CFLAGS) -DLONGCALLD_VERSION=\"$(LONGCALLD_VERSION)\"

# for gprof
ifneq ($(pg),)
	BUILD_MODE_NAME = gprof
	OPT_FLAGS = -O3 -g3 -fno-omit-frame-pointer
	THIRD_PARTY_OPT_FLAGS = -O3 -g3 -fno-omit-frame-pointer
	PG_FLAG   = -pg
	CFLAGS   += -pg
	HTSLIB_BUILD_CFLAGS += -pg
	WFA2_BUILD_CC_FLAGS += -pg
endif

LINK_LDFLAGS += $(LIB_PATH) $(EXTRA_LDFLAGS) $(PROFILE_LDFLAGS)
LINK_LIBS += $(PROFILE_LIBS)

ifneq ($(PREFIX),)
	OUT_PRE_DIR = $(PREFIX)
else
	OUT_PRE_DIR = .
endif

BIN_DIR = $(OUT_PRE_DIR)/bin
INC_DIR = ./include
SRC_DIR = ./src

HTS_ALL = hts_all
EDLIB_ALL = edlib_all
ABPOA_ALL = abpoa_all
WFA2_ALL = wfa2_all
CSOURCE    = $(wildcard ${SRC_DIR}/*.c)
CPPSOURCE  = $(wildcard $(SRC_DIR)/*.cpp)
CPPSOURCE += $(EDLIB_DIR)/src/edlib.cpp

OBJS    = $(CSOURCE:.c=.o) $(CPPSOURCE:.cpp=.o)
BIN     = $(BIN_DIR)/longcallD
ifneq ($(gdb),)
	BIN = $(BIN_DIR)/gdb_longcallD
endif

BUILD_CONFIG_FILE = .build-config.mk
ABPOA_BUILD_ARGS =
ifneq ($(pg),)
	ABPOA_BUILD_ARGS += pg=1
endif

UNIT_TEST_SRCS = test/unit/main.c test/unit/test_bam_utils.c test/unit/test_collect_var.c test/unity/unity.c
UNIT_TEST_OBJS = $(UNIT_TEST_SRCS:.c=.o)
UNIT_TEST_INCLUDE = $(INCLUDE) -I $(SRC_DIR) -I test/unity
UNIT_APP_OBJS = $(filter-out $(SRC_DIR)/main.o,$(OBJS))
UNITY_DEFS = -DUNITY_INCLUDE_DOUBLE -DUNITY_DOUBLE_PRECISION=1e-6

.PHONY: FORCE
FORCE:

$(BUILD_CONFIG_FILE): FORCE
	@tmp=$@.tmp; \
	{ \
		printf 'BUILD_MODE=%s\n' '$(BUILD_MODE_NAME)'; \
		printf 'CC=%s\n' '$(CC)'; \
		printf 'CXX=%s\n' '$(CXX)'; \
		printf 'CFLAGS=%s\n' '$(CFLAGS)'; \
		printf 'HTSLIB_BUILD_CFLAGS=%s\n' '$(HTSLIB_BUILD_CFLAGS)'; \
		printf 'HTSLIB_BUILD_LDFLAGS=%s\n' '$(HTSLIB_BUILD_LDFLAGS)'; \
		printf 'WFA2_BUILD_CC_FLAGS=%s\n' '$(WFA2_BUILD_CC_FLAGS)'; \
		printf 'ABPOA_BUILD_ARGS=%s\n' '$(ABPOA_BUILD_ARGS)'; \
		printf 'LINK_LDFLAGS=%s\n' '$(LINK_LDFLAGS)'; \
		printf 'LINK_LIBS=%s\n' '$(LINK_LIBS)'; \
		printf 'BIN=%s\n' '$(BIN)'; \
		printf 'ABPOA_LIB=%s\n' '$(ABPOA_LIB)'; \
		printf 'WFA2_LIB=%s\n' '$(WFA2_LIB)'; \
	} > $$tmp; \
	if [ ! -f $@ ] || ! cmp -s $$tmp $@; then mv $$tmp $@; else rm -f $$tmp; fi

$(SRC_DIR)/%.o: $(SRC_DIR)/%.c $(BUILD_CONFIG_FILE)
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all: $(HTS_ALL) $(EDLIB) $(ABPOA_LIB) $(WFA2_LIB) $(BIN)

$(HTS_ALL): $(HTSLIB)

$(HTSLIB): $(HTSLIB_DIR)/configure.ac $(BUILD_CONFIG_FILE)
	cd $(HTSLIB_DIR) && autoreconf -i && ./configure CC="$(CC)" CFLAGS="$(HTSLIB_BUILD_CFLAGS)" LDFLAGS="$(HTSLIB_BUILD_LDFLAGS)" && $(MAKE) clean && $(MAKE) CC="$(CC)" CFLAGS="$(HTSLIB_BUILD_CFLAGS)" LDFLAGS="$(HTSLIB_BUILD_LDFLAGS)" libhts.a

$(EDLIB): $(EDLIB_DIR)/src/edlib.cpp $(EDLIB_DIR)/include/edlib.h $(BUILD_CONFIG_FILE)
	$(CXX) $(CFLAGS) $(CXXFLAGS) -c $< $(INCLUDE) -o $@
$(EDLIB_ALL): $(EDLIB)

$(ABPOA_GDB_LIB): $(BUILD_CONFIG_FILE)
	cd $(ABPOA_DIR); $(MAKE) clean; $(MAKE) libabpoa gdb=1 sse41=1 OPT_FLAGS="$(THIRD_PARTY_OPT_FLAGS)" $(ABPOA_BUILD_ARGS)
$(ABPOA_NOR_LIB): $(BUILD_CONFIG_FILE)
	cd $(ABPOA_DIR); $(MAKE) clean; $(MAKE) libabpoa OPT_FLAGS="$(THIRD_PARTY_OPT_FLAGS)" $(ABPOA_BUILD_ARGS)

$(ABPOA_ALL): $(ABPOA_LIB)

$(WFA2_LIB): $(BUILD_CONFIG_FILE)
	cd $(WFA2_DIR); $(MAKE) clean; $(MAKE) setup lib_wfa CC="$(CC)" CXX="$(CXX)" CC_FLAGS="$(WFA2_BUILD_CC_FLAGS)" 
$(WFA2_ALL): $(WFA2_LIB)


$(BIN): $(OBJS) $(ABPOA_LIB) $(HTSLIB) $(WFA2_LIB) $(BUILD_CONFIG_FILE)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CXX) $(OBJS) -o $@ $(LIB) $(PG_FLAG)
# 	$(CC) $(OBJS) -o $@ $(LIB) $(PG_FLAG)

$(SRC_DIR)/align.o: $(SRC_DIR)/align.c $(SRC_DIR)/align.h $(SRC_DIR)/utils.h $(SRC_DIR)/bam_utils.h $(SRC_DIR)/seq.h $(BUILD_CONFIG_FILE)
	$(CC) -c -DUSE_SIMDE -DSIMDE_ENABLE_NATIVE_ALIASES $(CFLAGS) $< $(INCLUDE) -o $@

$(SRC_DIR)/assign_aln_hap.o: $(SRC_DIR)/assign_aln_hap.c $(SRC_DIR)/assign_aln_hap.h $(SRC_DIR)/utils.h $(SRC_DIR)/bam_utils.h
$(SRC_DIR)/bam_utils.o: $(SRC_DIR)/bam_utils.c $(SRC_DIR)/bam_utils.h $(SRC_DIR)/utils.h
$(SRC_DIR)/cgranges.o: $(SRC_DIR)/cgranges.c $(SRC_DIR)/cgranges.h $(SRC_DIR)/khash.h
$(SRC_DIR)/collect_var.o: $(SRC_DIR)/collect_var.c $(SRC_DIR)/collect_var.h $(SRC_DIR)/bam_utils.h
$(SRC_DIR)/kalloc.o: $(SRC_DIR)/kalloc.c $(SRC_DIR)/kalloc.h
$(SRC_DIR)/kmedoids.o : $(SRC_DIR)/kmedoids.c $(SRC_DIR)/kmedoids.h
$(SRC_DIR)/kthread.o: $(SRC_DIR)/kthread.c
$(SRC_DIR)/main.o: $(SRC_DIR)/main.c $(SRC_DIR)/call_var_main.h
$(SRC_DIR)/call_var_main.o: $(SRC_DIR)/bam_utils.c $(SRC_DIR)/call_var_main.c $(SRC_DIR)/call_var_main.h $(SRC_DIR)/main.h $(SRC_DIR)/utils.h $(SRC_DIR)/seq.h \
                            $(SRC_DIR)/collect_var.h
$(SRC_DIR)/seq.o: $(SRC_DIR)/seq.c $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h
$(SRC_DIR)/sdust.o: $(SRC_DIR)/sdust.c $(SRC_DIR)/sdust.h $(SRC_DIR)/kdq.h $(SRC_DIR)/kvec.h
$(SRC_DIR)/utils.o: $(SRC_DIR)/utils.c $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h $(SRC_DIR)/kseq.h
$(SRC_DIR)/vcf_utils.o: $(SRC_DIR)/vcf_utils.c $(SRC_DIR)/vcf_utils.h $(SRC_DIR)/utils.h

test/unit/%.o: test/unit/%.c $(BUILD_CONFIG_FILE)
	$(CC) -c $(CFLAGS) $(UNIT_TEST_INCLUDE) $(UNITY_DEFS) $< -o $@

test/unity/%.o: test/unity/%.c $(BUILD_CONFIG_FILE)
	$(CC) -c $(CFLAGS) -I test/unity $(UNITY_DEFS) $< -o $@

test/unit/unit_tests: $(UNIT_APP_OBJS) $(UNIT_TEST_OBJS) $(ABPOA_LIB) $(HTSLIB) $(WFA2_LIB) $(BUILD_CONFIG_FILE)
	$(CXX) $(UNIT_APP_OBJS) $(UNIT_TEST_OBJS) -o $@ $(LIB) $(PG_FLAG)

.PHONY: help test test_unit gprof_build gprof_report

help:
	@printf "\n"
	@printf "Build Modes:\n"
	@printf "  make bin/longcallD              Release build\n"
	@printf "  make gdb=1 bin/longcallD        Debug build\n"
	@printf "  make profile=1 bin/longcallD    gperftools/pprof build\n"
	@printf "  make pg=1 bin/longcallD         gprof build\n"
	@printf "  make gprof=1 bin/longcallD      Alias for pg=1\n"
	@printf "\n"
	@printf "Profiling Helpers:\n"
	@printf "  make gprof_build                Clean rebuild with pg=1\n"
	@printf "  make gprof_report               Run gprof on $$(pwd)/$(GMON)\n"
	@printf "\n"
	@printf "Other Targets:\n"
	@printf "  make test\n"
	@printf "  make test_unit\n"
	@printf "  make clean\n"
	@printf "  make clean_all\n"
	@printf "\n"

test: test_unit

test_unit: test/unit/unit_tests
	./test/unit/unit_tests

gprof_build:
	$(MAKE) clean_all
	$(MAKE) pg=1 bin/longcallD

gprof_report:
	@test -x $(BIN_DIR)/longcallD || { \
		printf "missing %s\n" "$(BIN_DIR)/longcallD"; \
		printf "run 'make gprof_build' or 'make gprof=1 bin/longcallD' first\n"; \
		exit 1; \
	}
	@test -f $(GMON) || { \
		printf "missing %s\n" "$(GMON)"; \
		printf "run the gprof-instrumented binary to produce gmon.out first\n"; \
		exit 1; \
	}
	gprof $(BIN_DIR)/longcallD $(GMON)

.PHONY: all hts_all abpoa_all wfa2_all clean clean_all clean_hts clean_abpoa clean_wfa2

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)
clean_all:
	rm -f $(SRC_DIR)/*.o $(BIN) $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) $(UNIT_TEST_OBJS) test/unit/unit_tests $(BUILD_CONFIG_FILE)
clean_hts:
	rm -f $(HTSLIB)
clean_abpoa:
	rm -f $(ABPOA_LIB) $(ABPOA_DIR)/src/*.o
clean_wfa2:
	rm -f $(WFA2_LIB)
