# Check if the OS is macOS or linux
UNAME_S := $(shell uname -s)


# add -fno-tree-vectorize to avoid certain vectorization errors in O3 optimization
# right now, we are using -O3 for the best performance, and no vectorization errors were found
EXTRA_FLAGS = -Wall -Wno-unused-function -Wno-misleading-indentation -Wno-unused-variable -Wno-alloc-size-larger-than

# Define the version number
LONGCALLD_VERSION =0.0.4
# Get the Git commit hash
GIT_COMMIT := $(shell git rev-parse --short HEAD 2> /dev/null)
ifneq ($(GIT_COMMIT),)
	LONGCALLD_VERSION = 0.0.4-$(GIT_COMMIT)
endif

HTSLIB_DIR  = ./htslib
HTSLIB      = $(HTSLIB_DIR)/libhts.a

ABPOA_DIR   = ./abPOA
ABPOA_LIB   = ./abPOA/lib/libabpoa.a
ABPOA_INC_DIR = $(ABPOA_DIR)/include

WFA2_DIR    = ./WFA2-lib
WFA2_LIB    = $(WFA2_DIR)/lib/libwfa.a

LIB_PATH    =
LIB         = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) $(LIB_PATH) -lm -lz -lpthread -llzma -lbz2 -lcurl
INCLUDE     = -I $(HTSLIB_DIR) -I $(ABPOA_INC_DIR) -I $(WFA2_DIR)

# Try linking against libdeflate
ifeq ($(shell echo "int main() {return 0;}" | ${CC} -x c - $(LIB_PATH) -ldeflate >/dev/null 2>&1 && echo "yes"),yes)
	LIB += -ldeflate
endif

ifeq ($(UNAME_S),Linux) # Linux
	LIB += -lcrypto
	ifneq ($(portable),)
		LIB += -static-libgcc -static-libstdc++
		ifneq ($(opt_lib),)
			LIB = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) -static-libgcc -static-libstdc++ -L${opt_lib} -lm -lz -lpthread -llzma -lbz2 -lcurl -lssl -lcrypto -lssh2 -ldeflate -lzstd 
		else
			ifneq ($(OPT_LIB),)
				LIB = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) -static-libgcc -static-libstdc++ -L${OPT_LIB} -lm -lz -lpthread -llzma -lbz2 -lcurl -lssl -lcrypto -lssh2 -ldeflate -lzstd 
			endif
		endif
	endif
endif

# for debug
ifneq ($(debug),)
	EXTRA_FLAGS  += -D __DEBUG__
endif

ABPOA_GDB_LIB = ./abPOA/lib/libabpoa_sse41.a
ABPOA_NOR_LIB     = ./abPOA/lib/libabpoa.a
# for gdb
ifneq ($(gdb),)
	OPT_FLAGS = -g
	ABPOA_LIB = $(ABPOA_GDB_LIB)
else
	OPT_FLAGS = -O3
	ABPOA_LIB = $(ABPOA_NOR_LIB)
endif

CFLAGS = $(OPT_FLAGS) $(EXTRA_FLAGS) -DLONGCALLD_VERSION=\"$(LONGCALLD_VERSION)\"

# for gprof
ifneq ($(pg),)
	PG_FLAG  = -pg
	CFLAGS  += -pg
endif

ifneq ($(PREFIX),)
	OUT_PRE_DIR = $(PREFIX)
else
	OUT_PRE_DIR = .
endif

BIN_DIR = $(OUT_PRE_DIR)/bin
INC_DIR = ./include
SRC_DIR = ./src

HTS_ALL = hts_all
ABPOA_ALL = abpoa_all
WFA2_ALL = wfa2_all
CSOURCE    = $(wildcard ${SRC_DIR}/*.c) 
OBJS    = $(CSOURCE:.c=.o)
BIN     = $(BIN_DIR)/longcallD
ifneq ($(gdb),)
	BIN = $(BIN_DIR)/gdb_longcallD
endif

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all: $(HTS_ALL) $(ABPOA_LIB) $(WFA2_LIB) $(BIN)

$(HTS_ALL): $(HTSLIB)

$(HTSLIB): $(HTSLIB_DIR)/configure.ac
	cd $(HTSLIB_DIR); autoreconf -i; ./configure; make CC=${CC}

$(ABPOA_GDB_LIB): 
	cd $(ABPOA_DIR); make libabpoa gdb=1 sse41=1
$(ABPOA_NOR_LIB):
	cd $(ABPOA_DIR); make libabpoa

$(ABPOA_ALL): $(ABPOA_LIB)

$(WFA2_LIB):
	cd $(WFA2_DIR); make setup lib_wfa
$(WFA2_ALL): $(WFA2_LIB)


$(BIN): $(OBJS)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(OBJS) -o $@ $(LIB) $(PG_FLAG)

$(SRC_DIR)/align.o: $(SRC_DIR)/align.c $(SRC_DIR)/align.h $(SRC_DIR)/utils.h $(SRC_DIR)/bam_utils.h $(SRC_DIR)/seq.h
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

.PHONY: all hts_all abpoa_all wfa2_all clean clean_all clean_hts clean_abpoa clean_wfa2

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)
clean_all:
	rm -f $(SRC_DIR)/*.o $(BIN) $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB)
clean_hts:
	rm -f $(HTSLIB)
clean_abpoa:
	rm -f $(ABPOA_LIB)
clean_wfa2:
	rm -f $(WFA2_LIB)
