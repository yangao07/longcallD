CC          = gcc
CXX         = g++
# add -fno-tree-vectorize to avoid certain vectorization errors in O3 optimization
# right now, we are using -O3 for the best performance, and no vectorization errors were found
EXTRA_FLAGS = -Wall -Wno-unused-function -Wno-misleading-indentation

HTSLIB_DIR  = ./htslib
HTSLIB      = $(HTSLIB_DIR)/libhts.a

EDLIB_DIR     = ./edlib
EDLIB_INC_DIR = ./edlib/include
EDLIB         = $(EDLIB_DIR)/src/edlib.o

ABPOA_DIR   = ./abPOA
ABPOA_LIB   = ./lib/libabpoa.a
ABPOA_INC_DIR = $(ABPOA_DIR)/include

WFA2_DIR    = ./WFA2-lib
WFA2_LIB    = $(WFA2_DIR)/lib/libwfa.a
LIB         = $(HTSLIB) $(ABPOA_LIB) $(WFA2_LIB) -lm -lz -lpthread #-llzma -lbz2 -lcurl
INCLUDE     = -I $(HTSLIB_DIR) -I $(EDLIB_INC_DIR) -I $(ABPOA_INC_DIR) -I $(WFA2_DIR)

# for debug
ifneq ($(debug),)
	EXTRA_FLAGS  += -D __DEBUG__
endif

# for gdb
ifneq ($(gdb),)
	OPT_FLAGS = -g
else
	OPT_FLAGS = -O3
endif

CFLAGS = $(OPT_FLAGS) $(EXTRA_FLAGS)

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
CSOURCE    = $(wildcard ${SRC_DIR}/*.c) 
CPPSOURCE  = $(wildcard $(SRC_DIR)/*.cpp)
CPPSOURCE += $(EDLIB_DIR)/src/edlib.cpp
OBJS    = $(CSOURCE:.c=.o) $(CPPSOURCE:.cpp=.o)
BIN     = $(BIN_DIR)/longcallD
ifneq ($(gdb),)
	BIN = $(BIN_DIR)/gdb_longcallD
endif

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all: $(HTS_ALL) $(EDLIB) $(ABPOA_LIB) $(WFA2_LIB) $(BIN)

# disable lzma, bz2 (CRAM), and libcurl (network protocol support)
$(HTS_ALL): $(HTSLIB)

$(HTSLIB): $(HTSLIB_DIR)/configure.ac
	cd $(HTSLIB_DIR); autoreconf -i; ./configure --disable-lzma --disable-bz2 --disable-libcurl --without-libdeflate; make;

# edlib
$(EDLIB): $(EDLIB_DIR)/src/edlib.cpp $(EDLIB_DIR)/include/edlib.h
	$(CXX) $(CFLAGS) -c $< $(INCLUDE) -o $@

$(ABPOA_LIB): 
	cd $(ABPOA_DIR); make libabpoa PREFIX=$(PWD) CC=gcc

$(WFA2_LIB):
	cd $(WFA2_DIR); make CC=gcc

$(BIN): $(OBJS)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CXX) $(OBJS) -o $@ $(LIB) $(PG_FLAG)

$(SRC_DIR)/align.o: $(SRC_DIR)/align.c $(SRC_DIR)/align.h $(SRC_DIR)/utils.h
	$(CC) -c $(CFLAGS) $< $(INCLUDE) -o $@

$(SRC_DIR)/assign_aln_hap.o: $(SRC_DIR)/assign_aln_hap.c $(SRC_DIR)/assign_aln_hap.h $(SRC_DIR)/utils.h $(SRC_DIR)/bam_utils.h
$(SRC_DIR)/bam_utils.o: $(SRC_DIR)/bam_utils.c $(SRC_DIR)/bam_utils.h $(SRC_DIR)/utils.h
$(SRC_DIR)/cgranges.o: $(SRC_DIR)/cgranges.c $(SRC_DIR)/cgranges.h $(SRC_DIR)/khash.h
$(SRC_DIR)/collect_var.o: $(SRC_DIR)/collect_var.c $(SRC_DIR)/collect_var.h $(SRC_DIR)/bam_utils.h
$(SRC_DIR)/kalloc.o: $(SRC_DIR)/kalloc.c $(SRC_DIR)/kalloc.h
$(SRC_DIR)/kthread.o: $(SRC_DIR)/kthread.c
$(SRC_DIR)/main.o: $(SRC_DIR)/main.c $(SRC_DIR)/call_var.h
$(SRC_DIR)/call_var.o: $(SRC_DIR)/bam_utils.c $(SRC_DIR)/call_var.c $(SRC_DIR)/call_var.h $(SRC_DIR)/main.h $(SRC_DIR)/utils.h $(SRC_DIR)/seq.h \
                        $(SRC_DIR)/collect_var.h
$(SRC_DIR)/phase_based_call_var.o: $(SRC_DIR)/phase_based_call_var.c $(SRC_DIR)/phase_based_call_var.h $(SRC_DIR)/align.h $(SRC_DIR)/utils.h
$(SRC_DIR)/seq.o: $(SRC_DIR)/seq.c $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h
$(SRC_DIR)/utils.o: $(SRC_DIR)/utils.c $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h $(SRC_DIR)/kseq.h
$(SRC_DIR)/vcf_utils.o: $(SRC_DIR)/vcf_utils.c $(SRC_DIR)/vcf_utils.h $(SRC_DIR)/utils.h

.PHONY: clean

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)