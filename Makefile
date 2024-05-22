CC          = gcc
# add -fno-tree-vectorize to avoid certain vectorization errors in O3 optimization
# right now, we are using -O3 for the best performance, and no vectorization errors were found
EXTRA_FLAGS = -Wall -Wno-unused-function -Wno-misleading-indentation

HTSLIB_DIR  = ./htslib
HTSLIB      = $(HTSLIB_DIR)/libhts.a
LIB         = $(HTSLIB) -lm -lz -lpthread -llzma -lbz2 -lcurl
INCLUDE     = -I $(HTSLIB_DIR)

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
SOURCE  = $(wildcard ${SRC_DIR}/*.c) 
OBJS    = $(SOURCE:.c=.o)
BIN     = $(BIN_DIR)/longcallD
ifneq ($(gdb),)
	BIN = $(BIN_DIR)/gdb_longcallD
endif

.c.o:
	$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

all: $(HTS_ALL) $(BIN)

$(HTS_ALL):
	cd $(HTSLIB_DIR); make;

$(BIN): $(OBJS)
	if [ ! -d $(BIN_DIR) ]; then mkdir $(BIN_DIR); fi
	$(CC) $(OBJS) -o $@ $(LIB) $(PG_FLAG)

$(SRC_DIR)/bam_utils.o: $(SRC_DIR)/bam_utils.c $(SRC_DIR)/bam_utils.h $(SRC_DIR)/utils.h
$(SRC_DIR)/collect_snps.o: $(SRC_DIR)/collect_snps.c $(SRC_DIR)/collect_snps.h $(SRC_DIR)/bam_utils.h
$(SRC_DIR)/kalloc.o: $(SRC_DIR)/kalloc.c $(SRC_DIR)/kalloc.h
$(SRC_DIR)/kthread.o: $(SRC_DIR)/kthread.c
$(SRC_DIR)/main.o: $(SRC_DIR)/main.c $(SRC_DIR)/phase_bam.h
$(SRC_DIR)/phase_bam.o: $(SRC_DIR)/phase_bam.c $(SRC_DIR)/phase_bam.h $(SRC_DIR)/main.h $(SRC_DIR)/utils.h $(SRC_DIR)/seq.h \
                        $(SRC_DIR)/collect_snps.h
$(SRC_DIR)/seq.o: $(SRC_DIR)/seq.c $(SRC_DIR)/seq.h $(SRC_DIR)/utils.h
$(SRC_DIR)/utils.o: $(SRC_DIR)/utils.c $(SRC_DIR)/utils.h $(SRC_DIR)/ksort.h $(SRC_DIR)/kseq.h

.PHONY: clean

clean:
	rm -f $(SRC_DIR)/*.o $(BIN)