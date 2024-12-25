
#####################################
## DEBUG MODE
## Just more printing
#####################################
#DEBUG_MODE=-D_DEBUG

#####################################
## Classify only this model.
## For debugging purposes.
#####################################
#PROCESS_ONLY='-D_PROCESS_ONLY="ENSBL"'

#####################################
## TRANSFORM KMEANS CLASSES TO KOPPEN EXTENDED
## Uncomment for making the transformation.
## The class centroids change a lot in Kmeans, especially on precipitation.
## The final colors obtained are the five colors (classes) of original Koppen,
## and not the expected 28 classes of Koppen extended.
#####################################
#TRANSFORM=-DTRANSFORM_KMEANS_TO_KOPPEN

#####################################
## Use Principal Components for KMeans classification instead of weather variables directly.
#####################################
USE_PC=-DKMEANS_USE_PC

#####################################
## Use OpenMP to speedup
#####################################
USE_OMP=-fopenmp

#################################
#################################

BIL_TARGET = bil.out
CLA_TARGET = classifier.out

OPT_FLAGS=${DEBUG_MODE} ${PROCESS_ONLY} ${TRANSFORM} ${USE_PC}

CC      = gcc
CFLAGS  = -std=c99 -Wall -O3 -g -Isrc/ ${OPT_FLAGS}
LFLAGS  = -L/home/andres/lapack-3.12.0/ -lnetcdf -lm -llapack -lrefblas -lgfortran ${USE_OMP}

SRCDIR  = src
BINDIR  = bin

BIL_SOURCES := bil.c
CLA_SOURCES := classifier.c koppen.c

# Add prefixes
BIL_SOURCES := $(addprefix $(SRCDIR)/, $(BIL_SOURCES))
CLA_SOURCES := $(addprefix $(SRCDIR)/, $(CLA_SOURCES))

rm      = rm -f

all: $(BINDIR)/$(BIL_TARGET) $(BINDIR)/$(CLA_TARGET)

$(BINDIR)/$(BIL_TARGET): $(BIL_SOURCES)
	$(CC) $(CFLAGS) $(BIL_SOURCES) $(LFLAGS) -o $@

$(BINDIR)/$(CLA_TARGET): $(CLA_SOURCES)
	$(CC) $(CFLAGS) $(CLA_SOURCES) $(LFLAGS) -o $@

.PHONY: clean
clean:
	@$(rm) $(BINDIR)/$(BIL_TARGET)
	@$(rm) $(BINDIR)/$(CLA_TARGET)
