CC      = g++
ifdef DEBUG
CFLAGS = -g
else
CFLAGS = -O3
endif
LIB_DIR = $(PWD)/../RooRarFit/tmp/
INCL_DIR = $(PWD)/../RooRarFit/tmp/RooRarFit
CFLAGS  += -Wall -std=c++11 -Wextra
CFLAGS  += -I $(INCL_DIR)
CFLAGS  += `root-config --cflags`
LDFLAGS = -L $(LIB_DIR) -lRooRarFit 
LDFLAGS += `root-config --glibs`
LDFLAGS += -lRooFit -lRooFitCore

OBJ = fit_table.cc fit_table.h

all: fit_table

# Turn off message: "Error in <TCling::RegisterModule>: cannot find dictionary module"
# ".pcm" should reside in LIB_DIR 
fit_table:  $(OBJ)
	$(CC) $(CFLAGS) -o fit_table fit_table.cc $(LDFLAGS) $(LIB)
	@ln -sf $(INCL_DIR)/RooRarFitCint_rdict.pcm $(LIB_DIR)

clean:
	\rm -f fit_table

show:
	echo "CFLAGS : $(CFLAGS)"
	echo "LDFLAGS: $(LDFLAGS)"
