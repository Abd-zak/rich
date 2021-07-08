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
SOFLAGS = -fPIC

##### CSEventData
CECLS = CSEventData
CEHDR = $(CECLS).h
CESRC = $(CECLS).cc
CEOBJ = $(CECLS).o 
CELINKDEF = $(CECLS)LinkDef.h
CEDICTC = dict$(CECLS).cxx
CEDICTO = dict$(CECLS).o
CSEVENT = libCSEvent.so
CSLIB = -L$(PWD) -lCSEvent

OBJ = fit_table.cc fit_table.h

all: fit_table

# Turn off message: "Error in <TCling::RegisterModule>: cannot find dictionary module"
# ".pcm" should reside in LIB_DIR 
fit_table:  $(OBJ) $(CSEVENT)
	$(CC) $(CFLAGS) -o fit_table fit_table.cc $(LDFLAGS) $(CSLIB) $(LIB)
	@ln -sf $(INCL_DIR)/RooRarFitCint_rdict.pcm $(LIB_DIR)

$(CSEVENT): $(CEOBJ) $(CEDICTO)
	g++ -shared $(CFLAGS) $(SOFLAGS) -o $(CSEVENT) $(CEOBJ) $(CEDICTO)

$(CEOBJ): $(CEHDR) $(CESRC) $(PHAST_INCFILE)
	g++ -c $(CFLAGS) $(SOFLAGS) $(PHAST_INCL) $(CESRC)

$(CEDICTO): $(CEHDR) $(CELINKDEF)
	rootcint -f $(CEDICTC) -c $(CEHDR) $(CELINKDEF)
	g++ -c $(CFLAGS) $(SOFLAGS) $(CEDICTC)

clean:
	\rm -f fit_table *.o

show:
	@echo "CFLAGS : $(CFLAGS)"
	@echo "LDFLAGS: $(LDFLAGS)"
