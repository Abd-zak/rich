CC      = g++
CFLAGS  = -Wall -Wextra -O2
CFLAGS  += `root-config --cflags`
LDFLAGS = -L $(PWD)/../RooRarFit/tmp -lRooRarFit 
LDFLAGS += `root-config --glibs`
LDFLAGS += -lRooFit

all: draw_rich_fig create_table

draw_rich_fig: draw_rich_fig.cc
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS)  -lHistPainter

create_table: create_table.cc
	$(CC) $(CFLAGS) -o $@ $< $(LDFLAGS) 

show:
	echo "CFLAGS : $(CFLAGS)"
	echo "LDFLAGS: $(LDFLAGS)"
