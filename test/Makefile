CFLAGS += -I../include
LDFLAGS += -L../src

profit: main.o
	$(CC) -o $@ $^ $(LDFLAGS) -lprofit -lgsl -lgslcblas -lRmath
