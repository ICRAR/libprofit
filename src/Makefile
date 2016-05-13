CFLAGS += -I../include -fPIC
LDFLAGS += -lm
OBJS = sersic.o profit.o sky.o psf.o convolve.o

.PHONY: all clean
all: libprofit.so profit-cli

clean:
	rm -f *.o libprofit.so profit-cli

libprofit.so: $(OBJS)
	$(LD) -shared $(LDFLAGS) -o $@ $^

profit-cli: profit-cli.o libprofit.so
	$(CC) $(CFLAGS) profit-cli.o -o $@ -lprofit $(LDFLAGS)
