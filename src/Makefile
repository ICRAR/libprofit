CFLAGS += -I../include -fPIC -DPROFIT_VERSION='"$(shell cat ../VERSION)"'
LDFLAGS += -L. -lm
OBJS = sersic.o profit.o sky.o psf.o convolve.o

.PHONY: all clean
all: libprofit.so profit-cli

clean:
	rm -f *.o libprofit.so profit-cli

libprofit.so: $(OBJS)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

profit-cli: profit-cli.o libprofit.so
	$(CC)  profit-cli.o -o $@ -lprofit $(LDFLAGS)
