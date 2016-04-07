CFLAGS += -I../include -fPIC
LDFLAGS += -lm
OBJS = sersic.o profit.o sky.o

.PHONY: all clean
all: libprofit.so

clean:
	rm -f *.o libprofit.so

libprofit.so: $(OBJS)
	$(LD) -shared $(LDFLAGS) -o $@ $^
