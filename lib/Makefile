lib : common llsq fft filter interpolation rtt rng geometry

common:
	cd $@; $(MAKE) $@ $<

filter:
	cd $@; $(MAKE) $@ $<

fft:
	cd $@; $(MAKE) $@ $<

llsq:
	cd $@; $(MAKE) $@ $<

interpolation:
	cd $@; $(MAKE) $@ $<

rtt:
	cd $@; $(MAKE) $@ $<

rng:
	cd $@; $(MAKE) $@ $<

geometry:
	cd $@; $(MAKE) $@ $<

clean:
	cd common; $(MAKE) $@
	cd filter; $(MAKE) $@
	cd fft; $(MAKE) $@
	cd llsq; $(MAKE) $@
	cd interpolation; $(MAKE) $@
	cd rtt; $(MAKE) $@
	cd rng; $(MAKE) $@
	cd geometry; $(MAKE) $@
	rm -rf *.mod
	rm -rf *.o

.PHONY : common filter fft llsq interpolation rtt rng geometry clean
