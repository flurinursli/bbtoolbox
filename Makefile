all: lib main

debug: lib dbg

lib:
	cd lib; $(MAKE) $@

main:
	cd main; $(MAKE) $@

dbg:
	cd main; $(MAKE) $@

clean:
	cd lib; $(MAKE) $@
	cd main; $(MAKE) $@
	rm -rf inveta.exe
	rm -rf debug.exe

.PHONY:	lib main dbg
