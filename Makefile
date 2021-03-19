all: lib main

debug:
	make MODE=debug

lib:
	cd lib; $(MAKE) $@

main:
	cd main; $(MAKE) $@

dbg:
	cd main; $(MAKE) $@

clean:
	cd lib; $(MAKE) $@
	cd main; $(MAKE) $@
	rm -rf toolbox.exe

.PHONY:	lib main dbg
