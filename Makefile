#master makefile for iemmod

all: mkdirs
	cd iemmod		; $(MAKE) install
	cd ascii2bin	; $(MAKE) install

mkdirs:
	-mkdir -p bin

clean:
	cd iemmod		; $(MAKE) $@
	cd ascii2bin	; $(MAKE) $@

realclean:
	cd iemmod		; $(MAKE) $@
	cd ascii2bin	; $(MAKE) $@
	rm -f bin/*
