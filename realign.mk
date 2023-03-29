## Makefile for boosting the alignment of sequences
## make -f src/realign.mk                  makes re-aligned sequences from INITALIGN
## make -f src/realign.mk KEEPX=''         turns off --keepx
## make -j3 -f src/realign.mk              permits parallel make
## make -j3 -f src/realign.mk ID=20211024  makes re-aligned sequences based on date 20211024
## make -f src/realign.mk Latest           makes a new Latest.ipkl

ID := 20230326
INITALIGN := data/SPIKE.protein.Human.$(ID).fasta
KEEPX := --keepx
ifeq ($(KEEPX),--keepx)
XID := xx$(ID)
else
XID := $(ID)
endif
XIDF := $(XID).fasta
XIDGZ := $(XID).fasta.xz

.DEFAULT_GOAL := all

.PHONY: all clean compress Latest

fx-$(XIDF): $(INITALIGN)
	python -m fixalign $(KEEPX) -i $(INITALIGN) --fix $@

tweak-$(XID).out: fx-$(XIDF)
	echo "[Y145Q,H146-] [Y144-,H146Q]" > $@
	echo "[Y145K,H146-] [Y144-,H146K]" >> $@
	echo "[Y145Q,H146I,K147-] [Y144-,H146Q,K147I]" >> $@
	python -m matchfasta -i $< --showmut | perl -nle '/\b((.)(\d+)-,\+\3\2)\b/ and printf "[%s] []\n",$$1' | sort | uniq >> $@

tkfx-$(XIDF): fx-$(XIDF) tweak-$(XID).out
	parallel -k --header 2 --recstart '>' --blocksize 400M --pipe \
	python -m tweakalign -M tweak-$(XID).out -M . $(KEEPX) --jobno {#} -i - -o - < fx-$(XIDF) > $@

ztkfx-$(XIDF): tkfx-$(XIDF)
	python -m fixfasta $(KEEPX) --rmgapcols --random -i $< -o $@

Latest: Latest.fasta

Latest.fasta: ztkfx-$(XIDF)
	cp $< $@

## if done, compress intermediate fastas

fx-$(XIDGZ): fx-$(XIDF) tkfx-$(XIDF)
	xz -T8 $<

tkfx-$(XIDGZ): tkfx-$(XIDF) ztkfx-$(XIDF)
	xz -T8 $<

ztkfx-$(XIDGZ): ztkfx-$(XIDF) Latest.fasta
	xz -T8 $< 

compress: fx-$(XIDGZ) tkfx-$(XIDGZ) ztkfx-$(XIDGZ)

clean:
	/bin/rm -f fx-$(XIDF)       fx-$(XIDGZ)
	/bin/rm -f tkfx-$(XIDF)   tkfx-$(XIDGZ)
	/bin/rm -f ztkfx-$(XIDF) ztkfx-$(XIDGZ)
	/bin/rm -f tweak-$(XID).out

all: Latest.fasta compress
