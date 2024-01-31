## Makefile for boosting the alignment of sequences
## make -f src/realign.mk                  makes re-aligned sequences from INITALIGN
## make -f src/realign.mk KEEPX=''         turns off --keepx
## make -j3 -f src/realign.mk              permits parallel make
## make -j3 -f src/realign.mk ID=20211024  makes re-aligned sequences based on date 20211024
## make -f src/realign.mk Latest           makes a new Latest.ipkl

MAKEFLAGS := --jobs=3

ID := 20240124
#INITALIGN := data/SPIKE.protein.$(ID).fasta
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

# get mstring pairs associated with inconsistent alignment
# for now, this is advisory
msp-$(XID).out: $(INITALIGN)
	python -m fixalign $(KEEPX) -i $< -M $@

fx-$(XIDF): $(INITALIGN)
	cat $< | parallel -k --header '>[^>]*' --recstart '>' --blocksize 400M --pipe python -m fixfasta $(KEEPX) -i - --jobno {#} -o - > $@

# idea here is to fix [...,A99-,+99A,...], where A is deleted and then inserted
# shouldn't happen but we want to fix it if it does
tweak1-$(XID).out: $(INITALIGN)
	cat $< | parallel -k --header '>[^>]*' --recstart '>' --blocksize 400M --pipe python -m matchfasta -i - --showmut --jobno {#} -o - | perl -nle '/\b((.)(\d+)-,\+\3\2)\b/ and printf "[%s] []\n",$$1' | sort -u > $@

## here, the tweaks are of the form: [V143-,+143R]->[V143R]
## not sure if that's always an improvement, but it does reduce the number of terms in the mstring
tweak2-$(XID).out: $(INITALIGN)
	cat $< | parallel -k --header '>[^>]*' --recstart '>' --blocksize 400M --pipe python -m matchfasta -i - --showmut --jobno {#} -o - | perl -nle '/\b((.)(\d+)-,\+\3(.))\b/ and printf "[%s] [%s%s%s]\n",$$1,$$2,$$3,$$4' | sort -u > $@

tweak-$(XID).out: tweak1-$(XID).out tweak2-$(XID).out
	cat tweak1-$(XID).out tweak2-$(XID).out > $@

tkfx-$(XIDF): fx-$(XIDF) tweak-$(XID).out
	parallel -k --header 2 --recstart '>' --blocksize 400M --pipe \
	python -m tweakalign -M tweak-$(XID).out -M . $(KEEPX) --jobno {#} -i - -o - < fx-$(XIDF) > $@

ztkfx-$(XIDF): tkfx-$(XIDF)
	python -m fixfasta $(KEEPX) --rmgapcols --random -i $< -o $@

ztkfx-nox-$(XIDF): tkfx-$(XIDF)
	python -m fixfasta --skipx --rmgapcols --random -i $< -o $@

beep: Latest
	beep

Latest: Latest.fasta Latest-nox.fasta

Latest.fasta: ztkfx-$(XIDF)
	cp $< $@

Latest-nox.fasta: ztkfx-nox-$(XIDF)
	cp $< $@
	/bin/rm ztkfx-nox-$(XIDF)

## if done, compress intermediate fastas

fx-$(XIDGZ): fx-$(XIDF) tkfx-$(XIDF)
	xz -T8 $<

tkfx-$(XIDGZ): tkfx-$(XIDF) ztkfx-$(XIDF)
	xz -T8 $<

ztkfx-$(XIDGZ): ztkfx-$(XIDF) Latest.fasta
	xz -T8 $< 

compress: fx-$(XIDGZ) tkfx-$(XIDGZ) ztkfx-$(XIDGZ)

clean:
	/bin/rm -f fx-$(XIDF)   fx-$(XIDGZ)
	/bin/rm -f tkfx-$(XIDF)   tkfx-$(XIDGZ)
	/bin/rm -f ztkfx-$(XIDF) ztkfx-$(XIDGZ)
	/bin/rm -f tweak*$(XID).out
	/bin/rm -f ztkfx-nox-$(XIDF)

all: Latest compress
