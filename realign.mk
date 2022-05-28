## Makefile for boosting the alignment of sequences
## make -f src/realign.mk                  makes re-aligned sequences from INITALIGN
## make -f src/realign.mk KEEPX=''         turns off --keepx
## make -j3 -f src/realign.mk              permits parallel make
## make -j3 -f src/realign.mk ID=20211024  makes re-aligned sequences based on date 20211024
## make -f src/realign.mk Latest           makes a new Latest.ipkl

ID := 20220522
INITALIGN := data/REG_COMP.SPIKE.protein.Human.$(ID).fasta.xz
KEEPX := --keepx
ifeq ($(KEEPX),--keepx)
XID := xx$(ID)
else
XID := $(ID)
endif
XIDGZ := $(XID).fasta.gz

.DEFAULT_GOAL := ztkfx-$(XIDGZ)

.PHONY: Latest

fx-$(XIDGZ): $(INITALIGN)
	python -m fixalign $(KEEPX) -i $(INITALIGN) --fix $@

tkfx-$(XIDGZ): fx-$(XIDGZ)
	python -m tweakalign $(KEEPX) -i $<  -o $@

ztkfx-$(XIDGZ): tkfx-$(XIDGZ)
	python -m fixfasta $(KEEPX) --rmgapcols -i $< -o $@

Latest: Latest.ipkl

Latest.ipkl: $(.DEFAULT_GOAL)
	python -m fixfasta $(KEEPX) --random -i $< -o $@
