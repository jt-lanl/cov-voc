## Makefile for boosting the alignment of sequences
## make -f src/realign.mk                  makes re-aligned sequences from INITALIGN
## make -f src/realign.mk KEEPX=''         turns off --keepx
## make -j3 -f src/realign.mk              permits parallel make
## make -j3 -f src/realign.mk ID=20211024  makes re-aligned sequences based on date 20211024
## make -f src/realign.mk Latest           makes a new Latest.ipkl

ID := 20211108
INITALIGN := data/REG_COMP.SPIKE.protein.Human.$(ID).fasta.gz
KEEPX := --keepx
ifeq ($(KEEPX),--keepx)
XID := xx$(ID)
else
XID := $(ID)
endif
XIDGZ := $(XID).fasta.gz

.DEFAULT_GOAL := ztkfx-sk-$(XIDGZ)

tweakfile.txt:
	echo '[E156G,F157-,R158-] [E156-,F157-,R158G]' > $@
	echo '[+214AAG,D215Y] [D215A,+215AGY]' >> $@
	echo '[Y144_,Y145-] [Y144-,Y145_]' >> $@

fx-$(XIDGZ): $(INITALIGN)
	python -m fixalign $(KEEPX) -i $(INITALIGN) --fix $@

tkfx-$(XIDGZ): fx-$(XIDGZ) tweakfile.txt
	python -m tweakalign $(KEEPX) -i $< -M tweakfile.txt  -o $@

sk-$(XIDGZ): $(INITALIGN)
	python -m xalign $(KEEPX) --dedash -i $(INITALIGN) -o $@

tksk-$(XIDGZ): sk-$(XIDGZ) tweakfile.txt
	python -m tweakalign $(KEEPX) -i $< -M tweakfile.txt -o $@

tkfx-sk-$(XID).msp: tksk-$(XIDGZ) tkfx-$(XIDGZ)
	python -m cfalign -a $^ --msp $@

tkfx-sk-$(XIDGZ): tkfx-sk-$(XID).msp tkfx-$(XIDGZ) 
	python -m tweakalign $(KEEPX) -M tkfx-sk-$(XID).msp -i tkfx-$(XIDGZ) -o $@

ztkfx-sk-$(XIDGZ): tkfx-sk-$(XIDGZ)
	python -m fixfasta $(KEEPX) --rmgapcols --random -i $< -o $@

Latest: Latest.ipkl

Latest.ipkl: ztkfx-sk-$(XIDGZ)
	python -m fixfasta $(KEEPX) -i $< -o $@

