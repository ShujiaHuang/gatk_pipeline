# Copyright (C) 2013 DNAnexus, Inc.
#
# This file is part of gatk_pipeline (DNAnexus platform app).

all: bamtools samtools

bamtools:
	mkdir -p bamtools/build
	cd bamtools/build; cmake ..; make
	cp -a bamtools/bin/bamtools-* resources/usr/bin

samtools:
	make -C samtools samtools
	cp -a samtools/samtools resources/usr/bin/samtools

clean:
	rm -f resources/usr/bin/samtools resources/usr/bin/bamtools*

.PHONY: all bamtools samtools clean
