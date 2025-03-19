#! /bin/bash

# Clean temp files
rm -rf tmp/*

# Clean debreak intermediary outputs (lots of files)
rm -rf data/GCA_*/debreak_minimap2/
rm -rf data/GCA_*/debreak_ngmlr/

#rm -rf data-lewontin/GCA_*/debreak_minimap2/
#rm -rf data-lewontin/GCA_*/debreak_ngmlr/
