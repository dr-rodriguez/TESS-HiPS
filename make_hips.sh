#!/usr/bin/env bash

#java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in=Cutouts out=TESS-S1-HiPS creator_did=ivo://STScI/TESS
# Removing live so fits aren't saved
java -Xmx16g -jar AladinBeta.jar -hipsgen maxthread=20 in=Cutouts out=TESS-HiPS creator_did=ivo://STScI/TESS 
#hips_pixel_bitpix=16 hips_pixel_cut='100 8500 sqrt'

# Add images to it
# java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in=Cutouts out=TESS-S1-HiPS APPEND
