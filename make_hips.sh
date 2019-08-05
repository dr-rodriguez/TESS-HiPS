#!/usr/bin/env bash

# # Removing live so fits aren't saved
# java -Xmx16g -jar AladinBeta.jar -hipsgen maxthread=20 in=Cutouts out=TESS-HiPS creator_did=ivo://STScI/TESS 
#hips_pixel_bitpix=16 hips_pixel_cut='100 8500 sqrt' (not needed with new scaling in cuts)

# Live version that allows appending
java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in="Cutouts/s0001" out="./TESS-HiPS" creator_did=ivo://STScI/TESS
# Add more sectors to it
java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in="Cutouts/s0002" out="./TESS-HiPS" creator_did=ivo://STScI/TESS APPEND
java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in="Cutouts/s0003" out="./TESS-HiPS" creator_did=ivo://STScI/TESS APPEND
java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in="Cutouts/s0004" out="./TESS-HiPS" creator_did=ivo://STScI/TESS APPEND
