#!/usr/bin/env bash

java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in=Cutouts out=TESS-S1-HiPS creator_did=ivo://STScI/TESS

# Add images to it
# java -Xmx16g -jar AladinBeta.jar -hipsgen -live maxthread=20 in=Cutouts out=TESS-S1-HiPS APPEND
