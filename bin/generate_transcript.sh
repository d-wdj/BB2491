#!/bin/bash
awk ' $3 == "transcript" { print $0 } ' ../reference/Homo_sapiens.GRCh38.90.gtf > ../reference/transcript.gtf

###ARCHIVED 2018.01.15
