#!/bin/bash

module load java 

java -jar dxCompiler-2.11.4.jar compile $1 -extras extraOptions.json -project project-GJbvyPjJy3Gy01jz4x8bXzgv -folder /deep_phewas/$2 -streamFiles all -f
