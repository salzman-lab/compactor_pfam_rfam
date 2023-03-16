#!/bin/sh


Rfam=/oak/stanford/groups/horence/zheludev/preprint_060722/pipeline/Rfam


cmscan --notextw --rfam --cut_ga --cpu 32 \
    --tblout RFAM.tblout --fmt 2 --clanin $Rfam/Rfam.clanin $Rfam/Rfam.cm $1 > RFAM.stdout