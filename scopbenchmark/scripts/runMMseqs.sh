#!/bin/bash

## conda environment: plass 
conda activate plass

## mmseqs search module
mmseqs search /storage/databases/alphafold2/alphafolddb /storage/databases/alphafold2/alphafolddb /data3/stephanie/prj-foldseek/alphafoldDB/mmseqs/mmseqs.m8 tmp -a 

## mmseqs convertalis module

