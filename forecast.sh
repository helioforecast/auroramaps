#!/bin/bash

cd /usr/local/app
conda init
source ~/.bashrc
conda activate aurora1
python aurora.py
