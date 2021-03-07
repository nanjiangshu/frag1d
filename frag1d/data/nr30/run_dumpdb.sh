#!/bin/bash

my_formatdb.py -l train.idlist -dbname train_Qij -datapath train_Qij -dataext .Qijbin -format text
my_formatdb.py -l train.idlist -dbname train_modm -datapath train_modm -dataext .modmbin -format text
my_formatdb.py -l train.idlist -dbname train_fragacc -datapath train_fragacc -dataext .fragaccbin -format text
