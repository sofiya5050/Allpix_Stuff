#!/bin/sh
# Author: Shravan
# Iteratively runs the CONF file for a voltage sweep from 50V-900V
	# pull config files from CERNbox	
echo "Run folder name: "
read RUNFOLDER
echo "Config file name: "
read FILENAME
CERNBOX="/eos/user/c/cmahajan/lxplus/allpix"
CONFPATH="$CERNBOX/$RUNFOLDER/$FILENAME"
SOURCEFILE="/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/setup.sh"
	# modifies CONF file and runs it
for i in $(seq 0 8);
    do
		source $SOURCEFILE
        biasV=$((i*50))
			# Modifying bias voltage in CONF file
		echo "Running python script to change CONF file."
        python confchange.py $CONFPATH $biasV
			# Running CONF file
		echo "Running $CONFPATH with bias voltage $biasV V."
			# run scripts + specify output to CERNbox
		OUTPUTDIR="$CERNBOX/$RUNFOLDER"
		allpix -c $CONFPATH -o output_directory=$OUTPUTDIR
		echo "ROOT file created."
    done
unset SOURCEFILE
unset RUNFOLDER
unset FILENAME
unset OUTPUTDIR
#Doesn't create any folders.
#Assumes that the ROOT filename is in the form of {anything}bias{text with no numbers}{bias voltage}V{anything}
#Doesn't change the ROOT filename except the bias voltage part.
#The bias voltages for the sweep are hardcoded, with steps of 50V.
