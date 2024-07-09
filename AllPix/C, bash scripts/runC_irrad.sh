#!/bin/sh
# Author: Shravan
# Iteratively runs the C file to get CCE data on the CONF files

echo "Run folder name: "
read RUNFOLDER
echo "Output CSV file name: "
read CSVFILE
CERNBOX="/eos/user/c/cmahajan/lxplus/allpix"
RUNPATH="$CERNBOX/$RUNFOLDER"
SOURCEFILE="/cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/setup.sh"
source $SOURCEFILE
for FILE in $RUNPATH/*.root;
    do
        if [ "$FILE" != "$RUNPATH/modules.root" ];  # Space between [ and "$FILE" is important
            then
                pattern="bias.*[0-9]+V"
                # Extracting the bias voltage number out of the filename
                [[ $FILE =~ $pattern ]]
                biasV=${BASH_REMATCH[0]}
                biasV="${biasV#"bias"}"
                biasV="${biasV%"V"}"
                root -b -l $FILE <<EOF
                .L /cvmfs/clicdp.cern.ch/software/allpix-squared/3.0.3/x86_64-centos7-gcc12-opt/lib/libAllpixObjects.so
                .L processData.C
                processData(_file0, "dut", "$CSVFILE", $biasV, 300, 7)
EOF
            fi
    done
# Assumes that all the .ROOT files in a folder are part of the analysis
# processData.C is in the same directory