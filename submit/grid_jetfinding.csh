#!/bin/csh

# used to submit sequential jobs on the grid

# first make sure program is updated and exists
 make bin/jetFindAnalysis || exit

set ExecPath = `pwd`
set analysis = $1
set execute = './bin/jetFindAnalysis'
set exponent = 3
set xmldir = /nfs/rhi/STAR/software/pythia8/share/Pythia8/xmldoc

# Now Submit jobs for each data file
set i = 0
while ( $i < 20 )

# Create the output file base name
set OutBase = out/outfile_${i}

# Make the output names and path

set outName = ${OutBase}.root


# Logfiles. Thanks cshell for this "elegant" syntax to split err and out
set LogFile     = log/jetAnalysis_${i}.log
set ErrFile     = log/jetAnalysis_${i}.err

echo "Logging output to " $LogFile
echo "Logging errors to " $ErrFile

set arg = "$xmldir $exponent $outName"

qsub -V -q erhiq -l mem=2GB -o $LogFile -e $ErrFile -N jetfinderAnalysis -- ${ExecPath}/submit/qwrap.sh ${ExecPath} $execute $arg

@ i++

end
