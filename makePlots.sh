#!/bin/csh -f
setenv file $1
setenv job $2
source /group/clas12/gemc/environment.csh 4a.2.4
#convert2GeV  $file job$job.root
root -l -b -q "electronProtonFinal.cxx+($job)"
root -l -b -q "electronProtonAnaShort.cxx($job)"
