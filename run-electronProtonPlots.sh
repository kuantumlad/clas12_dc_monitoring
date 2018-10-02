#!/bin/bash 

job=$1



root -l -b -q "quickFit.cxx+($job)"
echo " done with quick fit"

echo " starting electronProtonFinal "
root -l -b -q "electronProtonFinal.cxx+($job)"

echo " done with electron proton final "
echo " starting to make pdf file "

root -l -b -q "electronProtonAnaShort.cxx+($job)"

cp monitoring_job_pdf/monitoringJob${job}.pdf /u/home/bclary
