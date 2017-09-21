##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 08/04/2017: Modified.
#####	08/04/2017: Previously modified.
#####	08/04/2017: Created.
#####	Description: 
#####		- Bash script for launching Julia jobs.
#####		- Julia v0.5.1 parallel (pmap). 
#####	Modifications:
#####		08/04/2017:
#####			- Duplicated from p1a1-051.sh.
#####			- Modify script for running p2a1.jl.
##########################################################################

#!/bin/bash
##### Define local variables.
program="p2"
version="a1"
julia_version=051
##### Set datetime variable.
datetime=`date +%Y%m%d%H%M%S`
##### Execute program.
./jwrap${julia_version} ${program}${version}.jl > \
../03_Log/${program}${version}-${julia_version}_${datetime}.log 2>&1
