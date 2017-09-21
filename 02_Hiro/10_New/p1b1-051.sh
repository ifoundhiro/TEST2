##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 08/03/2017: Modified.
#####	08/03/2017: Previously modified.
#####	08/03/2017: Created.
#####	Description: 
#####		- Bash script for launching Julia jobs.
#####		- Refer to "main.jl" for command line arguments.
#####		- Julia v0.5.1 parallel (parfor) with 30 workers. 
#####	Modifications:
##########################################################################

#!/bin/bash
##### Define local variables.
program="p1"
version="b1"
julia_version=051
exectype="parfor"
num_workers=30
T=50
distCall="true"
##### Set datetime variable.
datetime=`date +%Y%m%d%H%M%S`
##### Execute program.
./jwrap${julia_version} main.jl \
${program}${version} ${exectype} ${num_workers} ${T} ${distCall} > \
../03_Log/${program}${version}-${julia_version}_${datetime}.log 2>&1
