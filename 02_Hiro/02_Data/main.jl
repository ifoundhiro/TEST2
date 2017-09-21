##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 08/03/2017: Modified.
##### 08/03/2017: Previously modified.
##### 08/03/2017: Created.
##### Description: 
##### 	- Program to run Julia diagnostic tests.
##### Modifications:
##### 		08/03/2017: 
#####			- Duplicated from p1a1.jl.
##########################################################################
##### Clear Julia.
workspace()
##### Set relative directory paths.
datapath="../02_Data"
logpath="../03_Log"
graphpath="../04_Graph"
latexpath="../05_Latex"
docpath="../06_Document"
temppath="../07_Temp"
validpath="../08_Validation"
##########################################################################
##### Display system information.
##########################################################################
println("\n******************************")
println("***** System Information *****")
println("******************************")
println("***** User:          ",ENV["USER"])
println("***** Julia version: ",VERSION)
println("***** Node:          ",gethostname())
println("***** Directory:     ",pwd())
println("***** Program name:  ",ARGS[1])
println("***** Datetime:      ",Dates.today()," "
,Dates.format(now(),"HH:MM:SS"))
##########################################################################
##### Extract and display command line arguments.
##########################################################################
##### Extract command line arguments.
i=1; progname=ARGS[i];
i+=1; exectype=ARGS[i];
i+=1; num_workers=parse(Int64,ARGS[i]);
i+=1; T=parse(Int64,ARGS[i]);
i+=1; distCall=(ARGS[i]=="true"); distCall_str=ARGS[i];
##### Display extracted arguments.
println("\n***** Command line arguments")
println("*****   Program name: ",progname)
println("*****   Exec type: ",exectype)
println("*****   # of workers: ",num_workers)
println("*****   # of time steps: ",T)
println("*****   Use built-in dist fnc: ",distCall,"\n")
##########################################################################
##### Error handle command line arguments.
##########################################################################
##### Check for valid execution type.
if exectype!="serial" && exectype!="parfor" && exectype!="pmap"
	error("invalid command line argument - execution type")
end
##### Check for valid distCall value.
if distCall_str!="true" && distCall_str!="false"
	error("invalid command line argument - distCall")
end
##### Check for mismatch between execution type and number of workers.
if exectype=="serial" && num_workers>1
	error("command line mismatch between exec type and # of workers")
elseif (exectype=="parfor" || exectype=="pmap") && num_workers==1
	error("command line mismatch between exec type and # of workers")
end
##########################################################################
##### Execute test code.
##########################################################################
##### Add workers if more than 1 specified.
if num_workers>1
	##### Load cluster manager package.
	using ClusterManagers;
	##### Use sge version of "addprocs()" to sync with scheduler.
	addprocs_sge(num_workers,queue="background.q");
end
##### Load necessary modules.
using diagnostics;
##### Execute diagnostic test.
diagnostics.test1(exectype,T,distCall)
##########################################################################
##### Clean-up workspace.
##########################################################################
##### Remove temporary files.  Note: Shell expression "*" does not expand.
tempfiles=readdir(temppath)
for tempfile in tempfiles
	res=contains(tempfile,progname)
	if res==true
		rm(temppath*"/"*tempfile)
	end
end
##########################################################################
##### Display system information.
##########################################################################
println("\n******************************")
println("***** System Information *****")
println("******************************")
println("***** User:          ",ENV["USER"])
println("***** Julia version: ",VERSION)
println("***** Node:          ",gethostname())
println("***** Directory:     ",pwd())
println("***** Program name:  ",progname)
println("***** Datetime:      ",Dates.today()," "
,Dates.format(now(),"HH:MM:SS"))
