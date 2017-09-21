##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 08/04/2017: Modified.
##### 08/04/2017: Previously modified.
##### 08/04/2017: Created.
##### Description: 
##### 	- Program to run diagnostic tests using pmap batch.
##### Modifications:
#####		08/04/2017: 
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
##### Set program name.
program="p2"
version="a1"
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
println("***** Program name:  ",program,version)
println("***** Datetime:      ",Dates.today()," "
,Dates.format(now(),"HH:MM:SS"))
##########################################################################
##### Execute test code.
##########################################################################
##### Load cluster manager package.
using ClusterManagers;
##### Add workers.
addprocs_sge(30,queue="background.q");
##### Load diagnostic module.
using diagnostics2;
##### Create instance of mutation problem.
myprob=diagnostics2.the_problem();
##### Run pmap batch.
myprob.mutation_pmap_batch();
##########################################################################
##### Clean-up workspace.
##########################################################################
##### Remove temporary files.  Note: Shell expression "*" does not expand.
tempfiles=readdir(temppath)
for tempfile in tempfiles
	res=contains(tempfile,program*version)
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
println("***** Program name:  ",program,version)
println("***** Datetime:      ",Dates.today()," "
,Dates.format(now(),"HH:MM:SS"))
