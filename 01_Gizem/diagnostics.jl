##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 08/03/2017: Modified.
##### 08/02/2017: Previously modified.
##### 07/31/2017: Created.
##### Description: 
##### 	- Module to conduct diagnostic tests.
#####		- Refer to the_problem.jl in 
#####			https://github.com/rsarfati/DSGE-Private.jl/tree/pf/src/estimate
##### Modifications:
#####		08/01/2017: 
#####			- Work on setup().
#####		08/02/2017: 
#####			- Finish up with setup().
#####			- Work on main routine.
#####		08/03/2017:
#####			- Adjust display prompts and function parameters.
#####			- Change diagnose() function name to test1().
##########################################################################
##### Define module name.
module diagnostics
##### Load necessary modules.
using HDF5, JLD, DSGE;
using QuantEcon:solve_discrete_lyapunov;
##### Define objects to be exported.
export test1

########################################################################
##### EXTERNAL FUNCTION DEFINITIONS
########################################################################

########################################################################
##### Description: Function to execute mutation function.
#####	Input parameters:
#####		exectype: pmap, parfor, or serial.
#####		T: Number of time steps.
#####		distCall: Whether to use built-in distribution function or not.
#####	Output parameters:
#####		Cholesky of input matrix.
########################################################################
function test1(
	exectype="serial"::String,	
	T=50::Int64,
	distCall=true::Bool)
	##### Display prompt.
	println("\n**************************")
	println("***** Initiate Setup *****")
	println("**************************")
	##### Obtain parameters from setup function.
	m,system,TTT,sqrtS2,s0,P0,s_lag_tempered,ε,yt,nonmissing, N_MH,c, n_particles,deterministic, μ, cov_s,s_t_nontempered = setup();
	##### Display prompt.	
	println("\n****************************")
	println("***** Execute Mutation *****")
	println("****************************\n")
	
	println("***** Test parameters")
	println("*****   Julia version: ",VERSION)
	println("*****   Exec type: ",exectype)
	println("*****   # of workers: ",nworkers())
	println("*****   # of time steps: ",T)
	println("*****   Use built-in dist fnc: ",distCall,"\n")
	##### Begin timing for entire mutation routine.
	tic();	
	##### Loop over number of specified time steps.
	for t = 1:T       
		##### Begin timing for current time step.
		tic();
		##### Run mutation 10 times per time step.
		for i=1:10
			##### Initialize vector.
			acpt_vec=zeros(n_particles);
			##### Display prompt.
			print("Mutation ");
			##### Begin timing for current iteration.
			tic();
			##### Execute pmap if specified.
			if exectype=="pmap"
				##### Display prompt.
				print("(in parallel - pmap) ");
				##### Run pmap.
				out = pmap(i->mutation_problem(c,N_MH,deterministic,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing,distCall),1:n_particles);
			##### Execute parfor if specified.
			elseif exectype=="parfor"
				##### Display prompt.
				print("(in parallel - parfor) ");
				##### Run parallel forloop.
				out = @sync @parallel (hcat) for i=1:n_particles
					mutation_problem(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing,distCall)
				end
			##### Run sequentially if specified.
			elseif exectype=="serial"
				##### Display prompt.
				print("(not parallel) ")
				##### Run sequentially.
				out = [mutation_problem(c,N_MH,deterministic,system,yt,s_lag_tempered[:,i],ε[:,i],cov_s,nonmissing,distCall) for i=1:n_particles]
			end  
			##### Set prompt counter.
			procnt="Time step: "*string(t)*" Iteration: "*string(i)*" Elapsed time: "; 			
			##### Display prompt.
			println(procnt,round(toq(),2)," seconds");
			##### Disentangle three outputs of mutation and enter them into appropriate arrays.
			for j = 1:n_particles
				s_t_nontempered[:,j] = out[j][1]
				ε[:,j] = out[j][2]
				acpt_vec[j]=out[j][3]
			end
		end
		##### Display prompt.
		println("***** Completion of period ",t,"/",T)
		println("***** Elapsed time: ",round(toq(),2)," seconds\n")
	end   
	##### Display prompt.	
	println("*******************************")
	println("***** Execution Completed *****")
	println("*******************************\n")
	
	println("***** Test parameters")
	println("*****   Julia version: ",VERSION)
	println("*****   Exec type: ",exectype)
	println("*****   # of workers: ",nworkers())
	println("*****   # of time steps: ",T)
	println("*****   Use built-in dist fnc: ",distCall,"\n")
	
	println("***** Total elapsed time: ",round(toq(),2)," seconds")
##### Close function definition.
end

########################################################################
##### INTERNAL FUNCTION DEFINITIONS
########################################################################

########################################################################
##### Description: Function to calculate Cholesky of a matrix.
#####	Input parameters:
#####		mat: Matrix.
#####	Output parameters:
#####		Cholesky of input matrix.
########################################################################
function get_chol(
	mat::Array{Float64,2})
	##### Return Cholesky of matrix.
	return Matrix(chol(nearestSPD(mat)));
##### Close function definition.
end

########################################################################
##### Description: Function to return data.
#####	Input parameters:
#####		None.  Assume local load only needed as temporary workaround.
#####	Output parameters:
#####		data: Extracted data object.
########################################################################
function get_data()
	##### Define local data location.
	mydatapath="../02_Data"	
	##### Define dsge folder path.
	filesw="/data/dsge_data_dir/dsgejl/realtime/input_data/data";
	##### Define CSV filename.
	csvname="realtime_spec=smets_wouters_hp=true_vint=110110.csv";
	##### Attempt access to dsge folder.
	try
		##### Load CSV file from dsge folder.	
		data=readcsv("$filesw/"*csvname,header=true);
	##### Execute if error encountered.		
	catch err
		##### Print messages.		
		println("failed to access dsge path: $err");
		println("continuing with local load of csv file");
		##### Load CSV file from local location.
		data=readcsv(mydatapath*"/"*csvname,header=true);
	##### Close try/catch clause.
	end
##### Close function definition.
end

########################################################################
##### Description: Function to setup for calculation.
#####	Input parameters:
#####		None.
#####	Output parameters:
#####		Variables to be used for calculation.
########################################################################
function setup()
	##### Create new instance of model.
  m=SmetsWouters("ss1",testing=true);
	##### Get data.
	data=get_data();
	##### Display loaded data information.
	println("\nLoaded data type: ",typeof(data));
	##### Convert data into array format.
	data=convert(Array{Float64,2},data[1][:,2:end]);
	##### Transpose data.
	data=data';
	##### Set parameters.
	N_MH=10;
	c=0.1;
	n_particles=4000;
	##### Obtain path to system data.
	sysdatapath=Pkg.dir()*"/DSGE/test/reference/system.jld";
	##### Load system data.
	system=load(sysdatapath,"system");
	##### Set parameters.
	RRR=system.transition.RRR;
	TTT=system.transition.TTT;
	S2=system.measurement.QQ;
	sqrtS2=RRR*get_chol(S2)';
	s0=zeros(size(TTT)[1]);
	P0=nearestSPD(solve_discrete_lyapunov(TTT,RRR*S2*RRR'));
	n_errors=size(S2,1);
	n_states=size(system.measurement.ZZ,2);
	s_lag_tempered_rand_mat=randn(n_states,n_particles);
	ε=randn(n_errors, n_particles);	
	s_lag_tempered=repmat(s0,1,n_particles)+ 
	get_chol(P0)'*s_lag_tempered_rand_mat;
	yt=data[:,25];
	nonmissing=!isnan(yt);
	deterministic=false;
	μ=mean(ε,2);
	cov_s=(1/n_particles)*(ε-repmat(μ,1,n_particles))*(ε-repmat(μ,1,n_particles))'
	if !isposdef(cov_s)
		cov_s=diagm(diag(cov_s));
	end	
	s_t_nontempered=TTT*s_lag_tempered+sqrtS2*ε;
	##### Set return values.
	return m,system,TTT,sqrtS2,s0,P0,s_lag_tempered,ε,yt,nonmissing,
	N_MH,c,n_particles,deterministic,μ,cov_s, s_t_nontempered
##### Close function definition.
end

##### Close module definition.
end


