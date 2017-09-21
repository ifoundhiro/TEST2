##########################################################################
##### User: Hirotaka Miura.                                   
##### Position: Research Analytics Associate.                                            
##### Organization: Federal Reserve Bank of New York.
##########################################################################   
##### 08/04/2017: Modified.
##### 08/04/2017: Previously modified.
##### 08/04/2017: Created.
##### Description: 
##### 	- Module to conduct diagnostic tests.
#####		- Refer to the_problem.jl in 
#####			https://github.com/rsarfati/DSGE-Private.jl/tree/pf/src/estimate
##### Modifications:
#####		08/04/2017:
#####			- Duplicated from diagnostics.jl.
#####			- Define type definition for mutation problem.
##########################################################################
##### Define module name.
module diagnostics2
##### Load necessary modules.
using HDF5, JLD, DSGE;
using QuantEcon:solve_discrete_lyapunov;
##### Define objects to be exported.
export the_problem

########################################################################
##### EXTERNAL DEFINITIONS
########################################################################

########################################################################
##### Description: Type definition for mutation problem.
########################################################################
type the_problem
	######################################################################
	##### Declare member variables and functions..
	######################################################################
	##### Declare member variables.
	m::DSGE.SmetsWouters{Float64}
	data0::Tuple
	data::Array{Float64,2}
	N_MH::Int64
	c::Float64
	n_particles::Int64
	sysdatapath::String
	system::DSGE.System{Float64}
	RRR::Array{Float64,2}
	TTT::Array{Float64,2}
	S2::Array{Float64,2}
	sqrtS2::Array{Float64,2}
	s0::Array{Float64,1}
	P0::Array{Float64,2}
	n_errors::Int64
	n_states::Int64
	s_lag_tempered_rand_mat::Array{Float64,2}
	ε::Array{Float64,2}
	s_lag_tempered::Array{Float64,2}
	yt::Array{Float64,1}
	nonmissing::Array{Bool,1}
	deterministic::Bool
	μ::Array{Float64,2}
	cov_s::Array{Float64,2}
	s_t_nontempered::Array{Float64,2}
	distCall::Bool
	T::Int64
	batchsize::Int64
	##### Declare member functions.
	mutation_wrapper::Function
	mutation_pmap_batch::Function
	######################################################################
	##### Define main function.
	######################################################################
	function the_problem()
		##### Create a new instance.
		this=new();
		##### Display prompt.
		println("\n**************************")
		println("***** Initiate Setup *****")
		println("**************************\n")
		##### Create new instance of model.
		this.m=SmetsWouters("ss1",testing=true);
		##### Get data.
		this.data0=get_data();
		##### Display loaded data information.
		println("\nLoaded data type: ",typeof(this.data0));
		##### Convert data into array format.
		this.data=convert(Array{Float64,2},this.data0[1][:,2:end]);
		##### Transpose data.
		this.data=this.data';
		##### Obtain path to system data.
		this.sysdatapath=Pkg.dir()*"/DSGE/test/reference/system.jld";
		##### Load system data.
		this.system=load(this.sysdatapath,"system");
		##### Set parameters.
		this.N_MH=10;
		this.c=0.1;
		this.n_particles=4000;		
		this.RRR=this.system.transition.RRR;
		this.TTT=this.system.transition.TTT;
		this.S2=this.system.measurement.QQ;
		this.sqrtS2=this.RRR*get_chol(this.S2)';
		this.s0=zeros(size(this.TTT)[1]);
		this.P0=nearestSPD(solve_discrete_lyapunov(
		this.TTT,this.RRR*this.S2*this.RRR'));
		this.n_errors=size(this.S2,1);
		this.n_states=size(this.system.measurement.ZZ,2);
		this.s_lag_tempered_rand_mat=randn(this.n_states,this.n_particles);
		this.ε=randn(this.n_errors,this.n_particles);	
		this.s_lag_tempered=repmat(this.s0,1,this.n_particles)+ 
		get_chol(this.P0)'*this.s_lag_tempered_rand_mat;
		this.yt=this.data[:,25];
		this.nonmissing=!isnan(this.yt);
		this.deterministic=false;
		this.μ=mean(this.ε,2);
		this.cov_s=(1/this.n_particles)*(this.ε-repmat(
		this.μ,1,this.n_particles))*(this.ε-repmat(
		this.μ,1,this.n_particles))'
		if !isposdef(this.cov_s)
			this.cov_s=diagm(diag(this.cov_s));
		end	
		this.s_t_nontempered=this.TTT*this.s_lag_tempered+this.sqrtS2*this.ε;		
		##### Set default execution values.
		this.T=50;
		this.distCall=true;
		##### Set pmap batch size.
		if mod(this.n_particles,nworkers())>0
			this.batchsize=round(this.n_particles/nworkers())+1;
		else
			this.batchsize=this.n_particles/nworkers();
		end
		####################################################################
		##### Define wrapper for mutation problem.
		####################################################################
		this.mutation_wrapper=function(i::Int64)
			##### Define mutation problem.
			mutation_problem(
				this.c,
				this.N_MH,
				this.deterministic,
				this.system,
				this.yt,
				this.s_lag_tempered[:,i],
				this.ε[:,i],
				this.cov_s,
				this.nonmissing,
				this.distCall
			);
		end
		####################################################################
		##### Define function to execute mutations using pmap batch.
		####################################################################
		this.mutation_pmap_batch=function()
			##### Display prompt.	
			println("\n****************************")
			println("***** Execute Mutation *****")
			println("****************************\n")
			println("***** Test parameters")
			println("*****   Julia version: ",VERSION)
			println("*****   Exec type: pmap batch")
			println("*****   # of workers: ",nworkers())
			println("*****   # of time steps: ",this.T)
			println("*****   Use built-in dist fnc: ",this.distCall,"\n")		
			##### Begin timing for entire mutation routine.
			tic();	
			##### Loop over number of specified time steps.
			for t = 1:this.T       
				##### Begin timing for current time step.
				tic();
				##### Run mutation 10 times per time step.
				for i=1:10
					##### Initialize vector.
					acpt_vec=zeros(this.n_particles);
					##### Display prompt.
					print("Mutation ");
					##### Begin timing for current iteration.
					tic();
					##### Execute pmap with batch.
					out=pmap(
						this.mutation_wrapper,
						1:this.n_particles,
						batch_size=this.batchsize
					);
					##### Set prompt counter.
					procnt="Time step: "*string(t)*" Iteration: "*
					string(i)*" Elapsed time: "; 			
					##### Display prompt.
					println(procnt,round(toq(),2)," seconds");
					##### Disentangle three outputs of mutation and enter them into appropriate arrays.
					for j = 1:this.n_particles
						this.s_t_nontempered[:,j] = out[j][1]
						this.ε[:,j] = out[j][2]
						acpt_vec[j]=out[j][3]
					end
				end
				##### Display prompt.
				println("***** Completion of period ",t,"/",this.T)
				println("***** Elapsed time: ",round(toq(),2)," seconds")
				println("***** s_t_nontempered min: ",
				round(minimum(this.s_t_nontempered),2)," max: ",
				round(maximum(this.s_t_nontempered),2)," mean: ",
				round(mean(this.s_t_nontempered),2)," sum: ",
				round(sum(this.s_t_nontempered),2))
				println("***** ε min: ",round(minimum(this.ε),2),
				" max: ",round(maximum(this.ε),2)," mean: ",
				round(mean(this.ε),2)," sum: ",round(sum(this.ε),2),"\n")		
			end   
			##### Display prompt.	
			println("*******************************")
			println("***** Execution Completed *****")
			println("*******************************\n")
			
			println("***** Test parameters")
			println("*****   Julia version: ",VERSION)
			println("*****   Exec type: pmap batch")
			println("*****   # of workers: ",nworkers())
			println("*****   # of time steps: ",this.T)
			println("*****   Use built-in dist fnc: ",this.distCall,"\n")
			
			println("***** Total elapsed time: ",round(toq(),2)," seconds")
		##### Close function definition.		
		end
		##### Return instance.
		return this
	##### Close function definition.	
	end
##### Close type definition.	
end

########################################################################
##### INTERNAL DEFINITIONS
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

##### Close module definition.
end


