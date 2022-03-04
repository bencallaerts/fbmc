## Flow Based Market Coupling
# author: Ben Callaerts
# last update: March 3, 2022

## Step 0: Activate environment
cd("C:/Users/User/Documents/Thesis/Modeling")

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using CSV
using DataFrames
using YAML
using PrettyTables

##  Step 1: Input data
data = "data-basic"

bus = CSV.read(joinpath(data, "bus.csv"), DataFrame; delim=";")
branch = CSV.read(joinpath(data, "branch.csv"), DataFrame; delim=";")
plant = CSV.read(joinpath(data, "plant.csv"), DataFrame; delim=";")
load = CSV.read(joinpath(data, "load.csv"), DataFrame; delim=";")
incidence = CSV.read(joinpath(data, "incidence.csv"), DataFrame; delim=";")
susceptance = CSV.read(joinpath(data, "susceptance.csv"), DataFrame; delim=";")
@info "Data read"

## Step 2: Create model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))
@info "Model created"

# Define sets
function define_sets!(m::Model)
    m.ext[:sets] = Dict()

    N = m.ext[:sets][:N] = bus[:,:BusID]
    L = m.ext[:sets][:L] = branch[:,:BranchID]
    P = m.ext[:sets][:P] = plant[:,:GenID]
    Z = m.ext[:sets][:Z] = sort(unique(bus[:,:Zone]))

    return N, L, P, Z
end

# Build GSK
function get_gsk_flat()
	gsk_temp = zeros(Float64, length(N), length(Z))
	for n in N
		zone_temp = bus.Zone[bus[:,:BusID].==n][1]
		gsk_value_temp = 1/size(bus[bus[:,:Zone].==zone_temp,:],1)
		gsk_temp[findfirst(N .== n), findfirst(Z .== zone_temp)] = gsk_value_temp
	end
	return gsk_temp
end

# Build PTDF
function build_ptdf(m::Model)
	N = m.ext[:sets][:N]
	L = m.ext[:sets][:L]
	Z = m.ext[:sets][:Z]
	gsk = m.ext[:parameters][:gsk]

	MWBase = 380^2
	slack_node = 4
	slack_position = findfirst(N .== slack_node)

	# Build nodal PTDFs
	line_sus_mat = Matrix(susceptance)/MWBase*Matrix(incidence)
	node_sus_mat = transpose(Matrix(incidence))*Matrix(susceptance)/MWBase*Matrix(incidence)

	line_sus_mat_ = line_sus_mat[:, 1:end .!= slack_position]
	node_sus_mat_ = node_sus_mat[1:end .!= slack_position, 1:end .!= slack_position]

	ptdfN = line_sus_mat_*inv(node_sus_mat_)
	zero_column = zeros(Float64, length(L), 1)
	ptdfN = hcat(ptdfN[:,1:(slack_position-1)], zero_column, ptdfN[:,slack_position:end])
	# PTDF = transpose(PTDF)
	ptdfZ = ptdfN*gsk
	return ptdfN, ptdfZ
end

# Process input parameters
function process_parameters!(m::Model)
	m.ext[:parameters] = Dict()

	N = m.ext[:sets][:N]
    P = m.ext[:sets][:P]
    Z = m.ext[:sets][:Z]

	m.ext[:parameters][:gsk] = get_gsk_flat()
	m.ext[:parameters][:ptdfN], m.ext[:parameters][:ptdfZ] = build_ptdf(m)

	# # starten van nodale ptdf
	# # gsk als input
	# # zonale ptdfs berekenen
	# m.ext[:parameters][:ptdf_Z] = [0.5 0.25 0; 0.5 0.25 0; -0.5 -0.75 0; -0.5 0.25 0]
	# m.ext[:parameters][:ptdf_N] = [0.5 -0.25 0.25 0; 0.5 0.75 0.25 0; -0.5 -0.25 -0.75 0; -0.5 -0.25 0.25 0]

	#m.ext[:parameters][:f] = [333.33 133.33 -266.66 -200]
	#m.ext[:parameters][:f] = [350 150 -250 -250]
	m.ext[:parameters][:f] = [350 150 -250 -150]
	#m.ext[:parameters][:f] = [0 0 0 0]
	#m.ext[:parameters][:np] = [333.33 66.66 -400]
	#m.ext[:parameters][:np] = [350 0 -400]
	m.ext[:parameters][:np] = [300 100 -400]
	#m.ext[:parameters][:np] = [0 0 0]
	m.ext[:parameters][:p_in_n] = Dict(map(n -> n => [p for p in P if plant[plant[:,:GenID].==p, :OnBus][1] == n], N))
	m.ext[:parameters][:p_in_z] = Dict(map(z -> z => [p for p in P if plant[plant[:,:GenID].==p, :Zone][1] == z], Z))
	m.ext[:parameters][:n_in_z] = Dict(map(z -> z => [n for n in N if bus[bus[:,:BusID].==n, :Zone][1] == z], Z))

    return m
end

# call functions
define_sets!(m)
@info "Sets created"
process_parameters!(m)
@info "Parameters created"

## Step 2b: Additional Functions
# Get marginal cost of plant p
function get_mc(p)
	return plant[plant[:,:GenID].==p, :Costs][1]
end

# Get maximum capacity
function get_gen_up(p)
	return plant.Pmax[plant[:,:GenID].==p][1]
end

function get_dem(n)
	return load[1,n]
end

function get_line_cap(l)
	return branch.Pmax[branch[:,:BranchID].==l][1]
end

function find_maximum_mc()
	max_temp = 0
    for p in plant[:,:GenID]
        mc_temp = get_mc(p)
        if mc_temp > max_temp
            max_temp = mc_temp
        end
    end
	return max_temp
end

## Step 3: Construct your models

function model1!(m::Model)
	# Clear m.ext entries
	m.ext[:variables] = Dict()
	m.ext[:expressions] = Dict()
	m.ext[:constraints] = Dict()

	# Sets
	Z = m.ext[:sets][:Z]
	N = m.ext[:sets][:N]
	L = m.ext[:sets][:L]
	P = m.ext[:sets][:P]

	# Parameters
	ptdf_Z = m.ext[:parameters][:ptdf_Z]
	ptdf_N = m.ext[:parameters][:ptdf_N]
	np = m.ext[:parameters][:np]
	f = m.ext[:parameters][:f]
	p_in_n = m.ext[:parameters][:p_in_n]
	p_in_z = m.ext[:parameters][:p_in_z]

	# Variables
	GEN = m.ext[:variables][:GEN] = @variable(m, [p=P], lower_bound=0, upper_bound=get_gen_up(p), base_name="generation")
	NP = m.ext[:variables][:NP] = @variable(m, [z=Z], base_name="net position")

	# Expressions
	CG = m.ext[:expressions][:CG] = @expression(m,[z=Z],
        sum(GEN[p]*get_mc(p) for p in p_in_z[z])
    	)
	Fp = m.ext[:expressions][:Fp] = @expression(m, [l=L],
		sum(ptdf_N[l,n]*(GEN[n]-get_dem(n)) for n=N)
		)
	Fref = m.ext[:expressions][:Fref] = @expression(m, [l=L],
		f[l] - sum(ptdf_Z[l,z]*np[z] for z in Z)
		)
	RAMpos = m.ext[:expressions][:RAMpos] = @expression(m, [l=L],
		get_line_cap(l)-Fref[l]
		)
	RAMneg = m.ext[:expressions][:RAMneg] = @expression(m, [l=L],
		get_line_cap(l)+Fref[l]
		)
	Fc = m.ext[:expressions][:Fc] = @expression(m, [l=L],
		Fref[l]+sum(ptdf_Z[l,z]*NP[z] for z in Z)
		)

	# Objective
	m.ext[:objective] = @objective(m, Min,
		sum(CG[z] for z in Z)
		)

	# Constraints
	m.ext[:constraints][:con9f] = @constraint(m, [n=N],
		get_dem(n)
		==
		sum(GEN[p] for p in p_in_n[n])
		)

	@info "Model 1 executed"
	return m
end

function model2!(m::Model)
	# Model
	model1!(m)

	# Sets
	Z = m.ext[:sets][:Z]
	N = m.ext[:sets][:N]

	# Parameters
	n_in_z = m.ext[:parameters][:n_in_z]
	p_in_z = m.ext[:parameters][:p_in_z]

	# Variables
	GEN = m.ext[:variables][:GEN]

	# Constraints
	for n in N
		delete(m,m.ext[:constraints][:con9f][n])
	end

	m.ext[:constraints][:con9f] = @constraint(m, [z=Z],
		sum(get_dem(n) for n in n_in_z[z])
		==
		sum(GEN[p] for p in p_in_z[z])
		)

	@info "Model 2 executed"
	return m
end

function model3!(m::Model)
	# Model
	model2!(m)

	# Sets
	Z = m.ext[:sets][:Z]

	# Mapping
	n_in_z = m.ext[:parameters][:n_in_z]
	p_in_z = m.ext[:parameters][:p_in_z]

	# Variables
	GEN = m.ext[:variables][:GEN]
	NP = m.ext[:variables][:NP]

	# Constraints
	for z in Z
		delete(m,m.ext[:constraints][:con9f][z])
	end

	m.ext[:constraints][:con9f] = @constraint(m, [z=Z],
		sum(get_dem(n) for n in n_in_z[z])
		==
		sum(GEN[p] for p in p_in_z[z])
		- NP[z]
		)

	m.ext[:constraints][:con9i] = @constraint(m,
		sum(NP[z] for z in Z) == 0
		)

	@info "Model 3 executed"
	return m
end

function model4!(m::Model)
	# Model
	model3!(m)

	# Sets
	Z = m.ext[:sets][:Z]
	L = m.ext[:sets][:L]

	# Parameters
	ptdf_Z = m.ext[:parameters][:ptdf_Z]
	np = m.ext[:parameters][:np]
	f = m.ext[:parameters][:f]

	# Variables
	NP = m.ext[:variables][:NP]

	# Constraints
	m.ext[:constraints][:con9j] = @constraint(m,[l=L],
		get_line_cap(l) - f[l]
		>=
		sum(ptdf_Z[l,z]*(NP[z] - np[z]) for z in Z)
	 	)

	m.ext[:constraints][:con9k] = @constraint(m, [l=L],
		-get_line_cap(l) - f[l]
		<=
		sum(ptdf_Z[l,z]*(NP[z]-np[z]) for z in Z)
		)

	@info "Model 4 executed"
	return m
end

function model5!(n::Model)
	# Model
	n.ext[:variables] = Dict()
	n.ext[:expressions] = Dict()
	n.ext[:constraints] = Dict()

	# Sets
	N = n.ext[:sets][:N]
    L = n.ext[:sets][:L]
    P = n.ext[:sets][:P]
    Z = n.ext[:sets][:Z]

	# Parameters
	g = n.ext[:parameters][:g]
	n_in_z = n.ext[:parameters][:n_in_z]
	p_in_z = n.ext[:parameters][:p_in_z]
	p_in_n = n.ext[:parameters][:p_in_n]
	a=0.25
	ptdf_N = n.ext[:parameters][:ptdf_N]

	# Variables
	UP = n.ext[:variables][:UP] = @variable(n, [p=P], lower_bound=0, base_name="upward redispatch")
	DOWN = n.ext[:variables][:DOWN] = @variable(n, [p=P], lower_bound=0, base_name="downward redispatch")

	# Expressions
	CC = n.ext[:expressions][:CC] = @expression(n,[z=Z],
		sum((1+a)*UP[p]*get_mc(p)
		- (1-a)*DOWN[p]*get_mc(p) for p in p_in_z[z])
		)

	# Objective
	n.ext[:objective] = @objective(n, Min,
		sum(CC[z] for z in Z)
		)

	# Constraints
	n.ext[:constraints][:con15b] = @constraint(n,
		sum(UP[p] for p in P) == sum(DOWN[p] for p in P)
		)
	n.ext[:constraints][:con15c] = @constraint(n, [l=L],
		-get_line_cap(l)
		<=
		sum(ptdf_N[l,n]*
		(sum(g[p]+UP[p]-DOWN[p] for p in p_in_n[n]) - get_dem(n))
		for n in N)
		<=
		get_line_cap(l)
		)
	n.ext[:constraints][:con15d] = @constraint(n, [p=P],
		UP[p] <= get_gen_up(p) - g[p]
		)
	n.ext[:constraints][:con15e] = @constraint(n, [p=P],
		DOWN[p] <= g[p]
		)

	@info "Model 5 executed"
	return n
end

# function modelx!(m::Model)
	## Model
	# modelxmin1!()
	## Sets
	## Parameters
	## Variables
	## Expressions
	## Objective
	## Constraints
	#@info "Model x executed"
	#return m
# end

#model1!(m)
#model2!(m)
#model3!(m)
#model4!(m)

model4!(m)
optimize!(m)
g = value.(m.ext[:variables][:GEN])
# if redispatch nodig, doe dan redispatch
# anders niet
n = Model(optimizer_with_attributes(Gurobi.Optimizer))
define_sets!(n)
process_parameters!(n)
n.ext[:parameters][:g] = g

model5!(n)

## Step 4: Solve
optimize!(n)
@info "Model optimised"

println("Termination status: ", termination_status(n) )
println("Objective: ", value.(n.ext[:objective]))
println("")
print(n)
println("")

nodes = DataFrame(
	node = [n for n in m.ext[:sets][:N]],
	D = [get_dem(n) for n in m.ext[:sets][:N]],
	GMAX = [get_gen_up(p) for p in m.ext[:sets][:P]],
	MC = [get_mc(p) for p in m.ext[:sets][:P]],
	GEN = [value.(m.ext[:variables][:GEN][p]) for p in m.ext[:sets][:P]]
)

CSV.write("C:/Users/User/Documents/Thesis/Modeling/results/nodes.csv", nodes,  delim=';', decimal=',')

pretty_table(nodes, header=["node", "Dem[MW]", "Gmax[MW]","MC[€/MWh]","GEN[MW]"])

lines = DataFrame(
	line=["1-2","2-4","4-3","3-1"],
	cap=[get_line_cap(l) for l in m.ext[:sets][:L]],
	Fref=[value.(m.ext[:expressions][:Fref][l]) for l in m.ext[:sets][:L]],
	Fc=[value.(m.ext[:expressions][:Fc][l]) for l in m.ext[:sets][:L]],
	Fp=[value.(m.ext[:expressions][:Fp][l]) for l in m.ext[:sets][:L]],
	RAMpos=[value.(m.ext[:expressions][:RAMpos][l]) for l in m.ext[:sets][:L]],
	RAMneg=[value.(m.ext[:expressions][:RAMneg][l]) for l in m.ext[:sets][:L]]
)

CSV.write("C:/Users/User/Documents/Thesis/Modeling/results/lines.csv", lines,  delim=';', decimal=',')

pretty_table(lines, header=["line","cap[MW]","Fref[MW]","Fc[MW]","Fp[MW]","RAM+[MW]","RAM-[MW]"])

zones = DataFrame(
	zone=["A","B","C"],
	NP=[value.(m.ext[:variables][:NP][z]) for z in m.ext[:sets][:Z]],
	CG=[value.(m.ext[:expressions][:CG][z]) for z in m.ext[:sets][:Z]],
	#CC=[value.(m.ext[:expressions][:CC][z]) for z in m.ext[:sets][:Z]]
)

CSV.write("C:/Users/User/Documents/Thesis/Modeling/results/zones.csv", zones,  delim=';', decimal=',')

#pretty_table(zones, header=["zone", "NP[MW]", "CG[EUR]", "CC[EUR]"])


## Step 5: Interpretation
using Plots
