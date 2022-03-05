## Flow Based Market Coupling
# author: Ben Callaerts
# last update: March 3, 2022

## Activate environment
cd("C:/Users/User/Documents/Thesis/Modeling")

using Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()
using CSV
using DataFrames
using YAML
using PrettyTables
using JuMP
using Gurobi
using Plots

##  Input data
data = "data-basic"

bus = CSV.read(joinpath(data, "bus.csv"), DataFrame; delim=";")
branch = CSV.read(joinpath(data, "branch.csv"), DataFrame; delim=";")
plant = CSV.read(joinpath(data, "plant.csv"), DataFrame; delim=";")
load = CSV.read(joinpath(data, "load.csv"), DataFrame; delim=";")
incidence = CSV.read(joinpath(data, "incidence.csv"), DataFrame; delim=";")
susceptance = CSV.read(joinpath(data, "susceptance.csv"), DataFrame; delim=";")
@info "Data read"

## Define functions

function create_model()
	m = Model(optimizer_with_attributes(Gurobi.Optimizer))

	m.ext[:sets] = Dict()
	get_sets!(m::Model)
	m.ext[:parameters] = Dict()
	get_parameters!(m::Model)
	m.ext[:expressions] = Dict()
	m.ext[:variables] = Dict()
	m.ext[:constraints] = Dict()

	return m
end

function get_sets!(m::Model)
    N = m.ext[:sets][:N] = bus[:,:BusID]
    L = m.ext[:sets][:L] = branch[:,:BranchID]
    P = m.ext[:sets][:P] = plant[:,:GenID]
    Z = m.ext[:sets][:Z] = sort(unique(bus[:,:Zone]))
end

function get_parameters!(m::Model)
	N = m.ext[:sets][:N]
    P = m.ext[:sets][:P]
    Z = m.ext[:sets][:Z]

	m.ext[:parameters][:gsk] = create_gsk_flat(m)
	m.ext[:parameters][:ptdfN], m.ext[:parameters][:ptdfZ] = build_ptdf(m)

	# # starten van nodale ptdf
	# # gsk als input
	# # zonale ptdfs berekenen
	# m.ext[:parameters][:ptdfZ] = [0.5 0.25 0; 0.5 0.25 0; -0.5 -0.75 0; -0.5 0.25 0]
	# m.ext[:parameters][:ptdfN] = [0.5 -0.25 0.25 0; 0.5 0.75 0.25 0; -0.5 -0.25 -0.75 0; -0.5 -0.25 0.25 0]

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

function build_ptdf!(m::Model)
	N = m.ext[:sets][:N]
	L = m.ext[:sets][:L]
	Z = m.ext[:sets][:Z]
	gsk = m.ext[:parameters][:gsk]

	MWBase = 380^2
	slack_node = 4
	slack_position = findfirst(N .== slack_node)

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

function create_gsk_flat!(m::Model)
	N = m.ext[:sets][:N]
	Z = m.ext[:sets][:Z]

	gsk_temp = zeros(Float64, length(N), length(Z))
	for n in N
		zone_temp = bus.Zone[bus[:,:BusID].==n][1]
		gsk_value_temp = 1/size(bus[bus[:,:Zone].==zone_temp,:],1)
		gsk_temp[findfirst(N .== n), findfirst(Z .== zone_temp)] = gsk_value_temp
	end
	return gsk_temp
end

function get_mc(p)
	return plant.Costs[plant[:,:GenID].==p][1]
end

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

function evaluate_congestion(m::Model)
	flow = abs.(value.(m.ext[:expressions][:Fp]))
	capacity = [get_line_cap(l) for l in m.ext[:sets][:L]]
	overflow = flow.>capacity
	return any(x->x==1, overflow)
end

function write_to_CSV(variable, heading, csv)
	data = zeros(Float64, 1, length(heading))
	data[:] = value.(variable[:])
	dataframe = DataFrame(data, :auto)
	rename!(dataframe, Dict(names(dataframe)[i] => Symbol.(heading[i]) for i = 1:ncol(dataframe)))
	CSV.write(joinpath("results", csv), dataframe,  delim=';', decimal=',')
end

## Market coupling
m = create_model()
@info "Market coupling model created"

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

function marketcoupling1!(m::Model)
	# Sets
	Z = m.ext[:sets][:Z]
	N = m.ext[:sets][:N]
	L = m.ext[:sets][:L]
	P = m.ext[:sets][:P]

	# Parameters
	ptdfN = m.ext[:parameters][:ptdfN]
	ptdfZ = m.ext[:parameters][:ptdfZ]
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
		sum(ptdfN[l,n]*(GEN[n]-get_dem(n)) for n=N)
		)
	Fref = m.ext[:expressions][:Fref] = @expression(m, [l=L],
		f[l] - sum(ptdfZ[l,z]*np[z] for z in Z)
		)
	RAMpos = m.ext[:expressions][:RAMpos] = @expression(m, [l=L],
		get_line_cap(l)-Fref[l]
		)
	RAMneg = m.ext[:expressions][:RAMneg] = @expression(m, [l=L],
		get_line_cap(l)+Fref[l]
		)
	Fc = m.ext[:expressions][:Fc] = @expression(m, [l=L],
		Fref[l]+sum(ptdfZ[l,z]*NP[z] for z in Z)
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

	@info "Market coupling 1 executed"
	return m
end

function marketcoupling2!(m::Model)
	# Model
	marketcoupling1!(m)

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

	@info "Market coupling 2 executed"
	return m
end

function marketcoupling3!(m::Model)
	# Model
	marketcoupling2!(m)

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

	@info "Market coupling 3 executed"
	return m
end

function marketcoupling4!(m::Model)
	# Model
	marketcoupling3!(m)

	# Sets
	Z = m.ext[:sets][:Z]
	L = m.ext[:sets][:L]

	# Parameters
	ptdfZ = m.ext[:parameters][:ptdfZ]
	np = m.ext[:parameters][:np]
	f = m.ext[:parameters][:f]

	# Variables
	NP = m.ext[:variables][:NP]

	# Constraints
	m.ext[:constraints][:con9j] = @constraint(m,[l=L],
		get_line_cap(l) - f[l]
		>=
		sum(ptdfZ[l,z]*(NP[z] - np[z]) for z in Z)
	 	)

	m.ext[:constraints][:con9k] = @constraint(m, [l=L],
		-get_line_cap(l) - f[l]
		<=
		sum(ptdfZ[l,z]*(NP[z]-np[z]) for z in Z)
		)

	@info "Market coupling 4 executed"
	return m
end

marketcoupling4!(m)
optimize!(m)
@info "Market coupling optimised"

## Congestion management
congestion = evaluate_congestion(m)

if congestion
	c = create_model()
	c.ext[:parameters][:g] = zeros(Float64, 1, length(m.ext[:sets][:P]))
	c.ext[:parameters][:g][:] = value.(m.ext[:variables][:GEN])[:]
	@info "Congestion management created"

	function congestion1!(c::Model)
		# Model
		c.ext[:variables] = Dict()
		c.ext[:expressions] = Dict()
		c.ext[:constraints] = Dict()

		# Sets
		N = c.ext[:sets][:N]
	    L = c.ext[:sets][:L]
	    P = c.ext[:sets][:P]
	    Z = c.ext[:sets][:Z]

		# Parameters
		g = c.ext[:parameters][:g]
		n_in_z = c.ext[:parameters][:n_in_z]
		p_in_z = c.ext[:parameters][:p_in_z]
		p_in_n = c.ext[:parameters][:p_in_n]
		a=0.25
		ptdfN = c.ext[:parameters][:ptdfN]

		# Variables
		UP = c.ext[:variables][:UP] = @variable(c, [p=P], lower_bound=0, base_name="upward redispatch")
		DOWN = c.ext[:variables][:DOWN] = @variable(c, [p=P], lower_bound=0, base_name="downward redispatch")

		# Expressions
		CC = c.ext[:expressions][:CC] = @expression(c,[z=Z],
			sum((1+a)*UP[p]*get_mc(p)
			- (1-a)*DOWN[p]*get_mc(p) for p in p_in_z[z])
			)

		# Objective
		c.ext[:objective] = @objective(c, Min,
			sum(CC[z] for z in Z)
			)

		# Constraints
		c.ext[:constraints][:con15b] = @constraint(c,
			sum(UP[p] for p in P) == sum(DOWN[p] for p in P)
			)
		c.ext[:constraints][:con15c] = @constraint(c, [l=L],
			-get_line_cap(l)
			<=
			sum(ptdfN[l,n]*
			(sum(g[p]+UP[p]-DOWN[p] for p in p_in_n[n]) - get_dem(n))
			for n in N)
			<=
			get_line_cap(l)
			)
		c.ext[:constraints][:con15d] = @constraint(c, [p=P],
			UP[p] <= get_gen_up(p) - g[p]
			)
		c.ext[:constraints][:con15e] = @constraint(c, [p=P],
			DOWN[p] <= g[p]
			)

		@info "Congestion management 1 executed"
		return c
	end

	congestion1!(c)
	optimize!(c)
	@info "Congestion management optimised"
else
	@info "No congestion management needed"
end

## Exporting to CSV
write_to_CSV(m.ext[:variables][:GEN], m.ext[:sets][:P], "g.csv")

## Post-processing
nodes = DataFrame(
	node = [n for n in m.ext[:sets][:N]],
	D = [get_dem(c) for n in m.ext[:sets][:N]],
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
