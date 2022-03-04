## Backbone to structure your code
# author: Kenneth Bruninx, Sebastian Gonzato
# last update: October 26, 2020
# description: backbone to structure your code. You're not obliged to use this
# in your assignment, but you may.

## Step 0: Activate environment - ensure consistency accross computers
using Pkg
Pkg.activate(@__DIR__) # @__DIR__ = directory this script is in
Pkg.instantiate() # If a Manifest.toml file exist in the current project, download all the packages declared in that manifest. Else, resolve a set of feasible packages from the Project.toml files and install them.

##  Step 1: input data
using CSV
using DataFrames
using YAML

data = YAML.load_file(joinpath(@__DIR__, "data_gep.yaml"))
ts = CSV.read(joinpath(@__DIR__, "Profiles_12_reprdays.csv"), DataFrame)
repr_days = CSV.read(joinpath(@__DIR__, "Weights_12_reprdays.csv"), DataFrame)

## Step 2: create model & pass data to model
using JuMP
using Gurobi
m = Model(optimizer_with_attributes(Gurobi.Optimizer))

# Step 2a: create sets
function define_sets!(m::Model, data::Dict)
    # create dictionary to store sets
    m.ext[:sets] = Dict()

    # define the sets
    m.ext[:sets][:JH] = 1:data["nTimesteps"] # Timesteps
    m.ext[:sets][:JD] = 1:data["nReprDays"] # Representative days
    m.ext[:sets][:ID] = [id for id in keys(data["dispatchableGenerators"])] # dispatchable generators
    m.ext[:sets][:IV] = [iv for iv in keys(data["variableGenerators"])] # variable generators
    m.ext[:sets][:I] = union(m.ext[:sets][:ID], m.ext[:sets][:IV]) # all generators

    # return model
    return m
end

# Step 2b: add time series
function process_time_series_data!(m::Model, data::Dict, ts::DataFrame)
    # extract the relevant sets
    JH = m.ext[:sets][:JH] # Time steps
    JD = m.ext[:sets][:JD] # Days

    # create dictionary to store time series
    m.ext[:timeseries] = Dict()

    # example: add time series to dictionary
    m.ext[:timeseries][:D] = [ts.Load[jh+data["nTimesteps"]*(jd-1)] for jh in JH, jd in JD]

    # return model
    return m
end

# step 2c: process input parameters
function process_parameters!(m::Model, data::Dict, repr_days::DataFrame)
    # extract the sets you need
    I = m.ext[:sets][:I]

    # generate a dictonary "parameters"
    m.ext[:parameters] = Dict()

    # example: legacy capacity
    d = merge(data["dispatchableGenerators"],data["variableGenerators"])
    LC = m.ext[:parameters][:LC] = Dict(i => d[i]["legcap"] for i in I) # MW
    VC = m.ext[:parameters][:VC] = Dict(id => d[id]["fuelCosts"] for id in ID) # MW

    # return model
    return m
end

# call functions
define_sets!(m, data)
process_time_series_data!(m, data, ts)
process_parameters!(m, data, repr_days)

## Step 3: construct your model
# Greenfield GEP - single year (Lecture 3 - slide 25, but based on representative days instead of full year)
function build_GEP_model!(m::Model)
    # Clear m.ext entries "variables", "expressions" and "constraints"
    m.ext[:variables] = Dict()
    m.ext[:expressions] = Dict()
    m.ext[:constraints] = Dict()

    # Extract sets
    JH = m.ext[:sets][:JH]
    JD = m.ext[:sets][:JD]
    ID = m.ext[:sets][:ID]
    I = m.ext[:sets][:I]
    IV = m.ext[:sets][:IV]

    # Extract time series data
    D = m.ext[:timeseries][:D] # demand

    # Extract parameters
    LC = m.ext[:parameters][:LC]
    VC = m.ext[:parameters][:VC]

    # Create variables
    cap = m.ext[:variables][:cap] = @variable(m, [i=I], lower_bound=0, base_name="capacity")
    g = m.ext[:variables][:g] = @variable(m, [i=I,jh=JH,jd=JD], lower_bound=0, base_name="generation")

    # Create affine expressions (= linear combinations of variables)
    # dummy example
    dummy = m.ext[:expressions][:dummy] = @expression(m, [i=IV,jh=JH,jd=JD],
        cap[i] - g[i,jh,jd]
    )

    # Formulate objective 1a
    m.ext[:objective] = @objective(m, Min,
        sum(VC[i]*g[i,jh,jd] for i in ID, jh in JH, jd in JD)
    )

    # constraints
    # 3a1 - conventional
    m.ext[:constraints][:con3a1conv] = @constraint(m, [i=ID,jh=JH,jd=JD],
        g[i,jh,jd] <= (cap[i]+LC[i])
    )

    return m
end

# recall that you can, once you've built a model, delete and overwrite constraints using the appropriate reference:
# example:   delete(m,m.ext[:constraints][:con3a1conv][id,jh,jd]) (needs to be done in for-loop over id, jh, jd)

# Build your model
build_GEP_model!(m)

## Step 4: solve
# current model is incomplete, so all variables and objective will be zero
optimize!(m)

# check termination status
print(
    """

    Termination status: $(termination_status(m))

    """
)


# print some output
@show value(m.ext[:objective])
@show value.(m.ext[:variables][:cap])

## Step 5: interpretation
using Plots
using StatsPlots

# examples on how to access data

# sets
JH = m.ext[:sets][:JH]
JD = m.ext[:sets][:JD]
I = m.ext[:sets][:I]

# parameters
D = m.ext[:timeseries][:D]
LC = m.ext[:parameters][:LC]

# variables/expressions
g = value.(m.ext[:variables][:g])
# λ = dual.(m.ext[:constraints][:con2a]) # con2a is not defined, will not work

# create arrays for plotting
# λvec = [λ[jh,jd] for jh in JH, jd in JD]
capvec = [cap[i] for  i in I]
