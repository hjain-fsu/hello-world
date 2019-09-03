using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("DiffEqBase")
Pkg.add("DiffEqSensitivity")
Pkg.add("Test")
Pkg.add("OrdinaryDiffEq")
Pkg.add("ParameterizedFunctions")
Pkg.add("IterableTables")
Pkg.add("FileIO")
Pkg.add("Distributed")
Pkg.add("SharedArrays")
Pkg.add("Distributions")
Pkg.add("DataFrames")
Pkg.add("Tables")
Pkg.add("CSVFiles")
Pkg.add("Plots")

using Plots
using DataFrames
using Tables
using CSVFiles
using Distributions
using OnlineStats
using Distributed
using SharedArrays
using OrdinaryDiffEq, ParameterizedFunctions
using DifferentialEquations
using DiffEqBase
using DiffEqSensitivity
using Test
using IterableTables, FileIO
using Base.Threads


#######################################################################
#
# General model equations
#
#######################################################################

function DNA_damage_de(dx,x,p,t)

    ############################################
    #
    # DNA damage differential equations
    #
    # t: time
    # x: vector or variables
    # p: parameters
    # dx: differential equations
    #############################################


    ########################
    # Parameter values
    ########################


Vo = 2*10^(-3)
Vi = 7*10^(-6)
Kd = 0.0143;

k1t = 0.0440; #Gammaoi
k2t = 29.086; #Gammaio

k1 = 0.1270; #Damage + N7 -> C7
k2 = 0.0973; #Damage + N7 <- C7
k5 = 0.0130; # C7 -> A7

k3 = 0.5003; #Damage + O6 -> C6
k4 = 0.0457; #Damage + O6 <- C6
k6 = 0.0159; # C6 -> A6

########################
# Equations
########################


# TMZ Extracellular
dx[1] = -k1t*x[1]+k2t*(Vi/Vo)*x[2]-Kd*x[1]

# TMZ Intracellular
dx[2] = k1t*(Vo/Vi)*x[1] - k2t*x[2]-Kd*x[2]-k1*x[2]*x[3]-k3*x[2]*x[5]+k2*x[4]+k4*x[6]

# N7 sites on DNA
dx[3] = -k1*x[2]*x[3]+k2*x[4]

# DNA adducts at N7 sites [C7]
dx[4] = k1*x[2]*x[3]-k2*x[4]-k5*x[4]

# O6 sites on DNA
dx[5] = -k3*x[2]*x[5]+k4*x[6]

# DNA adducts at O6 sites [C6]
dx[6] = k3*x[2]*x[5]-k4*x[6]-k6*x[6]

dx[7] = k5*x[4]
dx[8] = k6*x[6]
end



function DNA_repair_de(dx,x,pars,t)

    ############################################
    #
    # DNA repair differential equations
    #
    # t: time
    # x: vector or variables
    # p: parameters
    # dx: differential equations
    #############################################

    ########################
    # Parameter values
    ########################

kfa = 0.025-0.001; # D7 + APNG -> C7
kra = 0.001; # D7 + APNG <- C7
kpa = 0.35; # C7 -> DNAr

kfm = 0.0916-0.001; # D6 + MGMT -> C6
krm = 0.001; # D6 + MGMT <- C6
kpm = 0.1717; # C6 -> DNAr


lambdaM = 0.0043;#0.0043
lambdaA = 0.00577;#0.1233

Sa = lambdaA*pars[1]
Sm = lambdaM*pars[2]

########################
# Equations
########################

# N7 sites on DNA
dx[1] = kpa*x[5]


# O6 sites on DNA
dx[2] = kpm*x[6]

# D7 [3]
dx[3] = -kfa*x[3]*x[8]+kra*x[5]

# D6 [4]
dx[4] = -kfm*x[4]*x[7]+krm*x[6]

# C7 [5]
dx[5] = kfa*x[3]*x[8]-kra*x[5]-kpa*x[5]

# C6 [6]
dx[6] = kfm*x[4]*x[7]-krm*x[6]-kpm*x[6]

# MGMT [7]
dx[7] = Sm - lambdaM*x[7]-kfm*x[4]*x[7]+krm*x[6]

# APNG [8]
dx[8] = Sa - lambdaA*x[8]-kfa*x[3]*x[8]+kra*x[5]

end


function DNA_damage_soln(tdata,dt,TMZ)

    ############################################
    #
    # DNA damage solution of differential equations
    #
    # tdata: time
    # dt: time step size
    # TMZ: drug dose
    #############################################

# Initial Conditions
x0 = [TMZ; 0.0; 13.95*TMZ; 0.0; 7.3196*TMZ; 0.0; 0.0; 0.0];

# Create ODE problem
Problem = ODEProblem(DNA_damage_de,x0,tdata)

# Solve ODE
F1 = solve(Problem,AutoTsit5(Rosenbrock23()),reltol=1e-8,abstol=1e-8,saveat=collect(tdata[1]:dt:tdata[2]));

# Report variables of interest
N = F1[3,:]+F1[4,:]+F1[5,:]+F1[6,:]+F1[7,:]+F1[8,:]
Y = [F1[7,:], F1[8,:], N]
return F1[7,:], F1[8,:], N
end


function DNA_repair_soln(pars,Ninit,Oinit,tdata,dt)

    ############################################
    #
    # DNA repair solution of differential equations
    #
    # tdata: time
    # dt: time step size
    # pars: parameters
    # Ninit: initial amount of N7 adducts
    # Oinit: initial amount of O6 adducts
    #############################################

# Initial conditions
x0 = [0.0; 0.0; Ninit;  Oinit; 0.0; 0.0; pars[1]; pars[2]];

# Create and solve ODE system
Problem = ODEProblem(DNA_repair_de,x0,tdata,pars)
F1 = solve(Problem,AutoTsit5(Rosenbrock23()),saveat=collect(tdata[1]:dt:tdata[2]));

return F1[1,:], F1[2,:], F1[7,:], F1[3,:]
end


function Population_model(p)

    ############################################
    #
    # Monoclonal in vitro simulations
    #
    # p: parameter values
    #############################################


########################
# Parameters varied
########################
alpha = p[1]
APNg = p[2]
MGMt = p[3]


########################
# Parameter values
########################

K = 5.8*10^6;

con = 500.0;
Na = 10
ha = 1/Na

ht = (1/3)*(ha)^2
Mt = 300#Int((1/ht))

t = (0.0,120.0)
age = zeros(1,Na)
a0 = age[Int(Na/2)];

for m=1:Na
    age[m] = (m-1.0)*ha;            #Time nodes position
end


mu1 = 174.2872;
rho1 = 2.9288;


b1 = 2.6411;# O6 death
b2 = 7.2196;# N7 death

c1 = 24.383;# Death  age
c2 = 55.0068;# O6 death exp MGMT
c4 = 0.0112;# N7 death exp A7


########################
# Equations
########################

# Age-Structured Variables

A = zeros(Na,Mt);     # Arrested Cells

# Only Time Dependent Variables

P = zeros(1,Mt);      # Proliferative cells
N = zeros(1,Mt);      # Total number of cells

# Initial Conditions
P[1] = 100000.0
N[1] = 100000.0
#println("break")
F2, F3, F4 = DNA_damage_soln(t,0.4)
mu = mu1*((F2+F3)./F4)

# Explicit recursion
for m=1:(Mt-1)
    totalD = F2[m]+F3[m]
    A7, A6, mg, D7 = DNA_repair_soln([APNg, MGMt],F2[m],F3[m],(0.0,9.0),1)
    if totalD==0.0
        repair=0.0*A7
    else
        repair = rho1*((A7+A6)./totalD)
    end

    P[m+1] = P[m]+ht*(alpha*P[m]*(1-(P[m]/(K)))-mu[m]*P[m]+sum(repair.*A[:,m])*ha)

    A[1,m+1] = mu[m]*P[m]


   for l=2:Na-1

       delta6 = b1*(1.0/(1+exp(-c1*(age[l]-a0))))*exp(-c2*mg[l]); #Death by
       delta7 = b2*(1.0/(1+exp(-c1*(age[l]-a0))))*(1.0/(1+exp(-c4*D7[l]))); #Death by N7

       repairA = repair[l]

       A[l,m+1]= A[l,m]+ht*(-(1.0/(2.0*ha))*(A[l+1,m]-A[l-1,m])-delta6*A[l,m]-delta7*A[l,m]-repairA*A[l,m])

   end
    A[end,m+1]=A[end-1,m+1]


    N[m+1] = P[m+1]+sum(A[:,m+1])*ha
    end
return P

end

#######################################################################
#
# Sobol Sensitivity Analysis
#
#######################################################################

################################################
# Create parameter vectors
################################################

function give_rand_p(p_range,p_fixed=nothing,indices=nothing)

################################################
# Creates parameter vectors
#
# p_range: range of parameter values
# p_fixed: parameters that will be fixed
# indices: indices of fixed parameters
################################################

    if p_fixed == nothing # If we are not fixing any parameters, create random vector within range
        p = [(p_range[j][2] - p_range[j][1])*rand() + p_range[j][1] for j in 1:length(p_range)]
    else
        p =  zeros(length(p_range))
        j = 1
        for i in 1:length(p_range)
            if i in indices
                p[i] = p_fixed[j]
                j += 1
            else
                p[i] = (p_range[i][2] -p_range[i][1])*rand() + p_range[i][1]
            end
        end
    end
    p
end



function calc_mean_var(f,p_range,N)

    ################################################
    # Calculate mean variance
    #
    # f: function
    # p_range: range of parameter values
    # N: sample size
    ################################################

    y1 = Array(f(give_rand_p(p_range)))
    y0 = zero(y1)
    v = zero(y1)
    for i in 1:N
        y1 = Array(f(give_rand_p(p_range)))
        @. y0 += y1
        @. v += y1^2
    end
    y0 = @. y0/N
    y0_sq = [i.^2 for i in y0]
    v = @. v/N - y0_sq
    y0,v
end

function first_order_var(f,p_range,N,y0)

    ################################################
    # Calculates first order variance
    #
    # f: function
    # p_range: range of parameter values
    # N: sample size
    # y0: initial mean value
    ################################################


    ys = Array{typeof(y0)}(undef,length(p_range))
    for i in 1:length(p_range)
        y = zero(y0)
        for j in 1:N
            p2 = give_rand_p(p_range)
            p1 = give_rand_p(p_range,[p2[i]],[i])
            yer =  Array(f(p1)) .* Array(f(p2))
            @. y += yer
        end
        y = @. y/N - y0^2
        ys[i] = copy(y)
    end
    ys
end

function second_order_var(f,p_range,N,y0)

    ################################################
    # Calculates second order variance
    #
    # f: function
    # p_range: range of parameter values
    # N: sample size
    # y0: initial mean value
    ################################################

    ys = Array{typeof(y0)}(undef,Int((length(p_range)*(length(p_range)-1))/2))
    curr = 1
    for i in 1:length(p_range)
        for j in i+1:length(p_range)
            y = zero(y0)
            for k in 1:N
                p2 = give_rand_p(p_range)
                p1 = give_rand_p(p_range,[p2[i],p2[j]],[i,j])
                y .+=  Array(f(p1)) .* Array(f(p2))
            end
            y = @. y/N - y0^2
            ys[curr] = copy(y)
            curr += 1
        end
    end
    ys_frst_order = first_order_var(f,p_range,N,y0)
    j = 1
    for i in 1:length(p_range)
        for k in i+1:length(p_range)
            ys[j] = @. ys[j] - ( ys_frst_order[i] + ys_frst_order[k] )
            j += 1
        end
    end
    ys
end


function total_var(f,p_range,N,y0)

    ################################################
    # Calculates total variance
    #
    # f: function
    # p_range: range of parameter values
    # N: sample size
    # y0: initial mean value
    ################################################

    ys = Array{typeof(y0)}(undef,length(p_range))
    for i in 1:length(p_range)
        y = zero(y0)
        for j in 1:N
            p_fixed_all = []
            p_fixed_indices = []
            p2 = give_rand_p(p_range)
            for j in 1:length(p2)
                if j != i
                    push!(p_fixed_all,p2[j])
                    push!(p_fixed_indices,j)
                end
            end
            p1 = give_rand_p(p_range,p_fixed_all,p_fixed_indices)
            yer =  Array(f(p1)) .* Array(f(p2))
            @. y += yer
        end
        y = @. y/N - y0^2
        ys[i] = copy(y)
    end
    ys
end

function sobol_sensitivity(f,p_range,N,order=2)

    ################################################
    # Calculates Sobol indices
    #
    # f: function
    # p_range: range of parameter values
    # N: sample size
    # order: desired Sobol index
    ################################################


    y0,v = calc_mean_var(f,p_range,N)
    if order == 1
        first_order = first_order_var(f,p_range,N,y0)
        for i in 1:length(first_order)
            first_order[i] = @. first_order[i] / v
        end
        first_order
    elseif order == 2
        second_order = second_order_var(f,p_range,N,y0)
        for i in 1:length(second_order)
            second_order[i] = @. second_order[i] / v
        end
        second_order
    else
        total_indices = total_var(f,p_range,N,y0)
        for i in 1:length(total_indices)
            total_indices[i] = @. 1 - (total_indices[i] / v)
        end
        total_indices
    end
end

################################################
# Sobol sensitivity in our case
################################################

p_range = [[1e-7,1.5],[1e-7,1500],[1e-7,1500],[1e-7,10]]
N = 1000
@time total = sobol_sensitivity(Population_model,p_range,N,0)
@time first_order = sobol_sensitivity(Population_model,p_range,N,1)
@time second_order = sobol_sensitivity(Population_model,p_range,N,2)






#######################################################################
#
# Virtual Cohorts
#
#######################################################################


################################################################
# Xenograft phenotypic evolution under conventional treatment
################################################################

function give_rand_p_energy(p_range)

    ################################################
    # Create an array of parameters within a range
    # that satisfy the energy equation
    #
    # p_range: range of parameter values
    ################################################

        p =  zeros(length(p_range))
        j = 0
        while j==0
            p = [(p_range[j][2] -p_range[j][1])*rand() + p_range[j][1] for j in 1:length(p_range)]
        if (p[1]*3.5+p[2]*1.58+p[3]*1.58)-2370.0>=0 # energy equation
               j=0        # If cell line does not follow energy equation continue
           else
               j=1
        end
    end
    p
end


function multiCell(ncl)
    ################################################
    # Create multiple cell lines randomly
    #
    # ncl: number of cell lines
    ################################################

    # Set range of parameters
    p_range = [[1e-7,1.55],[1e-7,1500],[1e-7,1500]]
    p = Array{Array{Float64, 4}}(undef, ncl)

    # Create parameters for each cell line
    p = [give_rand_p_energy(p_range)' for k=1:ncl]
    M = zeros(300,ncl)

    # Create a tumor with such cell lines under conventional treatment
    for i=1:ncl
        par = [p[i][1], p[i][2], p[i][3], 5]
        Y = Population_model(par)
        M[:,i] = Y'
    end
    return M,p

end


function AllCells(Nrun,ncl)

    ################################################
    # Create heteroclonal tumors and report their
    # phenotypic composition
    #
    # Nrun: number of tumors
    # ncl: number of cell lines per tumor
    ################################################

    Cell_Line = zeros(Nrun*ncl)
    Percentage = zeros(Nrun*ncl)
    Alpha = zeros(Nrun*ncl)
    APNG = zeros(Nrun*ncl)
    MGMT =zeros(Nrun*ncl)
    Tumor = zeros(Nrun*ncl)
    global cont
    cont = 1
    for k=1:Nrun # Create each tumor
        M, p = multiCell(ncl)
    for i=1:ncl # Save tumor information
        Cell_Line[cont]=i
        Percentage[cont]=M[end,i]/sum(M[end,:])
        Alpha[cont]=p[i][1]
        APNG[cont]=p[i][2]
        MGMT[cont]=p[i][3]
        Tumor[cont]=k
        cont+=1
    end
    end
    return Cell_Line, Percentage, Alpha, APNG, MGMT, Tumor
end

################################################
# Export data to a csv
################################################

Cell_Line, Percentage, Alpha, APNG, MGMT, Tumor = AllCells(500,10)
df = DataFrame(Cell_line=Cell_Line[:], Percentage=Percentage[:], Alpha=Alpha[:],
APNG=APNG[:], MGMT=MGMT[:], Tumor=Tumor[:])
out_file = open("Xenograft_composition.csv", "w")
df |> save("Xenograft_compostion.csv")

# Note: data set for Most_Least_resistant cells is obtained by outputting
# only cells in percentales 1 and 4 (respectively).


################################################
# Monoclonal treatment
################################################

function modelTx(p,tmin,tmax)

    ############################################
    #
    # Monoclonal in vitro simulations with
    # flexible treatment options
    #
    # p: parameter values
    # tmin: initial time
    # tmax: final time
    #############################################


########################
# Parameters varied
########################

alpha = p[1]
APNg = p[2]
MGMt = p[3]
TMZ = p[4]
P0 = p[5]
N0 = p[6]

########################
# Parameter values
########################

K = 5.8*10^6;

Na = 10
ha = 1/Na

Mt = 400
ht = (tmax-tmin)/400
t = (tmin,tmax)
age = range(0,stop=5,length=Na+1)
a0 = age[Int(Na/2)];

mu1 = 174.28;
rho1 = 2.9288;

b1 = 2.6411;# O6 death
b2 = 7.2196;# N7 death

c1 = 24.383;# Death  age
c2 = 55.0068;# O6 death exp MGMT
c4 = 0.0112;# N7 death exp A7


# Age-Structured Variables

A = zeros(Na+1,Mt+1);     # Arrested Cells

# Only Time Dependent Variables

P = zeros(1,Mt+1);      # Proliferative cells
N = zeros(1,Mt+1);      # Total number of cells


# Initial Conditions
P[1] = P0
N[1] = N0
#println("break")
F2, F3, F4 = DNA_damage_soln(t,ht,TMZ)
mu = mu1*((F2+F3)./F4)'
if TMZ == 0.0
    for m=1:(Mt)
        P[m+1] = P[m]+ht*(alpha*P[m]*(1-(N[m]/K)))
        N[m+1] = P[m+1]
    end
else
#Explicit recursion
for m=1:(Mt)
    totalD = F2[m]+F3[m]
    A7, A6, mg, D7 = DNA_repair_soln([APNg, MGMt],F2[m],F3[m],(0.0,5.0),5/Na)
    if totalD==0.0
        repair=0.0*A7
    else
        repair = rho1*((A7+A6)./totalD)
    end

    P[m+1] = P[m]+ht*(alpha*P[m]*(1-(N[m]/K))-mu[m]*P[m]+sum(repair.*A[:,m])*ha)

   A[1,m+1] = A[1,m]+ht*mu[m]*P[m]
   for l=2:Na

       delta6 = b1*(1.0/(1+exp(-c1*(age[l]-a0))))*exp(-c2*mg[l]); #Death by O6
       delta7 = b2*(1.0/(1+exp(-c1*(age[l]-a0))))*(1.0/(1+exp(-c4*D7[l]))); #Death by N7

       repairA = repair[l]

       A[l,m+1]= A[l,m]+ht*(-(1.0/(2.0*ha))*(A[l+1,m]-A[l-1,m])-delta6*A[l,m]-delta7*A[l,m]-repairA*A[l,m])

   end
    A[end,m+1]= 0
    N[m+1] = P[m+1]+sum(A[:,m+1])*ha
    end
end
return P, N

end



function OptimalTx(iAPNG,iMGMT,dose,tfinal,pInit)

    ############################################
    #
    # Treatment of in vitro simulations
    #
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # tfinal: total time of treatment
    # pInit: tumor parameter values
    #############################################


    alpha = pInit[1]
    APNG = pInit[2]
    MGMT = pInit[3]
    P0 = 10^5
    N0 = 10^5
    Txtime = 10.0

# Treatment cycles

day = 0.0*24.0

    pars = [alpha APNG*iAPNG[1] MGMT*iMGMT[1] dose[1] P0 N0]
    P, N = modelTx(pars,day,day+Txtime);
    Tumor = P;
    P0 = P[end];
    N0 = N[end];

    pars = [alpha APNG MGMT 0.0 P0 N0]
    P2, N2 = modelTx(pars,day+Txtime,day+24.0);
    Tumor = hcat(Tumor,P2);
    P0 = P2[end];
    N0 = N2[end];

day = 1.0*24.0
        pars = [alpha APNG*iAPNG[2] MGMT*iMGMT[2] dose[2] P0 N0]
        P, N = modelTx(pars,day,day+Txtime);
        Tumor = hcat(Tumor,P);
        P0 = P[end];
        N0 = N[end];

        pars = [alpha APNG MGMT 0.0 P0 N0]
        P2, N2 = modelTx(pars,day+Txtime,day+24.0);
        Tumor = hcat(Tumor,P2);
        P0 = P2[end];
        N0 = N2[end];

day = 2.0*24.0
        pars = [alpha APNG*iAPNG[3] MGMT*iMGMT[3] dose[3] P0 N0]
        P, N = modelTx(pars,day,day+Txtime);
        Tumor = hcat(Tumor,P);
        P0 = P[end];
        N0 = N[end];

        pars = [alpha APNG MGMT 0.0 P0 N0]
        P2, N2 = modelTx(pars,day+Txtime,day+24.0);
        Tumor = hcat(Tumor,P2);
        P0 = P2[end];
        N0 = N2[end];

day = 3.0*24.0
        pars = [alpha APNG*iAPNG[4] MGMT*iMGMT[4] dose[4] P0 N0]
        P, N = modelTx(pars,day,day+Txtime);
        Tumor = hcat(Tumor,P);
        P0 = P[end];
        N0 = N[end];

        pars = [alpha APNG MGMT 0.0 P0 N0]
        P2, N2 = modelTx(pars,day+Txtime,day+24.0);
        Tumor = hcat(Tumor,P2);
        P0 = P2[end];
        N0 = N2[end];

day = 4.0*24.0
        pars = [alpha APNG*iAPNG[5] MGMT*iMGMT[5] dose[5] P0 N0]
        P, N = modelTx(pars,day,day+Txtime);
        Tumor = hcat(Tumor,P);
        P0 = P[end];
        N0 = N[end];

# End of the cycle: no treatment until tfinal
        pars = [alpha APNG MGMT 0.0 P0 N0]
        P2, N2 = modelTx(pars,day+Txtime,day+48.0);
        Tumor = hcat(Tumor,P2);
        P0 = P2[end];
        N0 = N2[end];
        pars = [alpha APNG MGMT 0.0 P0 N0]
        P2, N2 = modelTx(pars,day+48.0,tfinal);
        Tumor = hcat(Tumor,P2);
        P0 = P2[end];
        N0 = N2[end];


return Tumor

end





#######################################################################
#
# Genetic Algorithm
#
#######################################################################

# Initial Population

function InitialPopulation(PopSize)

    ############################################
    #
    # Create Initial Treatment strategies
    #
    # Popsize: number of strategies
    #############################################

# Create Array
# Each population contains:
# iAPNG: boolean vector where a 1 in the ith position means APNG inhibitor is given in day i, 0 it is not given
# iMGMT: boolean vector where a 1 in the ith position means MGMT inhibitor is given in day i, 0 it is not given
# TMZ: real valued vector where dosages range from 0 to 500.
Population = Array{Array{Array{Float64,1},1},1}(undef, PopSize)

for m in 1:PopSize
        iAPNG =  rand([0.0,1.0],5)
        iMGMT =  rand([0.0,1.0],5)
        TMZ = [(500.0-0.0)*rand()+0.0 for j in 1:5]
        while sum(TMZ)>500 # check that maximum dosage in the strategy is safe
                TMZ = [(500.0-0.0)*rand()+0.0 for j in 1:5]
        end
        Population[m] = [iAPNG, iMGMT, TMZ]
end
return Population
end


function FitnessTx(Tx,pInit)

    ############################################
    #
    # Fitness function: calculates the fitness of
    # a treatment strategy based on tumor cell
    # population at the end of treatment
    #
    # Tx: treatmentstrategy
    # pInit: parameter values for cell lines
    #############################################


        Tumor = OptimalTx(Tx[1],Tx[1],Tx[3],10.0*24,pInit)
        Fitness = Tumor[end]

return Fitness
end



function SelectionTx(initPopulation,popsize,pInit)

    ###################################################
    #
    # Selection function: selects the most fit
    # treatment strategies
    #
    # initPopulation: strategies that we are studying
    # popsize: total number of strategies
    # pInit: parameter values for cell lines
    ####################################################

        selectionPer = 0.2; # we chose the top 20% strategies
        fitness = Array{Float64,1}(undef,popsize)
        for m in 1:popsize # calculate fitness for each strategy
                fitness[m]=FitnessTx(initPopulation[m],pInit)
        end
        # rank strategies based on fitness
        fittingPopulation = initPopulation[sortperm(fitness)]
        fitPopulation = fittingPopulation[1:Int(round(popsize*selectionPer))]
        return fitPopulation
end



function CrossoverTx(selectedPopulation,popsize)
    ########################################################
    #
    # Crossover function: creates new treatment strategies
    # from the top 20% by pairing them
    #
    # selectedPopulation: top 20% strategies by fitness
    # popsize: total number of strategies
    #########################################################
        childrenPopulation = Array{Array{Array{Float64,1},1},1}(undef, popsize)
        for m in 1:popsize
            # Select parents
                Parent1 = selectedPopulation[rand(Vector(1:size(selectedPopulation,1)))]
                Parent2 = selectedPopulation[rand(Vector(1:size(selectedPopulation,1)))]

            # Select at random from which parent each treatment value comes from
                apngDice = rand(Vector(1:5))
                mgmtDice = rand(Vector(1:5))
                tmzDice = rand(Vector(1:5))

            # Create child strategy
                childrenAPNG = Parent2[1][:]
                childrenAPNG[1:apngDice] = Parent1[1][1:apngDice]

                childrenMGMT = Parent2[2][:]
                childrenMGMT[1:mgmtDice] = Parent1[2][1:mgmtDice]

                childrenTMZ = Parent2[3][:]
                childrenTMZ[1:tmzDice] = Parent1[3][1:tmzDice]

                childrenPopulation[m] = [childrenAPNG, childrenMGMT, childrenTMZ]

        end
        return childrenPopulation
end




function MutationTx(newPopulation, popsize)
    ########################################################
    #
    # Mutation function: randomly mutates a small percentage
    # of treatment strategies
    #
    # newPopulation: new strategies (children)
    # popsize: total number of strategies
    #########################################################

        ratemutation = 0.05
        mutatedPopulation = newPopulation

         for m in 1:popsize # until we reach the desired population size
                mut = rand(3) # random vector with a value for iAPNG, iMGMT, and TMZ schedules
                for k in 1:3 # loop through each part of a treatment: iAPNG, iMGMT, and TMZ doses
                        if mut[k] <= ratemutation # if the random value is smaller than the mutation rate, we mutate
                                gene = rand(Vector(1:5))
                                # Calculate mutated values
                                mutValues = [rand([0.0,1.0]), rand([0.0,1.0]), (500.0-0.0)*rand()+0.0]

                                TMZcurrent = mutatedPopulation[m][3][gene]
                                mutatedPopulation[m][k][gene] = mutValues[k]

                                # In the case of TMZ we need to check for maximum dosage
                                if k == 3
                                        if sum(mutatedPopulation[m][k][:])>500 # if the new total TMZ dose is too high, repeat
                                                mutatedPopulation[m][k][gene] = TMZcurrent
                                        end
                                end
                        end
                end

        end

        return mutatedPopulation


end

function ChoosingTx(popsize, generations,pInit)

    ########################################################
    #
    # Genetic Algorithm to find optimal treatment strategies
    #
    # popsize: total number of strategies
    # generations: total number of generations allowed
    # pInit: parameter values of tumor cells
    #########################################################

    # 1. Create an initial set of treatment strategies
        currentPop = InitialPopulation(popsize)

        for l in 1:generations
            # 2. Select top 20%
                selectedPop = SelectionTx(currentPop,popsize,pInit)
            # 3. Create childrens from pairing the top strategies
                crossPop = CrossoverTx(selectedPop,popsize)
            # 4. Allow mutations
                mutantsPop = MutationTx(crossPop,popsize)
            # 5. Set the result as your current population
                currentPop = mutantsPop
        end

        finalSelection = SelectionTx(currentPop,popsize,pInit)

        # Final. Select the best strategy
        Treatment = finalSelection[1]

        return Treatment

end




################################################################
# Monoclonal treatment optimization
################################################################


function MonoclonalData(p_range,nT, con_dose)
    ########################################################
    #
    # Treatment Optimization monoclonal tumors
    #
    # p_range: range of parameter values
    #          - this range can be modified to create 4 cohorts
    #            such that we have: APNG+/MGMT+, APNG-/MGMT+
    #            APNG+/MGMT-, and APNG-/MGMT-
    # nT: total number of tumors
    # con_dose: conventional TMZ dose
    #########################################################
    
    # Create cell lines
    p = Array{Array{Float64, 4}}(undef, nT)
    p = [give_rand_p_energy(p_range)' for k=1:nT]
    local df1, df2

    iAPNG_on = [1	1	1	1	1]
    iMGMT_on = [1	1	1	1	1]

    for i in 1:nT
        pInit = p[i]
        # Find optimal treatment strategy for each tumor
        treatment = ChoosingTx(20,20,pInit)

        # Calculate tumor time courses under treatment
        tfinal = 10.0*24
        Tumor1 = OptimalTx(iAPNG_on,iMGMT_on,con_dose,tfinal,pInit)
        Tumor2 = OptimalTx(iAPNG_on,iMGMT_on,treatment[3],tfinal,pInit)
        Tumor3 = OptimalTx(treatment[1],treatment[2],treatment[3],tfinal,pInit)

        if i==1 # Record all the data
            df1 = DataFrame(Conventional_1=Tumor1[:], OnlyTMZ_1=Tumor2[:], Optimal_1=Tumor3[:])
            df2 = DataFrame(Alpha=pInit[1],
            APNG=pInit[2], MGMT=pInit[3], APNG_day1=treatment[1][1], APNG_day2=treatment[1][2],
            APNG_day3=treatment[1][3], APNG_day4=treatment[1][4],APNG_day5=treatment[1][5],
            MGMT_day1=treatment[2][1],MGMT_day2=treatment[2][2],MGMT_day3=treatment[2][2],
            MGMT_day4=treatment[2][4],MGMT_day5=treatment[2][5],
            TMZ_day1=treatment[3][1],TMZ_day2=treatment[3][2],TMZ_day3=treatment[3][3],
            TMZ_day4=treatment[3][4],TMZ_day5=treatment[3][5])
        else
            df1[:,Symbol("Conventional",'_',i)] = Tumor1[:]
            df1[:,Symbol("OnlyTMZ",'_',i)] = Tumor2[:]
            df1[:,Symbol("Optimal",'_',i)] = Tumor3[:]
            newdf2 = DataFrame(Alpha=pInit[1],
            APNG=pInit[2], MGMT=pInit[3], APNG_day1=treatment[1][1], APNG_day2=treatment[1][2],
            APNG_day3=treatment[1][3], APNG_day4=treatment[1][4],APNG_day5=treatment[1][5],
            MGMT_day1=treatment[2][1],MGMT_day2=treatment[2][2],MGMT_day3=treatment[2][2],
            MGMT_day4=treatment[2][4],MGMT_day5=treatment[2][5],
            TMZ_day1=treatment[3][1],TMZ_day2=treatment[3][2],TMZ_day3=treatment[3][3],
            TMZ_day4=treatment[3][4],TMZ_day5=treatment[3][5])
            df2 = append!(df2,newdf2)
        end

    end

out_file = open("OptimalTxCells.csv", "w")
df1 |> save("OptimalTxCells.csv")

out_file2 = open("OptimalTx.csv", "w")
df2 |> save("OptimalTx.csv")

end




################################################################
# Heteroclonal treatment comparison
################################################################


function OptimalTx_Multicell(iAPNG,iMGMT,dose,tfinal,pInit,ncl)

    ########################################################
    #
    # Extension of OptimalTx to simulate treatment of
    # heteroclonal xenografts
    #
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # tfinal: total time of treatment
    # pInit: tumor parameter values
    # ncl: number of cell lines per tumor
    #########################################################

    # First cell line
    Tumor = OptimalTx(iAPNG,iMGMT,dose,tfinal,pInit[1][:])
    TumorFinal = zeros(size(Tumor,2),ncl)
    Total = zeros(size(Tumor,2),1)
    TumorFinal[:,1] = Tumor'
    Total[:] =+ Tumor'

for i=2:ncl # all other cell lines
        Tumor = OptimalTx(iAPNG,iMGMT,dose,tfinal,pInit[i][:])
        TumorFinal[:,i] = Tumor'
        Total[:] = Total+Tumor'
end

return TumorFinal,Total

end


function multiCell_Tx(ncl,iAPNG,iMGMT,dose,tfinal)

    ########################################################
    #
    # Creates heteroclonal xenografts and simulates their
    # treatment under a strategy
    #
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # tfinal: total time of treatment
    # ncl: number of cell lines per tumor
    #########################################################

    # Create parameter values
    p_range = [[1e-7,1.55],[1e-7,1500],[1e-7,1500]]
    p = Array{Array{Float64, 4}}(undef, ncl)
    p = [give_rand_p_energy(p_range)' for k=1:ncl]

    # Simulation of treatment
    Y,Total = OptimalTx_Multicell(iAPNG,iMGMT,dose,tfinal,p,ncl)

    return Y,p,Total

end


#######################################################################
#
# Pre-clinical trial: Human simulations
#
#######################################################################


function human_modelTx(p,tmin,tmax)

    ############################################
    #
    # Human (in vivo) simulations
    #
    # p: parameter values
    # tmin: initial time
    # tmax: finish time
    #############################################

    ########################
    # Parameters varied
    ########################

alpha = p[1]/(70.0)
APNg = p[2]
MGMt = p[3]
TMZ = p[4]
P0 = p[5]
N0 = p[6]


########################
# Parameter values
########################

Na = 10
ha = 1/Na

Mt = 400#Int((1/ht))
ht = (tmax-tmin)/400
t = (tmin,tmax)
age = range(0,stop=5,length=Na+1)
a0 = age[Int(Na/2)];

# Parameters

mu1 = 174.28;
rho1 = 2.9288;

b1 = 2.6411;# O6 death
b2 = 7.2196;# N7 death

c1 = 24.383;# Death  age
c2 = 55.0068;# O6 death exp MGMT
c4 = 0.0112;# N7 death exp A7


# Age-Structured Variables

A = zeros(Na+1,Mt+1);     # Arrested Cells

# Only Time Dependent Variables

P = zeros(1,Mt+1);      # Proliferative cells
N = zeros(1,Mt+1);      # Total number of cells


# Initial Conditions
P[1] = P0
N[1] = N0
#println("break")
F2, F3, F4 = DNA_damage_soln(t,ht,TMZ)
mu = mu1*((F2+F3)./F4)'
if TMZ == 0.0
    for m=1:(Mt)
        P[m+1] = P[m]+ht*(alpha*P[m])
        N[m+1] = P[m+1]
    end
else
#Explicit recursion
for m=1:(Mt)
    totalD = F2[m]+F3[m]
    A7, A6, mg, D7 = DNA_repair_soln([APNg, MGMt],F2[m],F3[m],(0.0,5.0),5/Na)
    if totalD==0.0
        repair=0.0*A7
    else
        repair = rho1*((A7+A6)./totalD)
    end

    P[m+1] = P[m]+ht*(alpha*P[m]-mu[m]*P[m]+sum(repair.*A[:,m])*ha)

   A[1,m+1] = A[1,m]+ht*mu[m]*P[m]
   for l=2:Na

       delta6 = b1*(1.0/(1+exp(-c1*(age[l]-a0))))*exp(-c2*mg[l]); #Death by O6
       delta7 = b2*(1.0/(1+exp(-c1*(age[l]-a0))))*(1.0/(1+exp(-c4*D7[l]))); #Death by N7

       repairA = repair[l]

       A[l,m+1]= A[l,m]+ht*(-(1.0/(2.0*ha))*(A[l+1,m]-A[l-1,m])-delta6*A[l,m]-delta7*A[l,m]-repairA*A[l,m])

   end
    A[end,m+1]= 0
    N[m+1] = P[m+1]+sum(A[:,m+1])*ha

end
end
return P, N

end


function cycle_Tx_human(cycle,Palpha, APNG, MGMT, iAPNG, iMGMT, dose, P0, N0, Txtime)

    ########################################################################################
    #
    # Application of one treatment cycle to a patient
    #
    # cycle: cycle number
    # Palpha: proliferation rate
    # APNG: expression of APNG
    # MGMT: expression of MGMT
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # P0: initial number of proliferating cells
    # N0: initial total number of cells
    # Txtime: length of treatment
    #########################################################################################


day = cycle*28.0*24.0+0.0*24.0

          pars = [Palpha APNG*iAPNG[1] MGMT*iMGMT[1] dose[1] P0 N0]
          P, N = human_modelTx(pars,day,day+Txtime);
          Tumor = P;
          P0 = P[end];
          N0 = N[end];

          pars = [Palpha APNG MGMT 0.0 P0 N0]
          P2, N2 = human_modelTx(pars,day+Txtime,day+24.0);
          Tumor = hcat(Tumor,P2);
          P0 = P2[end];
          N0 = N2[end];

day = cycle*28.0*24.0+1.0*24.0
              pars = [Palpha APNG*iAPNG[2] MGMT*iMGMT[2] dose[2] P0 N0]
              P, N = human_modelTx(pars,day,day+Txtime);
              Tumor = hcat(Tumor,P);
              P0 = P[end];
              N0 = N[end];

              pars = [Palpha APNG MGMT 0.0 P0 N0]
              P2, N2 = human_modelTx(pars,day+Txtime,day+24.0);
              Tumor = hcat(Tumor,P2);
              P0 = P2[end];
              N0 = N2[end];

day = cycle*28.0*24.0+2.0*24.0
              pars = [Palpha APNG*iAPNG[3] MGMT*iMGMT[3] dose[3] P0 N0]
              P, N = human_modelTx(pars,day,day+Txtime);
              Tumor = hcat(Tumor,P);
              P0 = P[end];
              N0 = N[end];

              pars = [Palpha APNG MGMT 0.0 P0 N0]
              P2, N2 = human_modelTx(pars,day+Txtime,day+24.0);
              Tumor = hcat(Tumor,P2);
              P0 = P2[end];
              N0 = N2[end];

day = cycle*28.0*24.0+3.0*24.0
              pars = [Palpha APNG*iAPNG[4] MGMT*iMGMT[4] dose[4] P0 N0]
              P, N = human_modelTx(pars,day,day+Txtime);
              Tumor = hcat(Tumor,P);
              P0 = P[end];
              N0 = N[end];

              pars = [Palpha APNG MGMT 0.0 P0 N0]
              P2, N2 = human_modelTx(pars,day+Txtime,day+24.0);
              Tumor = hcat(Tumor,P2);
              P0 = P2[end];
              N0 = N2[end];

day = cycle*28.0*24.0+4.0*24.0
              pars = [Palpha APNG*iAPNG[5] MGMT*iMGMT[5] dose[5] P0 N0]
              P, N = human_modelTx(pars,day,day+Txtime);
              Tumor = hcat(Tumor,P);
              P0 = P[end];
              N0 = N[end];

              pars = [Palpha APNG MGMT 0.0 P0 N0]
              P2, N2 = human_modelTx(pars,day+Txtime,day+48.0);
              Tumor = hcat(Tumor,P2);
              P0 = P2[end];
              N0 = N2[end];

      # From day 7 until 28
              pars = [Palpha APNG MGMT 0.0 P0 N0]
              P2, N2 = human_modelTx(pars,day+48.0,(cycle+1)*28.0*24.0);
              Tumor = hcat(Tumor,P2);
              P0 = P2[end];
              N0 = N2[end];
return Tumor, P0, N0
end




function human_TX_cycle(iAPNG,iMGMT,dose,months,noTxmonths,pInit)
    ##############################################################
    #
    # Human (in vivo) simulation of multiple treatment cycles
    #
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # months: total number of months under treatment
    # noTxmonths: total number of months without treatment
    # pInit: parameter values
    #############################################################

    ########################
    # Parameter values
    ########################

    alpha = pInit[1]
    APNG = pInit[2]
    MGMT = pInit[3]
    P0 = 5*10^7
    N0 = 5*10^7
    Txtime = 10.0


    ########################
    # Treatment
    ########################


    Tumor1, P01, N01, = cycle_Tx_human(0,alpha,APNG,MGMT,iAPNG,iMGMT,dose,P0,N0,Txtime)
    Tumor = Tumor1

    for i = 1:months # treat for several months
        Tumor2, P02, N02, = cycle_Tx_human(i,alpha,APNG,MGMT,iAPNG,iMGMT,dose,P01,N01,Txtime)
        Tumor = hcat(Tumor,Tumor2)
        P01 = P02
        N01 = N02
    end


    #################################
    # Study growth after treatment
    #################################

    if noTxmonths > 0

        for k = 1:noTxmonths
        TumorLast, P02, N02, = cycle_Tx_human(k,alpha,APNG,MGMT,iAPNG,iMGMT,0.0*dose,P01,N01,Txtime)
        Tumor = hcat(Tumor,TumorLast)
        P01 = P02
        N01 = N02
        end

    end

    return Tumor

end



function human_OptimalTx_Multicell(iAPNG,iMGMT,dose,tfinal,noTxmonths,pInit,ncl)

    ########################################################
    #
    # Simulation of heteroclonal human tumor treatment
    #
    #
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # tfinal: total number of months under treatment
    # noTxmonths: total number of months without treatment
    # pInit: parameter values
    # ncl: number of cell lines per tumor
    #######################################################


    # First cell lines
    Tumor = human_TX_cycle(iAPNG,iMGMT,dose,tfinal,noTxmonths,pInit[1][:])
    TumorFinal = zeros(size(Tumor,2),ncl)
    Total = zeros(size(Tumor,2),1)
    TumorFinal[:,1] = Tumor'
    Total[:] =+ Tumor'
    for i=2:ncl # All other cell lines
        Tumor = human_TX_cycle(iAPNG,iMGMT,dose,tfinal,noTxmonths,pInit[i][:])
        TumorFinal[:,i] = Tumor'
        Total[:] = Total+Tumor'
    end

    return TumorFinal,Total

end

function human_Outcome(ncl,nT,tfinal,noTxmonths, con_dose, opt_dose)

    ############################################
    #
    # Export data of human simulations
    #
    # ncl: number of cell lines per tumor
    # nT: number of patients
    # tfinal: number of months under treatment
    # noTxmonths: total number of months of study
    # con_dose: TMZ conventional dose schedule
    # opt_dose: TMZ optimal dose schedule
    #############################################

# Decision about inhibitors on/off
        iAPNG_on = [1	1	1	1	1]
        iMGMT_on = [1	1	1	1	1]
        iAPNG_off = [0	0	0	0	0]
        iMGMT_off = [0	0	0	0	0]

local df1a, df1b, df1c
for i in 1:nT
        p_range = [[1e-7,1.55],[1e-7,1500],[1e-7,1500]]
        p = Array{Array{Float64, 4}}(undef, ncl)
        p = [rand_p_energy(p_range)' for k=1:ncl]

        Y1,Total1 = human_OptimalTx_Multicell(iAPNG_on,iMGMT_on,con_dose,tfinal,noTxmonths,p,ncl)
        Y2,Total2 = human_OptimalTx_Multicell(iAPNG_on,iMGMT_on,opt_dose,tfinal,noTxmonths,p,ncl)
        Y3,Total3 = human_OptimalTx_Multicell(iAPNG_off,iMGMT_ff,opt_dose,tfinal,noTxmonths,p,ncl)
if i==1
df1a = DataFrame(Conventional_1=Total1[:])
df1b = DataFrame(Optimal_TMZ_1=Total2[:])
df1c = DataFrame(Treatment_D_1=Total3[:])

else
      df1a[:,Symbol("Conventional",'_',i)] = Total1[:]
      df1b[:,Symbol("Optimal_TMZ",'_',i)] = Total2[:]
      df1c[:,Symbol("Treatment_D",'_',i)] = Total3[:]
end

end

return df1a, df1b, df1c

end



function human_info_cells(iAPNG,iMGMT,dose,tfinal,noTxmonths,pInit,ncl)
    ################################################
    #
    # Export phenotype information in human trial
    #
    # iAPNG: vector with treatment decisions about APNG inhibitor (0 or 1 per day)
    # iMGMT: vector with treatment decisions about MGMT inhibitor (0 or 1 per day)
    # dose: vector with TMZ doses
    # tfinal: total number of months under treatment
    # noTxmonths: total number of months without treatment
    # pInit: parameter values
    # ncl: number of cell lines per tumor
    #################################################

    # Initialize variables
    Tumor = human_Tx_cycle(iAPNG,iMGMT,dose,tfinal,noTxmonths,pInit[1][:])
    TumorFinal = zeros(size(Tumor,2),ncl)
    Total = zeros(size(Tumor,2),1)

    Alpha = zeros(size(Tumor,2)+1,1)
    APNG = zeros(size(Tumor,2)+1,1)
    MGMT = zeros(size(Tumor,2)+1,1)

    All_alpha = zeros(ncl,1)
    All_APNG = zeros(ncl,1)
    All_MGMT = zeros(ncl,1)

    TumorFinal[:,1] = Tumor'
    Total[:] =+ Tumor'


for i=2:ncl # create each cell
        Tumor = human_Tx_cycle(iAPNG,iMGMT,dose,tfinal,noTxmonths,pInit[i][:])
        TumorFinal[:,i] = Tumor'
        Total[:] = Total+Tumor'
end

# Obtain information about cell line with highest proliferation rate,
# cell line with highest APNG expression, and cell line with highest
# MGMT expression

for k=1:ncl
    All_alpha[k] = pInit[k][1]
    All_APNG[k]  = pInit[k][2]
    All_MGMT[k]  = pInit[k][3]
end

fastAlpha = maximum(All_alpha)
positionFast = [i for (i, x) in enumerate(All_alpha) if x == fastAlpha]
fastAPNG = p[positionFast[1]][2]
fastMGMT = p[positionFast[1]][3]

maxAPNG = maximum(All_APNG)
maxMGMT = maximum(All_MGMT)

InfoCells = [fastAlpha, fastAPNG, fastMGMT, maxAPNG, maxMGMT]


# Calculate weighted average value of each parameter
for m=1:ncl
    Alpha[1] = Alpha[1] + pInit[m][1]/ncl
    APNG[1]  = APNG[1] + pInit[m][2]/ncl
    MGMT[1]  = MGMT[1] + pInit[m][3]/ncl
end

for j=2:size(Tumor,2)+1

        for m=1:ncl
            Alpha[j] = Alpha[j] + pInit[m][1]*(TumorFinal[j-1,m])/(Total[j-1])
            APNG[j]  = APNG[j] + pInit[m][2]*(TumorFinal[j-1,m])/(Total[j-1])
            MGMT[j]  = MGMT[j] + pInit[m][3]*(TumorFinal[j-1,m])/(Total[j-1])
        end

end
return TumorFinal,Total,Alpha,APNG,MGMT, InfoCells

end


function human_cells_Outcome(ncl,nT,tfinal,noTxmonths,con_dose, opt_dose)

    ################################################################
    #
    # Export data of phenotypic composition in human simulations
    #
    # ncl: number of cell lines per tumor
    # nT: number of patients
    # tfinal: number of months under treatment
    # noTxmonths: total number of months of study
    # con_dose: TMZ conventional dose schedule
    # opt_dose: TMZ optimal dose schedule
    ##############################################################

# Decision about inhibitors on/off
        iAPNG_on = [1	1	1	1	1]
        iMGMT_on = [1	1	1	1	1]
        iAPNG_off = [0	0	0	0	0]
        iMGMT_off = [0	0	0	0	0]

local df1a, df1b, df1c, df2a, df2b, df2c, df3a, df3b, df3c, df4a, df4b, df4c

# Create the date
for i in 1:nT
        p_range = [[1e-7,1.55],[1e-7,1500],[1e-7,1500]]
        p = Array{Array{Float64, 4}}(undef, ncl)
        p = [rand_p_energy(p_range)' for k=1:ncl]

        Y1,Total1, Alpha1, APNG1, MGMT1, InfoCells1 = human_info_cells(iAPNG_on,iMGMT_on,con_dose,tfinal,noTxmonths,p,ncl)
        Y2,Total2, Alpha2, APNG2, MGMT2, InfoCells2 = human_info_cells(iAPNG_on,iMGMT_on,opt_dose,tfinal,noTxmonths,p,ncl)
        Y3,Total3, Alpha3, APNG3, MGMT3, InfoCells3 = human_info_cells(iAPNG_off,iMGMT_off,opt_dose,tfinal,noTxmonths,p,ncl)

# Export the data
if i==1
df1a = DataFrame(Conventional_1=Alpha1[:])
df1b = DataFrame(Optimal_TMZ_1=Alpha2[:])
df1c = DataFrame(Treatment_D_1=Alpha3[:])


df2a = DataFrame(Conventional_1=APNG1[:])
df2b = DataFrame(Optimal_TMZ_1=APNG2[:])
df2c = DataFrame(Treatment_D_1=APNG3[:])


df3a = DataFrame(Conventional_1=MGMT1[:])
df3b = DataFrame(Optimal_TMZ_1=MGMT2[:])
df3c = DataFrame(Treatment_D_1=MGMT3[:])

df4a = DataFrame(Conventional_1=InfoCells1[:])
df4b = DataFrame(Optimal_TMZ_1=InfoCells2[:])
df4c = DataFrame(Treatment_D_1=InfoCells3[:])

else
      df1a[:,Symbol("Conventional",'_',i)] = Alpha1[:]
      df1b[:,Symbol("Optimal_TMZ",'_',i)] = Alpha2[:]
      df1c[:,Symbol("Treatment_D",'_',i)] = Alpha3[:]

      df2a[:,Symbol("Conventional",'_',i)] = APNG1[:]
      df2b[:,Symbol("Optimal_TMZ",'_',i)] = APNG2[:]
      df2c[:,Symbol("Treatment_D",'_',i)] = APNG3[:]


      df3a[:,Symbol("Conventional",'_',i)] = MGMT1[:]
      df3b[:,Symbol("Optimal_TMZ",'_',i)] = MGMT2[:]
      df3c[:,Symbol("Treatment_D",'_',i)] = MGMT3[:]

      df4a[:,Symbol("Conventional",'_',i)] = InfoCells1[:]
      df4b[:,Symbol("Optimal_TMZ",'_',i)] = InfoCells2[:]
      df4c[:,Symbol("Treatment_D",'_',i)] = InfoCells3[:]
end

end

return df1a, df1b, df1c, df2a, df2b, df2c, df3a, df3b, df3c, df4a, df4b, df4c

end
