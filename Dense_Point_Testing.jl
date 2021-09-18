using Distributions
using LinearAlgebra
using Plots

gr()
# simulation function returns the trajectory if corresponding initial conditions
function simulation(coef, pair, mode="random")


    α, α1,δV,δC,δCV, β,γ1,ν1,K = coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7],coef[8],coef[9]


    δt = 0.001

    u = pair
    #name list is [V,C,Cᵥ,Cᵥ*]

    L = [-δV 0. 0. α1;
         0. (α-δC) 0. 0. ;
         0. 0. (-ν1-δCV) 0.;
         0. 0. ν1*β -β]

    f = collect(0. for i in 1:4)

    u_temp = zeros(4)
    cauchy = 1.
    i = 1
    check_lis = NaN
    trajectory = [u]

    while cauchy >0

        #if cauchy ≥ 100
        #    return [Inf for k in 1:8]
        #end


        i += 1

        f[2] = -α*u[2]*(sum(u[2:end]))/K - γ1 *u[2] *u[1]
        f[3] = γ1*u[2]*u[1]
        du = L * u + f

        u += δt * du

        if minimum(u) < 0
            println("Negative Value Detected")
            println(u)
            break
        end




        if mod(i, 50) ==0
            cauchy = norm(u - u_temp)
            push!(trajectory, u)
            if cauchy > 10000
                return "Diverges"
                break
            end
            u_temp = u
        end

        if cauchy == 0
            push!(trajectory,u)



    end
    return trajectory[end]
end


function coef_generator()
    Z = 100
    α, α1,δV,δC,δCV, β,γ1,ν1,K = 0. ,0.,0.,0.,0.,0.,0.,0.,0.
    while Z ≥ 1
        δC = rand(Uniform(0,1))

        α1 = rand(Uniform(0,1))
        δV = rand(Uniform(0,1))

        δCV = rand(Uniform(0,1))
        β = rand(Uniform(0,1))
        γ1 = rand(Uniform(0,1))
        ν1 = rand(Uniform(0,1))
        K = rand(Uniform(0,10))
        Z = δV*(δCV +ν1) / (γ1*α1*K*ν1)



        α  = (δC + δC / (1 - Z)) / 2
    end
    coef = [α, α1,δV,δC,δCV, β,γ1,ν1,K]
    println(coef)
    return coef

end


coef = [2.1800800128881375, 0.6513320017438788, 0.6586493297127192, 0.8700888066161852, 0.18583475923361203, 0.741477446883182, 0.5070311611578964, 0.6606320420811831, 6.367343311480813]
initial_list =[[0.0, 4.144234433030661, 2.996956052072762, 0.0], [0.0, 3.4214847054127606, 4.187433207855143, 0.0], [0.0, 3.011749153342076, 2.0705941839904094, 0.0], [0.0, 4.8058141751917205, 1.652041197057993, 0.0], [0.0, 0.88987672479378, 4.697002988212679, 0.0], [0.0, 4.21112650957401, 4.771169104112895, 0.0], [0.0, 1.399785355379639, 0.35963528828716673, 0.0]]
C_lis = []
CV_lis = []
V_lis = []
CVS_lis = []
α, α1,δV,δC,δCV, β,γ1,ν1,K = coef[1],coef[2],coef[3],coef[4],coef[5],coef[6],coef[7],coef[8],coef[9]

V_bar3 = (α*K*α1*ν1*γ1 - (α*(ν1+δCV)*δV +γ1*δC*K*α1*ν1))/(γ1*(δV*(1+ν1)*α+γ1*α1*ν1*K))
C_bar3 = δV*(δCV+ν1)/(α1*γ1*ν1)
CV_bar3 = V_bar3*δV/(α1*ν1)
CVS_bar3 = δV*V_bar3/(α1)

V_bar1 = 0.
C_bar1 = 0.
CV_bar1 = 0.
CVS_bar1 = 0.

V_bar2 = 0.
C_bar2 = K*(α-δC)/α
CV_bar2 = 0.
CVS_bar2 = 0.

for i in 1:10
    initial = initial_list[i]

    traj = simulation(coef, initial)
    push!(V_lis, collect(ele[1] for ele in traj))
    push!(C_lis,  collect(ele[2] for ele in traj))
    push!(CV_lis,  collect(ele[3] for ele in traj))
    push!(CVS_lis,  collect(ele[4] for ele in traj))
end
