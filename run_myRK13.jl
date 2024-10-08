#=---------------------------------------
 	CARBON DYNAMICS MODEL PER PARTICLE  
---------------------------------------
 Created by Trang Nguyen
 Levine Lab, Univ of Southern California =#
 using MAT
 using DifferentialEquations


 # Initialize model
include("ModelInitialization.jl")
include("myRK13.jl")
include("createVector.jl")
include("createStructure.jl")
include("update_r_values_myRK13.jl")
include("update_depth_params_myRK13.jl")


# Main model loop
#i=1;
#count=0;

#PATCH["time"] = [convert(Float64, PATCH["time"])]
#PATCH["pom"] = [convert(Float64, PATCH["pom"])]
#PATCH["mom"] = [convert(Float64, PATCH["mom"])]
#PATCH["z"] = [convert(Float64, PATCH["z"])]

function run_model(PATCH, timestep, ndays)
    i = 1
    count = 0
    try
        while PATCH["time"][i]<ndays || i<timestep
            C=createVector(PATCH, i);
            if i % 10 == 0
                println("Step $i completed")
            end
            u0 = C
            tspan = (PATCH["time"][i], PATCH["time"][i] + dt)

            prob = ODEProblem(myRK13, u0, tspan)
            sol = solve(prob, Tsit5()) #Tsit5 is similar to ode45

            mytime = sol.t
            y = sol.u

            #[mytime,y]=ode45(@myRK13,[0 dt], C); 
            #more efficient to preallocate
            PATCH["time"][i+1]=PATCH["time"][i]+dt;
            iAttach= newAttachedCells;
            PATCH=createStructure(PATCH, y[end,:], i+1);
            
        # Time i
            w_aggregate=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(PATCH["r"][i]/100).^(D_Omand-1)*rs_3_minus_d/(18*mu)*100;#cm/s, Omand 2020, POM Note4 Trang, Sinking rate based on density at certain depth and based on the change of fluid density at depth
            w_smallpa=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(PATCH["r"][i]/100).^(D_s-1)/(18*mu)*100;#cm/s, Sinking rate for solid particle fraction dimension D=3 aka a_s^(D-3)=1
            if PATCH["r"][i]<=cut_off
                Ws_cms=w_smallpa; #cm/s
            else
                Ws_cms=w_aggregate;#cm/s
            end
            Ws=Ws_cms*864; #m/day
            
            rey= Ws_cms * PATCH["r"][i] / viscosity; #unitless
            PATCH["Sh"][i]=1+0.619*rey^0.412*sc^(1/3);#Sherwood number
            
            # Calculate new RADIUS based on new carbon content of the particle
            new_c_avg_pa=PATCH["pom"][i]; #uM.particle-1 
            
            f_mult_vol=new_c_avg_pa/convert_Alldredge_Omand/(c_avg_pa_ref/12*1000);# f * vol
            pseudo_r1=(3/4*f_mult_vol/pi)^(1/3); #r calculated if f=1 aka small particles
            pseudo_r2=(f_mult_vol/(4/3*pi*Af*rs_3_minus_d*100^(3-D_Omand)))^(1/D_Omand); #r calcualted if f<1
            if pseudo_r1<= cut_off
                PATCH["r"][i+1]=pseudo_r1;
            else
                PATCH["r"][i+1]=pseudo_r2;
            end
    
            # Calculate sinking depth
            w_aggregate=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(PATCH["r"][i]/100).^(D_Omand-1)*rs_3_minus_d/(18*mu)*100;#cm/s, Omand 2020, POM Note4 Trang, Sinking rate based on density at certain depth and based on the change of fluid density at depth
            w_smallpa=g*Af*(c_avg_pa_ref-fluid_density_ref)*1000*(PATCH["r"][i]/100).^(D_s-1)/(18*mu)*100;#cm/s, Sinking rate for solid particle fraction dimension D=3 aka a_s^(D-3)=1
            if PATCH["r"][i]<=cut_off
                Ws_cms=w_smallpa; #cm/s
            else
                Ws_cms=w_aggregate;#cm/s
            end
            Ws=Ws_cms*864; #m/day
            
            PATCH["depth"][i+1]=PATCH["depth"][i] + Ws*dt; #depth is in m


            i = i + 1
        end

        return PATCH
    catch e
        println("Error at time $i")
        println(e)
        return PATCH
    end
end

PATCH = run_model(PATCH, timestep, ndays)

# SAVE OUTPUT to .mat files
filename = "MS1/" * runName * ".mat"
file = matopen(filename, "w")
write(file, "PATCH", PATCH)
close(file)
