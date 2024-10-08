using Plots
using MAT

# Load the data from test.mat and noAttach.mat files
test_data = matread("MS1/test.mat")
no_attach_data = matread("MS1/noAttach.mat")

# Extract the bo_pa values from the loaded data
test_bo_pa1 = test_data["PATCH"]["bo_pa"][1, :]
test_bo_pa2 = test_data["PATCH"]["bo_pa"][2, :]
no_attach_bo_pa1 = no_attach_data["PATCH"]["bo_pa"][1, :]
no_attach_bo_pa2 = no_attach_data["PATCH"]["bo_pa"][2, :]

# Create the scatter plot
scatter(test_bo_pa1, no_attach_bo_pa1, xlabel = "test.mat bo_pa", 
ylabel = "noAttach.mat bo_pa", legend = false)
scatter(test_bo_pa2, no_attach_bo_pa2, xlabel = "test.mat bo_pa", 
ylabel = "noAttach.mat bo_pa", legend = false)
Plots.abline!(1, 0, line=:dash)

scatter(1:length(test_bo_pa1), test_bo_pa1, xlabel = "time step", 
ylabel = "bo_pa", legend = false, color = :red)
scatter!(1:length(no_attach_bo_pa1), no_attach_bo_pa1, legend = false,
color = :blue)

scatter(1:length(test_bo_pa2), test_bo_pa2, xlabel = "time step", 
ylabel = "bo_pa", legend = false, color = :red)
scatter!(1:length(no_attach_bo_pa2), no_attach_bo_pa2, legend = false,
color = :blue)
