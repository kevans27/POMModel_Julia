using MAT

file = matopen("test.mat", "r")
PATCH = read(file, "PATCH")
close(file)