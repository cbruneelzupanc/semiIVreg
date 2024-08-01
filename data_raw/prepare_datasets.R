N = 100000; set.seed(1234)

# Specification
model_type = "heterogenous"
param_error = c(1, 1, 0.6, 0.5)
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9)
param_p = c(0, -0.7, 0.7, 0, 0, 0)
param_y0 = c(3.2, 0.8, 0, 0)
param_y1 = c(3.2+0.4, 0.5, 0, 0)
param_genX = c(0.4, 0, 2)

roydata = simul_data(N, model_type, param_y0, param_y1,
                     param_p, param_Z, param_genX, param_error)


# Specification
set.seed(1234)
model_type = "homogenous"
param_error = c(1, 1.5, -0.6)
param_Z = c(0, 0, 0, 0, 1.5, 1.5, 0.9)
param_p = c(0, -0.5, 0.5, 0, 0, 0)
param_y0 = c(3.2, 0.8, 0, 0)
param_y1 = c(3.2+0.4, 0.5, 0, 0)
param_genX = c(0.4, 0, 2)

roydata2 = simul_data(N, model_type, param_y0, param_y1, param_p, param_Z, param_genX, param_error)


usethis::use_data(roydata, roydata2, overwrite = T)
