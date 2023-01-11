## DEFINE CONSTANTS

# "The IPCC authors concluded that ECS is very likely to be greater than 1.5 °C
# (2.7 °F) and likely to lie in the range 2 to 4.5 °C (3.6 to 8.1 °F), with a
# most likely value of about 3 °C (5.4 °F)."
ecs <- 3   # deg
alpha <- -3.93 / ecs   # deg^-1

# Density of seawater
rho <- 1027   # kg / m^3

# Specific heat capacity of water
cp <- 4218  # J / kg / deg

# Depth of upper and deep layers
hu <- 100  # m  (range 25-100)
hd <- 900  # m  (range 500-4000)

# Vertical diffusivity
kappa <- 1  # cm^2 / s   (range 0.2 - 2)
kappa <- kappa * 10e-5   # convert to m^2 / s
kappa <- kappa * 60*60*24*365   # convert to m^2 / year

# Heat diffusion coefficient
gamma <- 2 * kappa * cp * rho / (hu + hd)   # J / deg / m^2 / year

# Thermal inertia for the upper and deep layers
cu <- rho * cp * hu   # J / m^2 / deg
cd <- rho * cp * hd   # J / m^2 / deg

## DO THE SIMULATION
# Equations from "The inconstancy of the transient climate response
# parameter under increasing CO2", Gregory et al., 2017.
simulate <- function(
    # A vector containing total forcing values for the time period of interest.
    forcing,
    # Climate feedback parameter.
    alpha,
    # Parameter to capture heat diffusion into the deep ocean.
    gamma,
    # Thermal inertia parameters.
    cu, cd
) {
    # Get the number of iterations.
    years <- length(forcing)
    # Set up the vectors for the results.
    tu <- td <- rep(NA, years + 1)
    # Anomaly for the first year is 0, this is the reference year.
    tu[1] <- td[1] <- 0
    # Iterate.
    for (i in 1:years) {
        # Update the upper and deep layers.
        tu[i + 1] <-
            tu[i] + (forcing[i] + alpha * tu[i] - gamma * (tu[i] - td[i])) / cu
        td[i + 1] <- td[i] + (gamma * (tu[i] - td[i])) / cd
    }
    return(data.frame(
        year = 1750 + 0:years,
        tu = tu,
        td = td
    ))
}

## RUN THE SIMULATION
# This is the forcing data for the SSP119 scenario from 1750-2500.
data <- read.csv("data/SSPs/ERF_ssp585_1750-2500.csv")

# Use the total forcing values
forcing <- data$total

# Run the simulation.
out <- simulate(forcing, alpha, gamma, cu, cd)

# Plot the temperatures of the upper layer.
plot(x = out$year, y = out$tu)
