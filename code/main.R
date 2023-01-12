## DEFINE CONSTANTS

# "The IPCC authors concluded that ECS is very likely to be greater than 1.5 °C
# (2.7 °F) and likely to lie in the range 2 to 4.5 °C (3.6 to 8.1 °F), with a
# most likely value of about 3 °C (5.4 °F)."
ecs <- 3   # deg
alpha <- -3.93 / ecs   # W / m^2 / deg

# Density of seawater
rho <- 1027   # kg / m^3

# Specific heat capacity of water
cp <- 4218  # J / kg / deg

# Depth of upper and deep layers
hu <- 100  # m  (range 25-100)
hd <- 900  # m  (range 500-4000)

# Vertical diffusivity
kappa <- 1  # cm^2 / s   (range 0.2 - 2)
kappa <- kappa * 10e-4   # convert to m^2 / s

# Heat diffusion coefficient
gamma <- 2 * kappa * cp * rho / (hu + hd)   # J / deg / m^2 / s

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
    cu, cd,
    # Temporal resolution -- default 5 minutes
    tr = 60 * 5
) {
    # Get the number of iterations.
    years <- length(forcing)
    # Set up the vectors for the results.
    tu <- td <- rep(NA, years + 1)
    # Anomaly for the first year is 0, this is the reference year.
    tu[1] <- td[1] <- 0
    # Define seconds per year.
    spy <- 60 * 60 * 24 * 365
    # Iterate.
    for (i in 1:years) {
        tu_cur <- tu[i]
        td_cur <- td[i]
        # Take steps for every window of the simulation resolution. We keep
        # the forcing constant across the year.
        for (j in 1:(spy / tr)) {
            # Update the upper and deep layers.
            tu_cur <-
                tu_cur + tr * (forcing[i] + alpha * tu_cur - gamma * (tu_cur - td_cur)) / cu
            td_cur <- td_cur + tr * (gamma * (tu_cur - td_cur)) / cd
        }
        tu[i + 1] <- tu_cur
        td[i + 1] <- td_cur
    }
    return(data.frame(
        year = 1750 + 0:years,
        tu = tu,
        td = td
    ))
}

## RUN THE SIMULATION
ssps <- c(
    "ssp119", "ssp126", "ssp245", "ssp370", "ssp434", "ssp460", "ssp585"
)

results <- list()
for (ssp in ssps) {
    # Read the forcing data
    data <- read.csv(paste(
        "data/SSPs/ERF_", ssp, "_1750-2500.csv",
        sep = ""
    ))
    # Use the total forcing values
    forcing <- data$total
    # Run the simulation.
    results[[ssp]] <- simulate(
        forcing, alpha, gamma, cu, cd, tr = 60 * 60 * 24 * 365
    )
}

# Plot the temperatures of the upper layer, and save the data.
pdf(file = "output/projections.pdf")
n_lines <- length(results)
colors <- rainbow(n_lines)
plot(
    1, type = "n",
    xlim = c(1750, 2501), ylim = range(results$ssp585$tu),
    xlab = "Year", ylab = "Temperature anomaly (relative to 1750)",
    main = "Model temperature anomaly projections, relative to 1750"
)
for (i in 1:n_lines) {
    lines(results[[i]]$year, results[[i]]$tu, type = "l", col = colors[i])
    # Save the data.
    write.csv(
        result, file = paste(
            "output/", names(results)[i], ".csv",
            sep = ""
        ), quote = FALSE, row.names = FALSE
    )
}
legend("topleft", legend = names(results), fill = colors)
dev.off()

## COMPARE THE MODEL TO MEASUREMENTS
out <- results[1]

# This is the MET office data.
validation_data <- read.csv("data/HadCRUT.5.0.1.0.summary_series.global.annual.csv")

# This is our data in the same range.
model_data <- out$tu[out$year %in% validation_data$Time]

# MET data is anomaly relative to 1961-1990, so let's correct our outputs to that.
baseline <- mean(out$tu[out$year %in% 1961:1990])
model_data <- model_data - baseline

# Plot.
pdf(file = "output/validate.pdf")
plot(
    1850:2022, validation_data$Anomaly..deg.C.,
    col = "red",
    xlab = "Year",
    ylab = "Anomaly (relative to 1961-1990)",
    main = "Temperature anomaly relative to 1961-1990"
)
lines(1850:2022, model_data, col = "blue")
legend(
    "topleft",
    legend = c("MET measurements", "Model prediction"),
    fill = c("red", "blue")
)
dev.off()