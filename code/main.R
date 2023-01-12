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

## SIMULATION
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

## RUN THE SIMULATION FOR SSPS
run_for_ssps <- function(
    # Title of the plot output.
    plot_title,
    # Function to extract the forcing values from the data. Defaults to just
    # using the total forcing values.
    get_forcings = function(x, scenario) { x$total },
    # Suffix to add to the output filenames. Defaults to an empty string.
    file_suffix = "",
    # Set a custom range for the y axis of the plot. Defaults to scaling with
    # the data.
    ylim = NULL
) {
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
        # Extract the values for the forcing.
        forcing <- get_forcings(data, ssp)
        # Run the simulation.
        results[[ssp]] <- simulate(
            forcing, alpha, gamma, cu, cd, tr = 60 * 60 * 24 * 365
        )
    }

    # Plot the temperatures of the upper layer, and save the data.
    pdf(file = paste("output/projections", file_suffix, ".pdf", sep = ""))
    n_lines <- length(results)
    colors <- rainbow(n_lines)
    plot(
        1, type = "n",
        xlim = range(results$ssp585$year),
        ylim = if (is.null(ylim)) range(results$ssp585$tu) else ylim,
        xlab = "Year", ylab = "Temperature anomaly (relative to 1750)",
        main = plot_title
    )
    for (i in 1:n_lines) {
        lines(results[[i]]$year, results[[i]]$tu, type = "l", col = colors[i])
        # Save the data.
        write.csv(
            results[[i]], file = paste(
                "output/", names(results)[i],
                file_suffix, ".csv",
                sep = ""
            ), quote = FALSE, row.names = FALSE
        )
    }
    legend("topleft", legend = names(results), fill = colors)
    dev.off()
    results
}

results <- run_for_ssps(
    "Model temperature anomaly projections, relative to 1750",
    file_suffix = "_default"
)

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

## RUN THE SCENARIO

# CO2 mol. mass 44x10^-3 kg/mol
# (https://en.wikipedia.org/wiki/Carbon_dioxide)
co2_mm <- 44e-3   # kg / mol

# approx. mass of atmosphere: 5.14×10^18 kg
# https://homework.study.com/explanation/the-mass-of-the-earth-s-atmosphere-is-estimated-as-5-14-x-1015-t-1-t-1000-kg-the-average-molar-mass-of-air-is-28-8-g-mol-how-many-moles-of-gas-are-in-the-atmosphere.html
atm_m <- 5.14e18  # kg

# mean mol. mass of air 28.8x10^-3 kg/mol
# https://homework.study.com/explanation/the-mass-of-the-earth-s-atmosphere-is-estimated-as-5-14-x-1015-t-1-t-1000-kg-the-average-molar-mass-of-air-is-28-8-g-mol-how-many-moles-of-gas-are-in-the-atmosphere.html
air_mm <- 28.8e-3 # kg / mol

# Total moles of molecules in the atmosphere
tot_molecules <- atm_m / air_mm

# Total moles of molecules in one ton of CO2
co2_molecules <- 1000 / co2_mm

# ppm of one ton of CO2 in the atmosphere
co2_ppm <- co2_molecules * 1e6 / tot_molecules

# Radiative forcing conversion
# https://en.wikipedia.org/wiki/Radiative_forcing#Carbon_dioxide
get_rf_delta <- function(tons_co2_equ) {
    co2_0 <- 278   # ppm (est. in 1750)
    5.35 * log((co2_0 + co2_ppm * tons_co2_equ) / co2_0)
}

# Vegan diet expected to save 0.4 - 2.1 tCO2eq per capita (per year?)
# from page 803 of:
# https://www.ipcc.ch/report/ar6/wg3/downloads/report/IPCC_AR6_WGIII_FullReport.pdf
veg_t_saved <- 1.2  # tCO2eq / capita / year

# Read population projections.
pop_data <- read.csv("data/iamc_db.csv")[-(11:15), ]
pop_data <- subset(
    pop_data,
    # Use the projections from OECD
    subset = Model == "OECD Env-Growth",
    select = c(Scenario, X2010:X2100)
)
pop_data <- data.frame(
    year = seq(2010, 2100, 5),
    ssp1 = as.numeric(pop_data[pop_data$Scenario == "SSP1", -1]),
    ssp2 = as.numeric(pop_data[pop_data$Scenario == "SSP2", -1]),
    ssp3 = as.numeric(pop_data[pop_data$Scenario == "SSP3", -1]),
    ssp4 = as.numeric(pop_data[pop_data$Scenario == "SSP4", -1]),
    ssp5 = as.numeric(pop_data[pop_data$Scenario == "SSP5", -1])
)

# Create a dataframe that has values per year.
pop_data <- pop_data[-(1:2),]
pop_data <- pop_data[rep(seq_len(nrow(pop_data)), each = 5),][-(1:3),]
# Extend this to 2500 -- by repeating the last value (i.e. assume
# population plateaus)
pop_data <- rbind(pop_data, pop_data[rep(nrow(pop_data), 396),])
pop_data$year <- 2023:2500

# Run, and account for lower radiative forcing from veganism.
run_for_ssps(
    "Model temperature anomaly projections, under a \nvegan people scenario",
    get_forcings = function(data, scenario) {
        # Define a date range to look at.
        rnge <- 2023:2104
        # Drop unnecessary columns.
        reduced <- data$total
        # Get future and past values.
        future <- reduced[data$year  %in% rnge]
        past <- reduced[data$year <= 2022]
        # Scale future emissions.
        future <- future - get_rf_delta(
            veg_t_saved * pop_data[
                pop_data$year %in% rnge,
                substr(scenario, start = 1, stop = 4)
            ] * 1e6   # the data is in millions
        )
        # Return the total values.
        c(past, future)
    },
    file_suffix = "_scenario",
    ylim = c(0,4)
)

# Get a reduced default, to compare.
results <- run_for_ssps(
    "Model temperature anomaly projections, relative to 1750",
    get_forcings = function(data, scen) {
        data$total[data$year < 2105]
    },
    file_suffix = "_default_red",
    ylim = c(0,4)
)
