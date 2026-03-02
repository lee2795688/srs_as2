library(sf)
library(raster)
library(ncdf4)
library(dplyr)
library(nlme)
library(ggplot2)
library(lmtest)

# 1. Data Preparation

uk1 <- st_read("GBR_ADM1/geoBoundaries-GBR-ADM1.shp")

ncin <- nc_open("temp_data_SRS.nc")
tmp_brick <- brick("temp_data_SRS.nc", varname = "tmp")
time_raw <- ncvar_get(ncin, "time")
dates <- as.Date("1900-01-01") + time_raw

extract_temp <- function(country_name) {
  poly <- uk1[uk1$shapeName == country_name, ]
  poly <- st_transform(poly, crs(tmp_brick))
  r_crop <- crop(tmp_brick, poly)
  r_mask <- mask(r_crop, poly)
  monthly <- cellStats(r_mask, stat = "mean", na.rm = TRUE)
  data.frame(
    date = dates,
    year = as.numeric(format(dates, "%Y")),
    month = as.numeric(format(dates, "%m")),
    temp = monthly,
    region = country_name
  )
}

df_uk <- bind_rows(
  extract_temp("Scotland"),
  extract_temp("England"),
  extract_temp("Wales"),
  extract_temp("Northern Ireland")
)

# Annual mean temperature by region
df_annual <- df_uk %>%
  group_by(region, year) %>%
  summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")

df_annual$year_c <- df_annual$year - mean(df_annual$year)
df_annual$region <- factor(df_annual$region)

# 2. Model Fitting (ML for comparison)

# Model 1: Linear trend
M1 <- gls(temp ~ year_c, data = df_annual, method = "ML")
summary(M1)


# Model 2: Quadratic trend (acceleration)
M2 <- gls(temp ~ year_c + I(year_c^2), data = df_annual, method = "ML")
summary(M2)


# Model 3 (no AR): Quadratic + regional differences
M3_no_ar <- gls(temp ~ year_c + I(year_c^2) * region,
                data = df_annual, method = "ML")
summary(M3_no_ar)

# Model 3: Quadratic + regional differences + AR(1)
M3 <- gls(temp ~ year_c + I(year_c^2) * region,
          data = df_annual, method = "ML",
          correlation = corAR1(form = ~year | region))
summary(M3)


# 3. Model Comparison

anova(M1, M2)        # Is acceleration significant?
anova(M2, M3_no_ar)  # Are regional differences significant?
anova(M3_no_ar, M3)  # Is AR(1) necessary?
AIC(M1, M2, M3_no_ar, M3)

# 4. Final Model (REML for accurate estimates)

M_final <- gls(temp ~ year_c + I(year_c^2) * region,
               data = df_annual,
               correlation = corAR1(form = ~year | region))
summary(M_final)

# 5. Residual Diagnostics

par(mfrow = c(1, 2))
plot(fitted(M_final), residuals(M_final, type = "normalized"),
     xlab = "Fitted", ylab = "Normalized Residuals",
     main = "Residuals vs Fitted")
abline(h = 0, lty = 2)

acf(residuals(M_final, type = "normalized"),
    main = "ACF: Final Model (GLS AR(1))")



#-----------------
# UK overall trend
#-----------------

ggplot(df_annual, aes(year, temp)) +
  geom_point(alpha=0.3) +
  stat_smooth(method="lm", formula = y ~ poly(x,2), se=TRUE) +
  theme_minimal()


#---------
# facet
#---------

ggplot(df_annual, aes(year, temp)) +
  geom_point(alpha=0.2) +
  stat_smooth(method="lm", formula = y ~ poly(x,2), se=FALSE) +
  facet_wrap(~region) +
  theme_minimal()


#--------------------------
# acceleration coefficient
#--------------------------

install.packages("broom.mixed")
library(broom.mixed)

coef_df <- tidy(M_final)

#######################################
library(ggplot2)
library(dplyr)
library(sf)
library(viridis)

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "right"
  )

#---------------------------------------
# Figure1-Quadratic GLS Trend by Region
#---------------------------------------
library(ggplot2)
library(dplyr)

# predicted values
df_annual$fit <- predict(M_final)

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

p1 <- ggplot(df_annual, aes(x = year, y = temp)) +
  geom_point(alpha = 0.25, size = 0.8, color = "grey50") +
  geom_line(aes(y = fit), linewidth = 1.2, color = "black") +
  facet_wrap(~region, ncol = 2) +
  labs(
    title = "Nonlinear Warming Trends across the United Kingdom",
    subtitle = "GLS model with AR(1) errors",
    x = "Year",
    y = "Annual Mean Temperature (°C)"
  ) +
  theme_pub

ggsave("figure1_trends.png", p1,
       width = 8, height = 5, dpi = 300)

#-----------------------------------------
# Figure 2 – Acceleration Coefficient Plot
#-----------------------------------------
library(MASS)

# extract coefficients and covariance matrix
coef_vec <- coef(M_final)
vcov_mat <- vcov(M_final)

# baseline acceleration
beta2 <- coef_vec["I(year_c^2)"]

# total acceleration for each regions
regions <- c("England", "Northern Ireland", "Scotland", "Wales")

accel_results <- lapply(regions, function(r) {
  
  if (r == "England") {
    # baseline
    est <- beta2
    var <- vcov_mat["I(year_c^2)", "I(year_c^2)"]
  } else {
    term <- paste0("I(year_c^2):region", r)
    est <- beta2 + coef_vec[term]
    
    var <- vcov_mat["I(year_c^2)", "I(year_c^2)"] +
      vcov_mat[term, term] +
      2 * vcov_mat["I(year_c^2)", term]
  }
  
  se <- sqrt(var)
  
  data.frame(
    region = r,
    estimate = est,
    lower = est - 1.96 * se,
    upper = est + 1.96 * se
  )
})

accel_df <- do.call(rbind, accel_results)

accel_df$region <- factor(accel_df$region,
                          levels = accel_df$region[order(accel_df$estimate)])

p2 <- ggplot(accel_df,
             aes(x = estimate,
                 y = region)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y",
                linewidth = 0.8) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             color = "red",
             linewidth = 0.8) +
  labs(
    title = "Estimated Acceleration by 4 Countries",
    subtitle = "Quadratic term (β₂) with 95% Confidence Intervals",
    x = "Acceleration Coefficient",
    y = ""
  ) +
  theme_pub

ggsave("figure2_acceleration.png", p2,
       width = 7, height = 4, dpi = 300)


######################################################


