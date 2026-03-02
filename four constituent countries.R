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

# Model 2: Quadratic trend (acceleration)
M2 <- gls(temp ~ year_c + I(year_c^2), data = df_annual, method = "ML")

# Model 3 (no AR): Quadratic + regional differences
M3_no_ar <- gls(temp ~ year_c + I(year_c^2) * region,
                data = df_annual, method = "ML")

# Model 3: Quadratic + regional differences + AR(1)
M3 <- gls(temp ~ year_c + I(year_c^2) * region,
          data = df_annual, method = "ML",
          correlation = corAR1(form = ~year | region))

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

