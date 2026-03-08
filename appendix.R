# =========
# Appendix
# =========

library(sf)
library(raster)
library(ncdf4)
library(dplyr)
library(nlme)
library(ggplot2)
library(MASS)

# ============================================================
# 1. Data Preparation (4 Countries)
# ============================================================

uk1 <- st_read("GBR_ADM1/geoBoundaries-GBR-ADM1.shp")
uk2 <- st_read("GBR_ADM2/geoBoundaries-GBR-ADM2.shp")

ncin <- nc_open("temp_data_SRS.nc")
tmp_brick <- brick("temp_data_SRS.nc", varname = "tmp")
time_raw <- ncvar_get(ncin, "time")
dates <- as.Date("1900-01-01") + time_raw

extract_temp_adm1 <- function(country_name) {
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

df_uk_4 <- bind_rows(
  extract_temp_adm1("Scotland"),
  extract_temp_adm1("England"),
  extract_temp_adm1("Wales"),
  extract_temp_adm1("Northern Ireland")
)

df_annual_4 <- df_uk_4 %>%
  group_by(region, year) %>%
  summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")
df_annual_4$year_c <- df_annual_4$year - mean(df_annual_4$year)
df_annual_4$region <- factor(df_annual_4$region)

# ============================================================
# 2. Model Fitting (4 Countries)
# ============================================================

M_final_4 <- gls(temp ~ year_c + I(year_c^2) * region,
                 data = df_annual_4,
                 correlation = corAR1(form = ~year | region))

# ============================================================
# 3. Figure 1: Trends + Acceleration Coefficients
# ============================================================

dir.create("Fig", showWarnings = FALSE)

df_annual_4$fit <- predict(M_final_4)

theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    strip.text = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

# Figure 1 Left: Trends
p1 <- ggplot(df_annual_4, aes(x = year, y = temp)) +
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

ggsave("Fig/figure1_trends.png", p1, width = 8, height = 5, dpi = 300)

# Figure 1 Right: Acceleration Coefficients
coef_vec <- coef(M_final_4)
vcov_mat <- vcov(M_final_4)
beta2 <- coef_vec["I(year_c^2)"]

regions <- c("England", "Northern Ireland", "Scotland", "Wales")

accel_results <- lapply(regions, function(r) {
  if (r == "England") {
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

p2 <- ggplot(accel_df, aes(x = estimate, y = region)) +
  geom_point(size = 3) +
  geom_errorbar(aes(xmin = lower, xmax = upper),
                orientation = "y", linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  labs(
    title = "Estimated Acceleration by 4 Countries",
    subtitle = expression("Quadratic term (" * beta[2] * ") with 95% CI"),
    x = "Acceleration Coefficient",
    y = ""
  ) +
  theme_pub

ggsave("Fig/figure2_acceleration.png", p2, width = 7, height = 4, dpi = 300)

# ============================================================
# 4. Table 1: Sensitivity (4 Countries vs 8 Sub-regions)
# ============================================================

# 8 Sub-regions
regions_8 <- c("Highland", "Aberdeenshire", "Dumfries and Galloway",
               "Cumbria", "Norfolk", "Devon", "Powys", "Fermanagh and Omagh")

extract_temp_adm2 <- function(region_name) {
  poly <- uk2[uk2$shapeName == region_name, ]
  poly <- st_transform(poly, crs(tmp_brick))
  r_crop <- crop(tmp_brick, poly)
  r_mask <- mask(r_crop, poly)
  monthly <- cellStats(r_mask, stat = "mean", na.rm = TRUE)
  data.frame(
    date = dates,
    year = as.numeric(format(dates, "%Y")),
    month = as.numeric(format(dates, "%m")),
    temp = monthly,
    region = region_name
  )
}

df_uk_8 <- bind_rows(lapply(regions_8, extract_temp_adm2))

df_annual_8 <- df_uk_8 %>%
  group_by(region, year) %>%
  summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")
df_annual_8$year_c <- df_annual_8$year - mean(df_annual_8$year)
df_annual_8$region <- factor(df_annual_8$region)

# 4-country models (ML)
M1_4 <- gls(temp ~ year_c, data = df_annual_4, method = "ML")
M2_4 <- gls(temp ~ year_c + I(year_c^2), data = df_annual_4, method = "ML")
M3_4_no_ar <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_4, method = "ML")
M3_4 <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_4, method = "ML",
            correlation = corAR1(form = ~year | region))

# 8-region models (ML)
M1_8 <- gls(temp ~ year_c, data = df_annual_8, method = "ML")
M2_8 <- gls(temp ~ year_c + I(year_c^2), data = df_annual_8, method = "ML")
M3_8_no_ar <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_8, method = "ML")
M3_8 <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_8, method = "ML",
            correlation = corAR1(form = ~year | region))

# Print table results
cat("\n===== Table 1: Model Comparison =====\n")
cat("\n--- 4 Countries ---\n")
anova(M1_4, M2_4)
anova(M2_4, M3_4_no_ar)
anova(M3_4_no_ar, M3_4)

cat("\n--- 8 Sub-regions ---\n")
anova(M1_8, M2_8)
anova(M2_8, M3_8_no_ar)
anova(M3_8_no_ar, M3_8)

# Final beta2 p-values
M_final_4_reml <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_4,
                      correlation = corAR1(form = ~year | region))
M_final_8_reml <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_8,
                      correlation = corAR1(form = ~year | region))
cat("\nFinal beta2 (4 Countries): p =",
    summary(M_final_4_reml)$tTable["I(year_c^2)", "p-value"], "\n")
cat("Final beta2 (8 Sub-regions): p =",
    summary(M_final_8_reml)$tTable["I(year_c^2)", "p-value"], "\n")

# ============================================================
# 5. Figure 2: Warming Rate Map (all ADM2)
# ============================================================

all_regions <- unique(uk2$shapeName)

extract_slope <- function(region_name) {
  poly <- uk2[uk2$shapeName == region_name, ]
  poly <- st_transform(poly, crs(tmp_brick))
  r_crop <- crop(tmp_brick, poly)
  r_mask <- mask(r_crop, poly)
  monthly <- cellStats(r_mask, stat = "mean", na.rm = TRUE)
  df <- data.frame(
    year = as.numeric(format(dates, "%Y")),
    temp = monthly
  )
  annual <- df %>%
    group_by(year) %>%
    summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")
  if (all(is.na(annual$temp))) return(NULL)
  slope <- coef(lm(temp ~ year, data = annual))[2] * 100
  data.frame(shapeName = region_name, slope = slope)
}

cat("Computing warming rates for all ADM2 regions...\n")
slopes_all <- do.call(rbind, lapply(all_regions, function(r) {
  tryCatch(extract_slope(r), error = function(e) NULL)
}))
cat("Done!\n")

uk2_map <- merge(uk2, slopes_all, by = "shapeName", all.x = TRUE)

ggplot() +
  geom_sf(data = uk2_map, aes(fill = slope), color = "grey80", size = 0.05) +
  geom_sf(data = uk1[uk1$shapeName == "England", ],
          fill = NA, color = "blue", linewidth = 0.8) +
  geom_sf(data = uk1[uk1$shapeName == "Scotland", ],
          fill = NA, color = "red", linewidth = 0.8) +
  geom_sf(data = uk1[uk1$shapeName == "Wales", ],
          fill = NA, color = "green", linewidth = 0.8) +
  geom_sf(data = uk1[uk1$shapeName == "Northern Ireland", ],
          fill = NA, color = "purple", linewidth = 0.8) +
  scale_fill_gradient(low = "lightyellow", high = "darkred",
                      na.value = "grey90",
                      name = "°C / century") +
  labs(title = "Warming Rate across UK Regions (1900-2011)") +
  theme_minimal()

ggsave("Fig/Warming_Rate.png", width = 5, height = 6, dpi = 300)

print(p1)
print(p2)

# four country
M1_reml <- gls(temp ~ year_c, data = df_annual_4)
M2_reml <- gls(temp ~ year_c + I(year_c^2), data = df_annual_4)
M3_no_ar_reml <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_4)

ss_tot_4 <- sum((df_annual_4$temp - mean(df_annual_4$temp))^2)

models <- list(M1 = M1_reml, M2 = M2_reml, M3 = M3_no_ar_reml, M3_AR1 = M_final_4)

for (name in names(models)) {
  m <- models[[name]]
  ss_res <- sum(residuals(m)^2)
  r2 <- 1 - ss_res / ss_tot_4
  rmse <- sqrt(mean(residuals(m)^2))
  cat(name, ": R² =", round(r2, 3), ", RMSE =", round(rmse, 3), "°C\n")
}

# eight sub-region
M1_8_reml <- gls(temp ~ year_c, data = df_annual_8)
M2_8_reml <- gls(temp ~ year_c + I(year_c^2), data = df_annual_8)
M3_8_no_ar_reml <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_8)
M_final_8 <- gls(temp ~ year_c + I(year_c^2) * region, data = df_annual_8,
                 correlation = corAR1(form = ~year | region))

ss_tot_8 <- sum((df_annual_8$temp - mean(df_annual_8$temp))^2)

models_8 <- list(M1 = M1_8_reml, M2 = M2_8_reml, M3 = M3_8_no_ar_reml, M3_AR1 = M_final_8)

for (name in names(models_8)) {
  m <- models_8[[name]]
  ss_res <- sum(residuals(m)^2)
  r2 <- 1 - ss_res / ss_tot_8
  rmse <- sqrt(mean(residuals(m)^2))
  cat(name, "(8-region): R² =", round(r2, 3), ", RMSE =", round(rmse, 3), "°C\n")
}

