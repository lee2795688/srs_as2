library(sf)
library(raster)
library(ncdf4)
library(dplyr)
library(ggplot2)

uk <- st_read("GBR_ADM1/geoBoundaries-GBR-ADM1.shp")
scotland <- uk[uk$shapeName == "Scotland", ]


tmp_brick <- brick("temp_data_SRS.nc", varname = "tmp")

scot_brick <- crop(tmp_brick, scotland)
scot_brick <- mask(scot_brick, scotland)

scot_monthly <- cellStats(scot_brick, stat = "mean", na.rm = TRUE)

time_raw <- ncvar_get(ncin, "time")
dates <- as.Date("1900-01-01") + time_raw

df_scot <- data.frame(
  date = dates,
  year = as.numeric(format(dates, "%Y")),
  month = as.numeric(format(dates, "%m")),
  temp = scot_monthly
)

head(df_scot)
tail(df_scot)

plot(df_scot$date, df_scot$temp, type = "l",
     xlab = "Year", ylab = "Temperature (°C)",
     main = "Scotland Monthly Mean Temperature")

##heat map

climatology <- df_scot %>%
  group_by(month) %>%
  summarise(clim = mean(temp, na.rm = TRUE))

print(climatology)

df_scot$clim <- NULL
df_scot$anomaly <- NULL

df_scot <- df_scot %>%
  left_join(climatology, by = "month") %>%
  mutate(anomaly = temp - clim)

head(df_scot)


ggplot(df_scot, aes(x = year, y = factor(month), fill = anomaly)) +
  geom_tile() +
  scale_fill_gradient2(low = "blue", mid = "white", high = "red",
                       midpoint = 0, name = "Anomaly (°C)") +
  scale_y_discrete(labels = month.abb) +
  labs(title = "Scotland Temperature Anomaly Heatmap",
       x = "Year", y = "Month") +
  theme_minimal()

df_scot$season <- factor(
  ifelse(df_scot$month %in% c(12,1,2), "Winter",
         ifelse(df_scot$month %in% c(3,4,5), "Spring",
                ifelse(df_scot$month %in% c(6,7,8), "Summer", "Autumn"))),
  levels = c("Winter","Spring","Summer","Autumn")
)

df_scot$season_year <- ifelse(df_scot$month == 12, 
                              df_scot$year + 1, df_scot$year)

seasonal <- df_scot %>%
  group_by(season_year, season) %>%
  summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop")
ggplot(seasonal, aes(x = season, y = temp, fill = season)) +
  geom_violin(alpha = 0.7, trim = FALSE) +
  geom_boxplot(width = 0.15, alpha = 0.5) +
  labs(title = "Temperature Distribution by Season in Scotland",
       x = "Season", y = "Temperature (°C)") +
  theme_minimal()

ggplot(seasonal, aes(x = season_year, y = temp, color = season)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE) +
  facet_wrap(~season, scales = "free_y") +
  labs(title = "Scotland: Seasonal Temperature Trends (1900-2011)",
       x = "Year", y = "Temperature (°C)") +
  theme_minimal()

##M1: temp ~ year_c
##M2: temp ~ year_c + season
##M3: temp ~ year_c * season
##M4: temp ~ (year_c + I(year_c²)) * season