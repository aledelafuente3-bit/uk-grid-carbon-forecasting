library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(car)
library(rstatix)

# Load data
df <- read_csv("Carbon_Intensity_Data.csv", show_col_types = FALSE)

# Rename columns to R-friendly names
df <- df %>%
  rename(
    datetime_utc = `Datetime (UTC)`,
    actual = `Actual Carbon Intensity (gCO2/kWh)`,
    forecast = `Forecast Carbon Intensity (gCO2/kWh)`,
    index = Index
  ) %>%
  mutate(
    datetime_utc = ymd_hms(datetime_utc, tz = "UTC"),
    hour = hour(datetime_utc) + minute(datetime_utc) / 60,
    dow = wday(datetime_utc, label = TRUE, week_start = 1),
    error = actual - forecast
  ) %>%
  arrange(datetime_utc)

# Sanity checks
glimpse(df)
summary(df$actual)
nrow(df)

library(tidyverse)
library(lubridate)

raw <- read_csv("Carbon_Intensity_Data.csv", show_col_types = FALSE)

bad <- raw %>%
  mutate(parsed = ymd_hms(`Datetime (UTC)`, tz = "UTC")) %>%
  filter(is.na(parsed)) %>%
  count(`Datetime (UTC)`, sort = TRUE)

bad
library(tidyverse)
library(lubridate)

df <- read_csv("Carbon_Intensity_Data.csv", show_col_types = FALSE) %>%
  rename(
    datetime_utc = `Datetime (UTC)`,
    actual = `Actual Carbon Intensity (gCO2/kWh)`,
    forecast = `Forecast Carbon Intensity (gCO2/kWh)`,
    index = Index
  ) %>%
  mutate(
    # datetime_utc is already <dttm>, just ensure timezone is UTC
    datetime_utc = with_tz(datetime_utc, tzone = "UTC"),
    hour = hour(datetime_utc) + minute(datetime_utc) / 60,
    dow = wday(datetime_utc, label = TRUE, week_start = 1),
    error = actual - forecast
  ) %>%
  arrange(datetime_utc)

# Checks
sum(is.na(df$datetime_utc))  # should be 0
nrow(df)                     # should be 1443
df %>%
  mutate(diff_mins = as.numeric(difftime(datetime_utc, lag(datetime_utc), units = "mins"))) %>%
  filter(!is.na(diff_mins)) %>%
  count(diff_mins, sort = TRUE)
# Create time series (half-hourly = 48 per day)
y <- ts(df$actual, frequency = 48)

# STL decomposition
fit_stl <- stl(y, s.window = "periodic")

# Extract components
t_index <- time(y)  # in "days" because frequency = 48 per day
obs     <- fit_stl$time.series[, "remainder"] + fit_stl$time.series[, "seasonal"] + fit_stl$time.series[, "trend"]
season  <- fit_stl$time.series[, "seasonal"]
trend   <- fit_stl$time.series[, "trend"]
rem     <- fit_stl$time.series[, "remainder"]

# Plot with clear labels
op <- par(mfrow = c(4,1), mar = c(3,4,2,1), oma = c(3,0,3,0))

plot(t_index, obs, type = "l",
     xlab = "", ylab = "gCO2/kWh",
     main = "Observed carbon intensity")

plot(t_index, season, type = "l",
     xlab = "", ylab = "gCO2/kWh",
     main = "Seasonal (diurnal, period = 48)")

plot(t_index, trend, type = "l",
     xlab = "", ylab = "gCO2/kWh",
     main = "Trend (multi-day)")

plot(t_index, rem, type = "l",
     xlab = "", ylab = "gCO2/kWh",
     main = "Remainder")

mtext("Time (days from start of series; UTC)", side = 1, outer = TRUE, line = 1)
mtext("STL decomposition of half-hourly grid carbon intensity (Actual)", side = 3, outer = TRUE, line = 1, font = 2)

par(op)


# Autocorrelation
Acf(y, main = "ACF: Actual Carbon Intensity")
Pacf(y, main = "PACF: Actual Carbon Intensity")

# Stationarity test
adf.test(y)

# Represent daily cycle smoothly
df <- df %>%
  mutate(
    hour_rad = 2 * pi * hour / 24,
    sin_hour = sin(hour_rad),
    cos_hour = cos(hour_rad)
  )

# Models
m1 <- lm(actual ~ forecast, data = df)
m2 <- lm(actual ~ forecast + sin_hour + cos_hour + dow, data = df)

anova(m1, m2)
summary(m2)

# Diagnostics
par(mfrow = c(2,2))
plot(m2)
par(mfrow = c(1,1))

df %>%
  summarise(
    MAE = mean(abs(error)),
    RMSE = sqrt(mean(error^2)),
    Bias = mean(error)
  )
# Group sizes
df %>% count(index)

# Variance test
leveneTest(actual ~ index, data = df)

# ANOVA
aov1 <- aov(actual ~ index, data = df)
summary(aov1)

# Effect size
eta_squared(aov1)

# Post-hoc
TukeyHSD(aov1)

kruskal.test(actual ~ index, data = df)

df %>%
  dunn_test(actual ~ index, p.adjust.method = "holm")

# Welch ANOVA
oneway.test(actual ~ index, data = df)

# Games-Howell post-hoc
rstatix::games_howell_test(df, actual ~ index)
pkgs <- c("tidyverse","lubridate","forecast","tseries","broom","car","rstatix",
          "sandwich","lmtest","knitr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))
pkgs <- c("tidyverse","lubridate","forecast","tseries","broom","car","rstatix",
          "sandwich","lmtest","knitr")
to_install <- pkgs[!pkgs %in% installed.packages()[,"Package"]]
if(length(to_install) > 0) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))
# Ensure datetime is POSIXct
df <- df %>% mutate(datetime_utc = as.POSIXct(datetime_utc, tz = "UTC"))

library(tidyverse)
library(lubridate)
library(forecast)
library(tseries)
library(broom)
library(car)
library(rstatix)
library(sandwich)
library(lmtest)
cat("R is responsive\n")
library(tseries)
library(forecast)
library(tibble)

y <- ts(df$actual, frequency = 48)

# ADF test (fast)
adf <- tseries::adf.test(y)

# ACF summary (fast): report the strongest autocorrelation lags
acf_vals <- acf(y, plot = FALSE, lag.max = 96)$acf[-1]   # exclude lag 0
top_lags <- order(abs(acf_vals), decreasing = TRUE)[1:5]
top_lags_tbl <- tibble(
  lag = top_lags,
  acf = round(acf_vals[top_lags], 3)
)

table2 <- tibble(
  metric = c(
    "Time step (minutes)",
    "Seasonal frequency (per day)",
    "ADF statistic",
    "ADF p-value",
    "Max |ACF| (lags 1–96)",
    "Lag of max |ACF|"
  ),
  value = c(
    30,
    48,
    round(unname(adf$statistic), 3),
    as.character(adf$p.value),
    round(max(abs(acf_vals)), 3),
    which.max(abs(acf_vals))
  )
)

print(table2)
print(top_lags_tbl)
library(broom)
library(sandwich)
library(lmtest)
library(dplyr)
library(lubridate)
library(rstatix)
library(tibble)

# Ensure time features exist
df <- df %>%
  mutate(
    hour = hour(datetime_utc) + minute(datetime_utc)/60,
    dow = wday(datetime_utc, label = TRUE, week_start = 1),
    hour_rad = 2*pi*hour/24,
    sin_hour = sin(hour_rad),
    cos_hour = cos(hour_rad)
  )

# ---- Table 3: Regression (robust SE)
m2 <- lm(actual ~ forecast + sin_hour + cos_hour + dow, data = df)
rob <- coeftest(m2, vcov = vcovHC(m2, type = "HC1"))

table3 <- broom::tidy(rob) %>%
  transmute(
    term,
    estimate = round(estimate, 3),
    se_robust = round(std.error, 3),
    t = round(statistic, 2),
    p = format.pval(p.value, digits = 3)
  )

fit_stats <- tibble(
  term = c("N", "R-squared", "Adj R-squared", "Residual SE"),
  estimate = round(c(nobs(m2), summary(m2)$r.squared, summary(m2)$adj.r.squared, summary(m2)$sigma), 3),
  se_robust = NA, t = NA, p = NA
)

table3_out <- bind_rows(table3, fit_stats)
print(table3_out)

# ---- Table 4A: Index summary
table4a <- df %>%
  group_by(index) %>%
  summarise(
    n = n(),
    mean_actual = round(mean(actual), 2),
    sd_actual = round(sd(actual), 2),
    median_actual = round(median(actual), 2),
    mean_forecast = round(mean(forecast), 2),
    mean_error = round(mean(actual - forecast), 2),
    .groups = "drop"
  ) %>%
  arrange(mean_actual)

print(table4a)

# ---- Table 4B: Welch + Games-Howell
welch <- oneway.test(actual ~ index, data = df)
table4b_welch <- tibble(
  test = "Welch one-way ANOVA (actual ~ index)",
  F = round(unname(welch$statistic), 2),
  df_num = unname(welch$parameter[1]),
  df_den = round(unname(welch$parameter[2]), 3),
  p = format.pval(welch$p.value, digits = 3)
)
print(table4b_welch)

gh <- games_howell_test(df, actual ~ index) %>%
  mutate(
    estimate = round(estimate, 2),
    conf.low = round(conf.low, 2),
    conf.high = round(conf.high, 2),
    p.adj = format.pval(p.adj, digits = 3)
  )
print(gh)
write.csv(table1, "Table1.csv", row.names = FALSE)
write.csv(table2, "Table2.csv", row.names = FALSE)
write.csv(table3_out, "Table3.csv", row.names = FALSE)
write.csv(table4a, "Table4A.csv", row.names = FALSE)
write.csv(gh, "Table4B_GamesHowell.csv", row.names = FALSE)
library(forecast)
library(tseries)
library(tibble)

# 1) Make the half-hourly time series
y <- ts(df$actual, frequency = 48)

# 2) ADF test (creates 'adf')
adf <- tseries::adf.test(y)

# 3) Fit ARIMA (FAST so it won't hang)  <-- creates 'fit_arima'
fit_arima <- forecast::auto.arima(
  y,
  seasonal = TRUE,
  stepwise = TRUE,
  approximation = TRUE
)
acc_arima <- as.data.frame(forecast::accuracy(fit_arima))

# 4) Fit seasonal ETS via STLF  <-- creates 'fit_stlf'
fit_stlf <- forecast::stlf(y, method = "ets", h = 48)
acc_stlf <- as.data.frame(forecast::accuracy(fit_stlf))

# 5) Now build Table 2.2  <-- uses the objects above
table2.2 <- tibble(
  metric = c(
    "Time step (minutes)",
    "Seasonal frequency (per day)",
    "ADF statistic",
    "ADF p-value",
    "ARIMA model (fast)",
    "ARIMA AICc",
    "ARIMA RMSE (training)",
    "STLF method",
    "STLF RMSE (training)"
  ),
  value = c(
    30,
    48,
    round(unname(adf$statistic), 3),
    as.character(adf$p.value),
    as.character(fit_arima),
    round(fit_arima$aicc, 2),
    round(acc_arima$RMSE[1], 2),
    "stlf(method='ets')",
    round(acc_stlf$RMSE[1], 2)
  )
)

print(table2.2)
library(forecast)

y <- ts(df$actual, frequency = 48)

# hold out last 2 days (96 half-hours)
h <- 96
train <- head(y, length(y) - h)
test  <- tail(y, h)

fit_arima <- auto.arima(train, seasonal = TRUE, stepwise = TRUE, approximation = TRUE)
fc_arima <- forecast(fit_arima, h = h)

fit_stlf <- stlf(train, method = "ets", h = h)
fc_stlf <- forecast(fit_stlf, h = h)

rmse <- function(actual, pred) sqrt(mean((actual - pred)^2))

rmse_arima <- rmse(as.numeric(test), as.numeric(fc_arima$mean))
rmse_stlf  <- rmse(as.numeric(test), as.numeric(fc_stlf$mean))

tibble(
  model = c("ARIMA (fast)", "STLF (ETS seasonal)"),
  RMSE_test = c(rmse_arima, rmse_stlf)
) |> print()
table2.2 <- tibble(
  metric = c(
    "Time step (minutes)",
    "Seasonal frequency (per day)",
    "ADF statistic",
    "ADF p-value",
    "ARIMA model (fast)",
    "ARIMA AICc (training)",
    "ARIMA RMSE (test, 2-day holdout)",
    "STLF method",
    "STLF RMSE (test, 2-day holdout)"
  ),
  value = c(
    30,
    48,
    -4.97,
    0.01,
    "ARIMA(0,0,3)(0,1,1)[48]",
    10931.17,
    38.9,
    "stlf(method='ets')",
    25.3
  )
)

print(table2.2)

library(forecast)

# Half-hourly series
y <- ts(df$actual, frequency = 48)

# ACF up to 2 days = 96 half-hours
acf_obj <- Acf(y, lag.max = 96, plot = FALSE)

# Convert lags to hours (each lag = 0.5 hours)
lags_hours <- acf_obj$lag[, , 1] * 24  # because Acf lag is in "days" for ts

plot(lags_hours, acf_obj$acf[, , 1], type = "h",
     xlab = "Lag (hours)", ylab = "Autocorrelation",
     main = "Autocorrelation function (ACF) of half-hourly carbon intensity\n(up to 48 hours)",
     ylim = c(min(acf_obj$acf), 1))

abline(h = 0, lty = 2)
# significance bounds (approx)
ci <- 1.96 / sqrt(length(y))
abline(h = c(-ci, ci), lty = 3)
pacf_obj <- Pacf(y, lag.max = 96, plot = FALSE)
lags_hours <- pacf_obj$lag[, , 1] * 24

plot(lags_hours, pacf_obj$acf[, , 1], type = "h",
     xlab = "Lag (hours)", ylab = "Partial autocorrelation",
     main = "Partial autocorrelation function (PACF) of half-hourly carbon intensity\n(up to 48 hours)")

abline(h = 0, lty = 2)
ci <- 1.96 / sqrt(length(y))
abline(h = c(-ci, ci), lty = 3)
pacf_obj <- Pacf(y, lag.max = 96, plot = FALSE)
lags_hours <- pacf_obj$lag[, , 1] * 24

plot(lags_hours, pacf_obj$acf[, , 1], type = "h",
     xlab = "Lag (hours)", ylab = "Partial autocorrelation",
     main = "Partial autocorrelation function (PACF) of half-hourly carbon intensity\n(up to 48 hours)")

abline(h = 0, lty = 2)
ci <- 1.96 / sqrt(length(y))
abline(h = c(-ci, ci), lty = 3)
# Your model
m2 <- lm(actual ~ forecast + sin_hour + cos_hour + dow, data = df)

library(ggplot2)

ggplot(df, aes(x = index, y = actual)) +
  geom_boxplot(outlier.alpha = 0.3) +
  labs(
    title = "Observed grid carbon intensity by Carbon Intensity Index category",
    x = "Index category",
    y = "Actual carbon intensity (gCO2/kWh)"
  )
# Model
m2 <- lm(actual ~ forecast + sin_hour + cos_hour + dow, data = df)

# Give extra outer top space for a main title
op <- par(mfrow = c(2,2),
          mar = c(4,4,2,1),   # inner margins
          oma = c(0,0,4,0))   # outer margins (top = 4 lines)

plot(m2)

mtext("OLS diagnostic plots for carbon-intensity regression")
outer = TRUE, side = 3, line = 1, cex = 1.0, font = 2)

par(op)

# Figure 4: Observed vs Forecast for a 2-day holdout (ARIMA vs STLF)
library(tidyverse)
library(forecast)
library(lubridate)

# --- 1) Create half-hourly ts
y <- ts(df$actual, frequency = 48)

# --- 2) Define 2-day holdout = 96 half-hours
h <- 96
n <- length(y)

train_y <- window(y, end = time(y)[n - h])
test_y  <- window(y, start = time(y)[n - h + 1])

# Matching datetime vector for plotting (assumes df is already ordered by datetime_utc)
test_time <- df$datetime_utc[(n - h + 1):n]

# --- 3) Fit models on training set
# ARIMA (fast + robust)
fit_arima <- auto.arima(train_y, seasonal = TRUE)

# STLF (ETS seasonal) — handles frequency=48 properly
fit_stlf <- stlf(train_y, method = "ets", h = h)

# --- 4) Forecast h steps ahead
fc_arima <- forecast(fit_arima, h = h)
fc_stlf  <- fit_stlf  # stlf() already returns a forecast object

# --- 5) Build tidy dataframe for ggplot
plot_df <- tibble(
  datetime_utc = test_time,
  observed = as.numeric(test_y),
  ARIMA = as.numeric(fc_arima$mean),
  STLF = as.numeric(fc_stlf$mean)
) %>%
  pivot_longer(cols = c(ARIMA, STLF), names_to = "model", values_to = "forecast")

# --- 6) Plot (journal-ish but simple)
ggplot(plot_df, aes(x = datetime_utc)) +
  geom_line(aes(y = observed), linewidth = 0.8) +
  geom_line(aes(y = forecast, linetype = model), linewidth = 0.8) +
  labs(
    title = "Observed vs forecast grid carbon intensity (2-day holdout)",
    subtitle = "Comparison of seasonal ARIMA and STLF (ETS seasonal) forecasts",
    x = "Time (UTC)",
    y = "Carbon intensity (gCO2/kWh)",
    linetype = "Forecast model"
  ) +
  theme_minimal()
# ---- Table 4: OLS regression with robust (HC1) standard errors ----
library(tidyverse)
library(lubridate)
library(broom)
library(sandwich)
library(lmtest)
library(knitr)

# 1) Ensure required columns exist
df <- df %>%
  mutate(
    datetime_utc = as.POSIXct(datetime_utc, tz = "UTC"),
    hour = hour(datetime_utc) + minute(datetime_utc)/60,
    dow  = wday(datetime_utc, label = TRUE, week_start = 1),
    hour_rad = 2*pi*hour/24,
    sin_hour = sin(hour_rad),
    cos_hour = cos(hour_rad)
  )

# 2) Fit OLS model (Model 2)
m2 <- lm(actual ~ forecast + sin_hour + cos_hour + dow, data = df)

# 3) Robust HC1 coefficient table
rob <- lmtest::coeftest(m2, vcov = sandwich::vcovHC(m2, type = "HC1"))

table4_coef <- broom::tidy(rob) %>%
  transmute(
    term = term,
    estimate = round(estimate, 3),
    se_HC1 = round(std.error, 3),
    t = round(statistic, 2),
    p = format.pval(p.value, digits = 3)
  )

# 4) Add model-fit rows (journal style)
table4_fit <- tibble(
  term = c("N", "R-squared", "Adj. R-squared", "Residual SE"),
  estimate = c(
    nobs(m2),
    summary(m2)$r.squared,
    summary(m2)$adj.r.squared,
    summary(m2)$sigma
  ),
  se_HC1 = NA, t = NA, p = NA
) %>%
  mutate(estimate = round(as.numeric(estimate), 3))

table4 <- bind_rows(table4_coef, table4_fit)

# ---- Table 5 (fixed): Combined Index category analysis ----
library(tidyverse)
library(car)
library(rstatix)
library(knitr)

# Safety: remove any old table5 lying around
rm(list = c("table5", "desc", "gh"), envir = .GlobalEnv)

# Ensure index is a factor with sensible order
df <- df %>%
  mutate(index = factor(index, levels = c("low", "moderate", "high", "very high")))

# 1) Descriptive statistics by Index
desc <- df %>%
  group_by(index) %>%
  summarise(
    n = n(),
    mean = mean(actual, na.rm = TRUE),
    sd = sd(actual, na.rm = TRUE),
    median = median(actual, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(across(where(is.numeric), ~round(.x, 2))) %>%
  mutate(section = "Descriptive statistics")

# 2) Levene’s test
lev <- car::leveneTest(actual ~ index, data = df)
levene_row <- tibble(
  section = "Levene’s test (variance equality)",
  index = "All groups",
  n = NA_real_,
  mean = NA_real_,
  sd = NA_real_,
  median = NA_real_,
  detail = paste0(
    "F = ", round(unname(lev[1, "F value"]), 2),
    ", df = ", lev[1, "Df"], ", ", lev[2, "Df"],
    ", p = ", format.pval(unname(lev[1, "Pr(>F)"]), digits = 3)
  )
)

# 3) Welch ANOVA
welch <- oneway.test(actual ~ index, data = df)
welch_row <- tibble(
  section = "Welch one-way ANOVA",
  index = "All groups",
  n = NA_real_,
  mean = NA_real_,
  sd = NA_real_,
  median = NA_real_,
  detail = paste0(
    "F = ", round(unname(welch$statistic), 2),
    ", df_num = ", round(unname(welch$parameter[1]), 2),
    ", df_den = ", round(unname(welch$parameter[2]), 2),
    ", p = ", format.pval(welch$p.value, digits = 3)
  )
)

# 4) Games–Howell post-hoc
gh <- rstatix::games_howell_test(df, actual ~ index) %>%
  mutate(
    index = paste(group1, "vs", group2),
    n = NA_real_,
    mean = round(estimate, 2),                      # mean difference
    sd = NA_real_,                                  # keep numeric column type
    median = NA_real_,
    detail = paste0(
      "95% CI [", round(conf.low, 2), ", ", round(conf.high, 2),
      "], p_adj = ", format.pval(p.adj, digits = 3)
    ),
    section = "Games–Howell post-hoc (pairwise)"
  ) %>%
  select(section, index, n, mean, sd, median, detail)

# Add empty detail column to desc so all pieces align
desc2 <- desc %>%
  mutate(detail = NA_character_) %>%
  select(section, index, n, mean, sd, median, detail)

# 5) Combine into ONE Table 5
table5_combined <- bind_rows(desc2, levene_row, welch_row, gh)

# 6) Print
print(table5_combined)

knitr::kable(
  table5_combined,
  caption = "Table 5. Differences in observed half-hourly grid carbon intensity across Carbon Intensity Index categories. Descriptive statistics are followed by Levene’s test for variance equality, Welch’s one-way ANOVA, and Games–Howell post-hoc pairwise comparisons."
)
# ---- Table 4 (fixed): Combined OLS regression (HC1) + fit stats + error metrics ----
library(tidyverse)
library(lubridate)
library(broom)
library(sandwich)
library(lmtest)
library(knitr)

# Safety: remove any old objects that might confuse things
rm(list = c("table4", "table4_combined", "coef_tbl", "fit_tbl", "err_tbl"), envir = .GlobalEnv)

# 1) Ensure variables exist + create diurnal terms
df <- df %>%
  mutate(
    datetime_utc = as.POSIXct(datetime_utc, tz = "UTC"),
    hour = hour(datetime_utc) + minute(datetime_utc)/60,
    dow  = factor(wday(datetime_utc, label = TRUE, week_start = 1), ordered = FALSE),
    hour_rad = 2*pi*hour/24,
    sin_hour = sin(hour_rad),
    cos_hour = cos(hour_rad),
    error = actual - forecast
  )

# 2) Fit OLS model
m2 <- lm(actual ~ forecast + sin_hour + cos_hour + dow, data = df)

# 3) Robust HC1 coefficient table
rob <- lmtest::coeftest(m2, vcov = sandwich::vcovHC(m2, type = "HC1"))

coef_tbl <- broom::tidy(rob) %>%
  transmute(
    section = "Model coefficients (OLS; HC1 robust SE)",
    term = term,
    estimate = round(estimate, 3),
    se_HC1 = round(std.error, 3),
    t = round(statistic, 2),
    p = format.pval(p.value, digits = 3),
    detail = NA_character_
  )

# 4) Model fit rows (keep same column types)
fit_tbl <- tibble(
  section = "Model fit",
  term = c("N", "R-squared", "Adj. R-squared", "Residual SE (sigma)"),
  estimate = c(
    nobs(m2),
    summary(m2)$r.squared,
    summary(m2)$adj.r.squared,
    summary(m2)$sigma
  ),
  se_HC1 = NA_real_,
  t = NA_real_,
  p = NA_character_,
  detail = NA_character_
) %>%
  mutate(estimate = round(as.numeric(estimate), 3))

# 5) Error metrics (Actual - Forecast) as rows (same column types)
err_vals <- tibble(
  MAE  = mean(abs(df$error), na.rm = TRUE),
  RMSE = sqrt(mean(df$error^2, na.rm = TRUE)),
  Bias = mean(df$error, na.rm = TRUE)
)

err_tbl <- tibble(
  section = "Forecast error metrics (Actual - Forecast)",
  term = c("MAE", "RMSE", "Mean bias"),
  estimate = round(c(err_vals$MAE, err_vals$RMSE, err_vals$Bias), 3),
  se_HC1 = NA_real_,
  t = NA_real_,
  p = NA_character_,
  detail = NA_character_
)

# 6) Combine into ONE Table 4
table4_combined <- bind_rows(coef_tbl, fit_tbl, err_tbl)

# 7) Print
print(table4_combined)

knitr::kable(
  table4_combined,
  caption = "Table 4. OLS regression of observed half-hourly grid carbon intensity on forecast carbon intensity and temporal structure (diurnal sine/cosine terms and day-of-week). Coefficients are reported with heteroskedasticity-robust (HC1) standard errors; model fit statistics and forecast error metrics (Actual − Forecast) are also shown."
)

 