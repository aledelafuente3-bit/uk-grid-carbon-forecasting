# UK Electricity-Grid Carbon-Intensity Forecasting

Time-series analysis of UK electricity-grid carbon intensity, developed to test how trend and seasonality affect short-term forecasts of power-system emissions.

## Analytical workflow

- cleaned and explored the carbon-intensity series in R;
- tested temporal structure with autocorrelation and partial-autocorrelation diagnostics;
- assessed seasonal patterns and stationarity;
- developed and compared forecasting specifications, including ARIMA-based approaches; and
- evaluated model performance using held-out observations and diagnostic plots.

## Repository contents

```text
code/LASTSCRIPT.R              Main analysis and forecasting workflow
data/Carbon_Intensity_Data.csv  Analytical dataset
figures/                        Diagnostic and model-output figures
```

## Relevance

Grid-carbon forecasts can support operational scheduling, energy procurement and location-specific assessment of electricity demand. The project demonstrates statistical model selection, diagnostics and translation of power-system data into decision-useful evidence.

## Tools

R · time-series analysis · ARIMA · forecast evaluation · data visualisation
