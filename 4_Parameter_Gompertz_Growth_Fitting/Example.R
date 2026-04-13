# Run growth curve fitting
result <- fit_growth_curve_csv(
  file_path = "growth_data.csv",  # Full CSV path
  time_col = "Time",              # Exact time column name
  od_start_col = "Rep1"           # Exact first OD column name
)

# Show the plot
print(result$growth_plot)

# View model summary
summary(result$fitted_model)
