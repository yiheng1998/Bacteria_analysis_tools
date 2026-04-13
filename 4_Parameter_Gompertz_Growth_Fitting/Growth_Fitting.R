# Load required R packages
library(drc)    
library(ggplot2) 

# ------------------------------------------------------------------------------
# Universal Microbial Growth Curve Fitting Function
# ------------------------------------------------------------------------------
# Parameters:
#   file_path    : Full path to CSV file
#   time_col     : Column name/index of time points (e.g., "Time" or 1)
#   od_start_col : First column name/index of OD replicates (e.g., "Rep1" or 2)
#   pred_step    : Step size for fitted curve (default = 0.5)
# ------------------------------------------------------------------------------
fit_growth_curve_csv <- function(file_path,
                                 time_col,
                                 od_start_col,
                                 pred_step = 0.5) {
  
  # Step 1: Read CSV data
  raw_data <- read.csv(file_path, stringsAsFactors = FALSE, strip.white = TRUE)
  
  # Clean column names (remove extra spaces)
  colnames(raw_data) <- trimws(colnames(raw_data))
  
  # Validate columns
  if (!time_col %in% colnames(raw_data)) {
    stop(sprintf("ERROR: Time column '%s' not found! Available columns: %s", 
                 time_col, paste(colnames(raw_data), collapse=", ")))
  }
  if (!od_start_col %in% colnames(raw_data)) {
    stop(sprintf("ERROR: OD start column '%s' not found! Available columns: %s", 
                 od_start_col, paste(colnames(raw_data), collapse=", ")))
  }
  
  # Extract data
  time_vec <- raw_data[[time_col]]
  od_cols <- colnames(raw_data)[which(colnames(raw_data) == od_start_col):ncol(raw_data)]
  od_data <- raw_data[, od_cols, drop = FALSE]
  
  # Step 2: Calculate mean and SD
  processed_df <- data.frame(
    time = time_vec,
    mean_od = rowMeans(od_data, na.rm = TRUE),
    sd_od = apply(od_data, 1, function(x) sd(x, na.rm = TRUE))
  )
  
  # Step 3: Gompertz model fitting
  gompertz_fit <- drm(mean_od ~ time,
                      data = processed_df,
                      fct = G.4(),
                      na.action = na.omit)
  
  # Step 4: Prediction data
  max_time <- max(time_vec)
  pred_time <- data.frame(time = seq(min(time_vec), max_time, by = pred_step))
  pred_od <- predict(gompertz_fit, newdata = pred_time)
  
  # Step 5: Growth phase division
  max_pred <- max(pred_od)
  lag_end <- pred_time$time[which.min(abs(pred_od - 0.37 * max_pred))]
  log_end <- pred_time$time[which.min(abs(pred_od - 0.90 * max_pred))]
  
  # Step 6: Adaptive plot
  growth_plot <- ggplot(processed_df, aes(x = time, y = mean_od)) +
    annotate("rect", xmin = min(time_vec), xmax = lag_end,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "gray50") +
    annotate("rect", xmin = lag_end, xmax = log_end,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "deepskyblue2") +
    annotate("rect", xmin = log_end, xmax = max_time,
             ymin = -Inf, ymax = Inf, alpha = 0.1, fill = "green4") +
    geom_errorbar(aes(ymin = mean_od - sd_od, ymax = mean_od + sd_od),
                  width = (max(time_vec)-min(time_vec))/50,
                  color = "#2E5A99", linewidth = 0.8) +
    geom_point(shape = 16, size = 3.5, color = "#2E5A99") +
    geom_line(data = pred_time, aes(y = pred_od),
              linewidth = 1.2, color = "#2E5A99") +
    annotate("text", x = (min(time_vec)+lag_end)/2, y = max(processed_df$mean_od)*0.9,
             label = "Lag Phase", fontface = 2, size = 4.5) +
    annotate("text", x = (lag_end+log_end)/2, y = max(processed_df$mean_od)*0.9,
             label = "Log Phase", fontface = 2, size = 4.5) +
    annotate("text", x = (log_end+max_time)/2, y = max(processed_df$mean_od)*0.9,
             label = "Stationary Phase", fontface = 2, size = 4.5) +
    labs(x = "Time (h)", y = "OD600") +
    theme_bw() +
    theme(panel.grid = element_line(color = "#E5E5E5"),
          panel.border = element_rect(linewidth = 0.8),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 13, face = "bold"))
  
  # Step 7: Print results
  cat("=============================================\n")
  cat("4-Parameter Gompertz Model Fitting Results\n")
  cat("=============================================\n")
  print(coef(gompertz_fit))
  cat("=============================================\n")
  
  return(list(processed_data = processed_df,
              fitted_model = gompertz_fit,
              growth_plot = growth_plot))
}
