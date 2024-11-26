library("ggplot2")
data <- read.table("Multistate_mcmc_samples_5M.txt", header=TRUE, sep="\t")
mean_data <- mean(data$qMN)
sd_data <- sd(data$qMN)

# Calculate confidence interval
conf_int <- mean_data + c(-1.96, 1.96) * sd_data / sqrt(length(data$qMN))

# Create a histogram with ggplot2
ggplot(data = data, aes(x = qMN)) + 
    geom_histogram(binwidth = 0.001, fill = "black", color = "black") + 
    geom_vline(xintercept = mean_data, color = "red", linetype = "dashed") + 
    # Highlight the CI range using geom_rect instead of geom_ribbon
    geom_rect(aes(xmin = conf_int[1], xmax = conf_int[2], ymin = 0, ymax = Inf),
              fill = "blue", alpha = 0.2) +  # Blue shaded CI area
    ggtitle("Histogram of Motile to Non-Motile Transition Rates (qMN) with 95% CI") + 
    xlab("qMN") + 
    ylab("Frequency") + 
    scale_x_continuous(breaks = c(0, 0.2, 0.4, 0.6, 1.0), limits = c(0, 1)) +  # Set limits from 0 to 1
    theme(
        panel.background = element_blank(),  # Remove background color
        plot.background = element_blank(),   # Remove plot background
        panel.grid = element_blank(),        # Remove grid lines
        axis.line = element_line(color = "black")  # Keep axis lines
    ) 
