require(ggplot2)
require(tidyverse)
data <- read.csv("LabK_17.csv") %>% select(contains(c("Current", "Intensity", "Angle", "Time"))) 
PI = 3.141579
A = 1.684 
B = 15.2e3
T0 = 293
RT0 = 0.93
alpha0 = 5.4e-3
m2nm_conv = 1e9

# Tidy dataset
full_data <- data %>% 
            pivot_longer(cols=everything()) %>% 
            mutate(Run = as.integer(str_extract(name, "(\\d+$)"))) %>%
            mutate(name = str_extract(name, "^(.*?)(?=Run\\.\\.\\d+$)")) %>%
            pivot_wider(id_expand = T, names_from = name, values_from = value, values_fn = list) %>%
            unnest(!contains("Run")) %>%
            drop_na() %>%
            rename("Current" = Output.Current.Ch.H..A..,
                   "Angle"= Angle..rad..,
                   "Intensity" = Relative.Intensity.....,
                   "Time" = Time..s..)


# Examining other data I found range 2-2.25 for our spectrumlight we find around 1.5-1.68 thus we 
deg_per_rev = 7.152770364
start_deg = 67

processed <- full_data %>%
  filter(Angle <= 25) %>%
  group_by(Run) %>%
  # Calculate Theta
  mutate(Rot_deg = start_deg-(Angle/(2*PI))*deg_per_rev) %>%
  mutate(Theta = Rot_deg*(PI/180)) %>%

  # Calculate refraction 
  mutate(refraction = sqrt(((2 / sqrt(3)) * sin(Theta) + 0.5)^2 + 0.75)) %>%
  # Filter out asymptotic refraction values
  # filter(refraction < 1.684) %>%
  
  # Calculate Lambda
  mutate(lambda = sqrt((B/(refraction - A))) / m2nm_conv)

# Tune A
while (any(is.na(processed$lambda))) {
  A = A - 0.01
  processed$lambda = sqrt((B / (processed$refraction - A)))  / m2nm_conv
}

write.csv(processed, file="processeddata.csv")

ggplot(processed, aes(x=Time, y=lambda, col=Run, group=Run)) + geom_line()


# Plot curves
ggplot(processed %>% filter(lambda < 1e-6), aes(x=lambda, y=Intensity, col=Run, group=Run)) + geom_line()
ggsave("curves.png")

# Get Resistance data
voltage_data <- read.csv("voltagedata.csv") %>%
                rename("LampVoltage" = Vlamp..V.,
                       "SetVoltage" = Vset..V.)
current_data <- processed %>%
                group_by(Run) %>%
                filter(Current > 0.67) %>%
                summarise(Current = max(Current))
resistance_data <- voltage_data %>%
                   left_join(current_data, by = "Run") %>%
                   mutate(resistance = LampVoltage/Current)
relevant_data <- processed %>% select(c("Run", "lambda", "Intensity"))

lambda_T <- relevant_data %>%
        # Filters peak intensity
            group_by(Run) %>%
            filter(Intensity == max(Intensity)) %>%
        # Calculate Temperature
            left_join(resistance_data, by = "Run") %>%
            select(!c("Current", "LampVoltage")) %>%
            1.
            mutate(temperature = T0 + (((resistance/RT0)-1)/alpha0)) %>%
            group_by(temperature) %>%
        # Get Mean and SD
            summarise(lambda_peak = mean(lambda), lambda_peak_error = sd(lambda)) %>%
            drop_na() %>%
        # Get 1/T
            mutate(Tinv = 1/temperature)


write.csv(lambda_T, file="peaklambda_temp.csv")


## Calculate k
lTinv_model = lm(lambda_peak ~ Tinv, data=lambda_T)
summary(lTinv_model)

slope = as.numeric(lTinv_model$coefficients[2])

constant0 = (1/4.965114)
h = 6.6e-34
c = 2.98e8

k_b = (constant0*h*c)/slope

ggplot(lambda_T, aes(x = Tinv, y = lambda_peak)) + 
  geom_point() + 
  geom_line(aes(x=Tinv, y=lTinv_model$fitted.values)) +
  geom_errorbar(aes(ymin = lambda_peak - lambda_peak_error, 
                    ymax = lambda_peak + lambda_peak_error)) +
  labs(x = "1/Temperature (K)", y = "Mean Peak Wavelength (m)")
ggsave("lambdatemp.png")

ggplot(lambda_T, aes(x=temperature, y=lambda_peak)) + geom_point()

## Total Intensity

auc <- function(x, y) {
  sum_val = 0
  for (i in 1:(length(x) - 1)) {
    sum_val = sum_val + (y[i] + 0.5 * abs(y[i] - y[i+1])) * abs(x[i+1] - x[i])
  }
  return(sum_val)
}

total_intensity_temp <- relevant_data %>%
        group_by(Run) %>%
        # Calculate Temperature
        left_join(resistance_data, by = "Run") %>%
        select(!c("Current", "LampVoltage")) %>%

        mutate(temperature = T0 + (((resistance/RT0)-1)/alpha0)) %>%
        # Final summary
  # Filter outliers
        group_by(temperature) %>%
        summarise(U = auc(lambda, Intensity)) %>%
        mutate(temperature4 = temperature^4) %>%
        filter(U > 4e-6) %>%
        mutate(U_error = U*0.1)
write.csv(total_intensity_temp, "Utemp.csv")

# Model relation
UT_model = lm(U~temperature4, data=total_intensity_temp)
summary(UT_model)
sigma_approx = as.numeric(UT_model$coefficients[2])

sigma = 5.670374419e-8

dev = abs(sigma_approx - sigma)

ggplot(total_intensity_temp, aes(x=temperature4, y=U)) + 
  geom_point() + 
  geom_errorbar(aes(ymin=U-U_error,ymax=U+U_error)) + 
  geom_line(aes(y=UT_model$fitted.values)) +
  labs(x="Temperature^4 (K^4)", y="Total Relative Intensity (%)")
ggsave("UT4.png")

k_b_true = 1.380649e-23