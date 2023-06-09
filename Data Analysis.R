
#Library
library(ggplot2)
library(ggpubr)
library(dplyr)
library(googlesheets4)
library(lme4)
library(lmerTest)
library(tidyr)
library(car)
library(ggfortify)
library(MuMIn)
library(emmeans)
library(rstatix)


#Depth interpolation calculations!!

# Load the data from google sheets
data_measurements = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "measurements")

## Depth calculations

# Create a hypotenuse
data_measurements$transect_pos_cm = data_measurements$knot*10-10


# Add missing depths
data_depth = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "depth")


# Here the depth is recalculated 
data_depth$depth_corrected = data_depth$`viva(cm)` + data_depth$`depth_measured (cm)`

# Here the depth_corrected is removed when there is no measurement
data_depth$depth_corrected[is.na(data_depth$`depth_measured (cm)`)]=NA
data_depth$transect_pos_cm = data_depth$knot*10-10


# Step 1: Interpolation within known points
data_depth <- data_depth %>%
  group_by(transectlabel) %>%
  arrange(transect_pos_cm) %>%
  mutate(depth_interpolated = approx(transect_pos_cm, depth_corrected, transect_pos_cm, rule = 1)$y)

# Step 1.5: Calculate slope between each pair of known points
data_depth <- data_depth %>%
  mutate(slope_interpolated = (depth_interpolated - lag(depth_interpolated)) / (transect_pos_cm - lag(transect_pos_cm)))

# Step 2: Extrapolation for the edges
data_depth <- data_depth %>%
  ungroup() %>%
  arrange(transectlabel, transect_pos_cm)

# For-loop for start points
for(i in unique(data_depth$transectlabel)) {
  # Calculate slope of the first non-NA pair
  slope_start <- data_depth$slope_interpolated[data_depth$transectlabel == i & !is.na(data_depth$slope_interpolated)][1]
  
  # Find the indices to be replaced
  indices <- which(data_depth$transectlabel == i & is.na(data_depth$depth_interpolated) & data_depth$transect_pos_cm < min(data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)]))
  
  # Replace the values
  data_depth$depth_interpolated[indices] <- data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][1] - slope_start * (data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][1] - data_depth$transect_pos_cm[indices])
}

# For-loop for end points
for(i in unique(data_depth$transectlabel)) {
  # Calculate slope of the last non-NA pair
  slope_end <- data_depth$slope_interpolated[data_depth$transectlabel == i & !is.na(data_depth$slope_interpolated)][length(data_depth$slope_interpolated[data_depth$transectlabel == i & !is.na(data_depth$slope_interpolated)])]
  
  # Find the indices to be replaced
  indices <- which(data_depth$transectlabel == i & is.na(data_depth$depth_interpolated) & data_depth$transect_pos_cm > max(data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)]))
  
  # Replace the values
  data_depth$depth_interpolated[indices] <- data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][length(data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)])] + slope_end * (data_depth$transect_pos_cm[indices] - data_depth$transect_pos_cm[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)][length(data_depth$depth_interpolated[data_depth$transectlabel == i & !is.na(data_depth$depth_interpolated)])])
}

# In these for-loops, we iterate over each unique transect. For each transect, we first compute the slope at the start (or end) using the first (or last) non-NA interpolated depth value. We then identify the indices of the depth_interpolated values to be replaced, and replace them with the extrapol

# Join the depth with the original data and round it to 2 decimals

data_measurements = left_join(data_measurements,data_depth, by =c("transectlabel","transect_pos_cm"))
data_measurements$depth_interpolated = round(data_measurements$depth_interpolated, 2)


######


#Make the interpolated depth positive

data_measurements$depth_interpolated = data_measurements$depth_interpolated*-1


#### Adult and germling distribution data!!
####Multi density chart
tra_dat = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "adult", col_types = "ccccccccnnncnncc")

# Adult depth data calculations
# Subset the data needed to interpolate the depth to the different points
# Correct the depth by the water level from the nearby station
depth_data <- 
  tra_dat %>% 
  mutate(depth_correct = (water_level_cm + depth_cm) ) %>%
  select(date, transect_id, position, depth_correct) %>%
  distinct()

# Problems with the start of transect 2, remove positions 0 to 4
# We do not have a starting depth
depth_data <- 
  depth_data %>%
  filter( !(transect_id == 2 & position %in% 0:4) )

# Split into a list
depth_list <- split(depth_data, depth_data$transect_id)

# Loop over all transects
depth_out <- vector("list", length = length(depth_list))
for(i in 1:length(depth_list)) {
  
  # Initialise a data.frame to work with
  df <- depth_list[[i]]
  
  # Get the dividers
  dividers <- which(!is.na(df$depth_correct) )
  
  # Duplicate middles for which the data are needed for multiple calculations
  if (length(dividers) > 2) {
    
    dups <- dividers[-c(1, length(dividers))]
    
    df <- 
      df[c(1:nrow(df), dups), ] %>%
      arrange(date, transect_id, position)
    
  }
  
  # Add an ID column
  x <- vector("list", length = (length(dividers)-1))
  for(j in 1:(length(dividers)-1) ) {
    
    x[[j]] <- rep(j, ( (dividers[j+1] - dividers[j])+1 ) )
    
  }
  
  # Write the ID column into the data.frame
  df$position_id <- unlist(x)
  
  # Split data by the position ID
  df_list <- split(df, df$position_id)
  
  # For each block, interpolate the depth between the points
  int_depth <- 
    lapply(df_list, function(y){
      
      y %>%
        group_by(position_id) %>%
        mutate(y_int = first(depth_correct),
               d_depth = (last(depth_correct) - first(depth_correct)),
               d_position = last(10*position) - first(10*position)) %>%
        ungroup() %>%
        mutate(slope = d_depth/d_position) %>%
        select(-d_depth, -d_position) %>%
        mutate(distance = 0:(length(position)-1)*10 ) %>%
        mutate(depth_interpolated = y_int + ((distance)*slope) )  %>%
        select(-distance) %>%
        mutate(depth_interpolated = if_else(!is.na(depth_correct), depth_correct, depth_interpolated)) %>%
        select(-y_int, -slope, -position_id)
      
    }) %>%
    bind_rows(.,) %>%
    distinct()
  
  depth_out[[i]] <- int_depth
  
} 

# Bind the list into a data.frame
depth_out <- bind_rows(depth_out)

# Join this back to the full dataset
tra_a <- 
  full_join(tra_dat, depth_out, by = c("date", "transect_id", "position")) %>%
  select(date, transect_id, site_code, time, position, water_level_cm, 
         depth_cm, depth_correct, depth_interpolated, binomial_code, length_cm, circum_cm, 
         field_observer, notes)

# Check the summary statistics
summary(tra_a)

# Check for unique values, especially for the binomial codes
lapply(tra_a, function(x) unique(x))

# Remove missing binomial codes and the missing values "", NA
tra_a <- 
  tra_a %>%
  filter( !(binomial_code %in% c("-9999", "") | is.na(binomial_code) | is.na(depth_interpolated) )  )

# Work with the transect data
ggplot(data = tra_a,
       mapping = aes(x = binomial_code, y = depth_interpolated)) +
  geom_point()

transect_summary <- 
  tra_a %>% 
  group_by(binomial_code) %>%
  summarise(mean_depth = mean(depth_interpolated, na.rm = TRUE),
            min_depth = min(depth_interpolated, na.rm = TRUE),
            max_depth = max(depth_interpolated, na.rm = TRUE),
            quant_20 = quantile(depth_interpolated, 0.20, na.rm = TRUE),
            quant_80 = quantile(depth_interpolated, 0.80, na.rm = TRUE))

transect_summary$data_id <- "transect"

transect_summary <- 
  transect_summary %>%
  select(data_id, binomial_code, mean_depth:quant_80)

# Load the raw allometric data
allo_dat = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "adult2", col_types = "ccccdcccccdcddddccdddccdccc")


# Remove the missing values from allo_dat
allo_dat <- 
  allo_dat %>%
  filter(!is.na(depth_cm)) %>%
  filter(!is.na(water_level_cm)) %>%
  mutate(depth_correct = depth_cm + water_level_cm)

nrow(allo_dat)
head(allo_dat)
names(allo_dat)

allo_dat %>%
  group_by(binomial_code) %>%
  summarise(n = n())

unique(allo_dat$site_code)

# Plot out the depth distribution
allo_dat %>%
  ggplot(data = ., 
         mapping = aes(x = binomial_code, y = depth_correct)) +
  geom_point()

# Get summary statistics
allo_summary <- 
  allo_dat %>%
  group_by(binomial_code) %>%
  summarise(mean_depth = mean(depth_correct, na.rm = TRUE),
            min_depth = min(depth_correct, na.rm = TRUE),
            max_depth = max(depth_correct, na.rm = TRUE),
            quant_20 = quantile(depth_correct, 0.20, na.rm = TRUE),
            quant_80 = quantile(depth_correct, 0.80, na.rm = TRUE))

# Add an allometry identifier
allo_summary$data_id <- "allometry_data"

allo_summary <- 
  allo_summary %>%
  select(data_id, binomial_code, mean_depth:quant_80)

# We used these summary statistics to choose our experimental depths

# Plot the depth data for each species
names(allo_dat)
names(tra_a)

# Bind these data together
all_depth <- 
  bind_rows(
    
    select(allo_dat, binomial_code, depth_correct ) %>%
      mutate(data_set = "allometry"),
    
    select(tra_a, binomial_code, depth_correct = depth_interpolated) %>%
      mutate(data_set = "transect")
    
  )


all_depth = all_depth %>% filter(binomial_code == "fu_ve" |binomial_code == "fu_se" |binomial_code == "fu_sp" |binomial_code == "as_no" )

all_depth$serrated = "no"
all_depth$serrated[all_depth$binomial_code == "fu_se"] = "yes"
all_depth$serrated[all_depth$binomial_code == "as_no"] = NA

#Plot of all 4 species
#Plot of morphology of germlings and adults

# Filter the relevant data
adult_data <- all_depth %>% 
  filter(!is.na(binomial_code) & !is.na(depth_correct) & data_set == "transect")

# Create a named vector for custom legend labels
legend_labels <- c("fu_ve" = "F. vesiculosus",
                   "fu_se" = "F. serratus",
                   "fu_sp" = "F. spiralis",
                   "as_no" = "A. nodosum")

# Create the multi-density plot with positive x-axis values
chart = ggplot(adult_data, aes(x = abs(depth_correct), fill = binomial_code)) +
  geom_density(alpha = 0.5) +
  labs(title = "",
       x = "Depth (cm)",
       y = "Density",
       fill = "Species") +
  scale_fill_manual(values = c("fu_sp" = "#99FF33",
                               "fu_ve" = "#6699FF",
                               "as_no" = "#FF99CC",
                               "fu_se" = "#FF6666"),
                    labels = legend_labels) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) +
  theme(legend.position = "bottom") +
  scale_x_continuous(labels = abs) +
  xlim(c(0,85))

print(chart)

adult_data$serrated[is.na(adult_data$serrated)] = "A. nodosum"

#Change the legends labels
legend_labels <- c("no" = "Not serrated",
                   "yes" = "Serrated",
                   "A. nodosum" = "A. nodosum")

#Adult morphology distribution - flip and reverse and stuff and add Aschophyllum
p_adults = ggplot(adult_data, aes(x = abs(depth_correct), fill = serrated)) +
  geom_density(alpha = 0.7) +
  xlab("Depth (cm)") +
  ylab("Density") +
  scale_fill_manual(values = c("#FF99CC", "#6699FF", "#FF6666"), drop = TRUE, labels = legend_labels) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  guides(fill = guide_legend(title = "Morphology",
                             label.theme = element_text(size = 20))) + 
  xlim(c(0,85))


#Exposed and Sheltered Germlings 

#Sheltered germlings
# Remove NA values from Morphology variable
data_sheltered = data_measurements %>% filter(Condition != "Exposed")
data_sheltered <- data_sheltered[!is.na(data_sheltered$Morphology), ]

p_germling_sheltered = ggplot(data_sheltered, aes(x = depth_interpolated, fill = Morphology)) +
  geom_density(alpha = 0.7) +
  xlab("Depth (cm)") +
  ylab("Density") +
  scale_fill_manual(values = c("#FF99CC", "#6699FF", "#FF6666"), drop = TRUE) +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  guides(fill = guide_legend(title = "Morphology",
                             label.theme = element_text(size = 20))) + 
  xlim(c(0,85))

#Exposed germlings - NA there?

#data_exposed = data_measurements %>% filter(Condition != "Sheltered")
#data_exposed <- data_exposed[!is.na(data_exposed$Condition), ]

data_exposed = data_measurements %>% 
  filter(Condition != "Sheltered" & !is.na(depth_interpolated) & !is.na(Morphology))

p_germling_exposed = ggplot(data_exposed, aes(x = depth_interpolated, fill = Morphology)) +
  geom_density(alpha = 0.7, na.rm = TRUE) +
  xlab("Depth (cm)") +
  ylab("Density") +
  scale_fill_manual(values = c("#6699FF", "#FF6666"), drop = TRUE) + theme_classic() +
  theme_classic() +
  theme(axis.text = element_text(size = 20),
        axis.title = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20)) + 
  guides(fill = guide_legend(title = "Morphology",
                             label.theme = element_text(size = 20))) + 
  xlim(c(0,85))

density_plots = ggarrange(chart, p_adults, p_germling_sheltered, p_germling_exposed,
                          nrow =4, ncol = 1,
                          legend = "none",
                          labels= c("A","B","C","D")) +
  theme(text = element_text(size = 16))

ggsave("C:/Users/hugom/OneDrive/Skrivbord/examensarbete/Graphs/density_plots.png",
       plot = density_plots, width = 20, height = 14, dpi = 300)


####Germling traits!!!

# Filter away A. nodosum - distinguishable as germlings

data_fucus = data_measurements %>% filter(Morphology != "A. nodosum")

#####Do traits of fucus change on a depth gradient

#TDMC
m1 = lm(TDMC ~ depth_interpolated, data = data_fucus)
summary(m1)
m1_length = lmer(TDMC ~ depth_interpolated + length_mm + (1 | transectlabel), data = data_fucus)
summary(m1_length)
#TDMC does not depend on depth when you control for length, therefore TDMC depends on size

#TDMC plot
plot_tdmc <- ggscatter(x = "depth_interpolated", y = "TDMC",
                       xlab = "Depth (cm)", ylab = "TDMC",
                       color = "Morphology",
                       palette = c("#3399FF", "#FF6666"),
                       data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m1$coefficients[1], slope = m1$coefficients[2], size = 1) +
  guides(color = guide_legend(title = "Morphology",
                              label.theme = element_text(size = 15)))


#SA:P
m2 = lm(SA_P ~ depth_interpolated, data = data_fucus)
summary(m2)

m2_length = lmer(SA_P ~ depth_interpolated + length_mm + (1 | transectlabel), data = data_fucus)
summary(m2_length)

#SA:P does not depend on depth when you control for length, therefore SA:P depends on length

#SA:P plot
plot_sap <- ggscatter(x = "depth_interpolated", y = "SA_P",
                      xlab = "Depth (cm)", ylab = "SA:P",
                      color = "Morphology",
                      palette = c("#3399FF", "#FF6666"),
                      data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m2$coefficients[1], slope = m2$coefficients[2], size = 1) +
  guides(color = FALSE)


#LW ratio
m3 = lm(LW_ratio ~ depth_interpolated, data = data_fucus)
summary(m3)
m3_length = lmer(LW_ratio ~ depth_interpolated + length_mm + (1 | transectlabel), data = data_fucus)
summary(m3_length)
#LW-ratio still depends on depth when you control for length

#LW ratio plot
plot_lw <- ggscatter(x = "depth_interpolated", y = "LW_ratio",
                     xlab = "Depth (cm)", ylab = "Length width ratio",
                     color = "Morphology",
                     palette = c("#3399FF", "#FF6666"),
                     data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m3$coefficients[1], slope = m3$coefficients[2], size = 1) +
  guides(color = FALSE)

#LP ratio
m4 = lm(LP_ratio ~ depth_interpolated, data = data_fucus)
summary(m4)
m4_length = lmer(LP_ratio ~ depth_interpolated + length_mm + (1 | transectlabel), data = data_fucus)
summary(m4_length)
#LP-ratio depends more on depth when you control for length, but also on depth

#LP ratio plot
plot_lp <- ggscatter(x = "depth_interpolated", y = "LP_ratio",
                     xlab = "Depth (cm)", ylab = "Length perimeter ratio",
                     color = "Morphology",
                     palette = c("#3399FF", "#FF6666"),
                     data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m4$coefficients[1], slope = m4$coefficients[2], size = 1) +
  guides(color = FALSE)


#STA(mm2g-1)
m5 = lm(STAmm2_g ~ depth_interpolated, data = data_fucus)
summary(m5)
m5_length = lmer(STAmm2_g ~ depth_interpolated + length_mm + (1 | transectlabel), data = data_fucus)
summary(m5_length)
#STA still depends on depth when controlling for length, which it also slightly depends on

#STA plot
plot_sta <- ggscatter(x = "depth_interpolated", y = "STAmm2_g",
                     xlab = "Depth (cm)",
                     color = "Morphology",
                     palette = c("#3399FF", "#FF6666"),
                     data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m5$coefficients[1], slope = m5$coefficients[2], size = 1) +
  guides(color = FALSE) +
  labs(y = bquote('STA' ~mm^2 ~g^-1))

#Scatter plots - arranged
scatter_plots = ggarrange(plot_tdmc, plot_sap, plot_lw, plot_lp, plot_sta,
                          nrow =2, ncol = 3, common.legend = TRUE,
                          legend = "top",
                          labels= c("A","B","C","D","E")) +
  theme(text = element_text(size = 50))

#Save to desktop
ggsave("C:/Users/hugom/OneDrive/Skrivbord/examensarbete/Graphs/scatter_plots.png",
       plot = scatter_plots, width = 20, height = 14,
       dpi = 300)


#Anova tests
m12 = lmer(TDMC ~ Morphology * Condition + length_mm + (1 | transect), data = data_fucus)
anova(m12)
m13 = lmer(SA_P ~ Morphology * Condition + length_mm + (1 | transect), data = data_fucus)
anova(m13)
m14 = lmer(LW_ratio ~ Morphology * Condition + length_mm + (1 | transect), data = data_fucus)
anova(m14)
m15 = lmer(LP_ratio ~ Morphology * Condition + length_mm + (1 | transect), data = data_fucus)
anova(m15)
m16 = lmer(STAmm2_g ~ Morphology * Condition + length_mm + (1 | transect), data = data_fucus)
anova(m16)


#summarize and get R2 values
library(jtools)
summ(m1_length)
summ(m2_length)
summ(m3_length)
summ(m4_length)
summ(m5_length)

#Pearson correlation
data_fucus %>% select(STAmm2_g,TDMC,LW_ratio,LP_ratio,SA_P,depth_interpolated,length_mm) %>% cor_mat()
m17 = data_fucus %>% select(TDMC, SA_P, LW_ratio, LP_ratio, STAmm2_g, depth_interpolated, length_mm) %>% cor_mat() %>% cor_get_pval()
m17

####Is there a difference between exposed and sheltered
#Boxplots

#TDMC
m6 = lm(TDMC ~ Condition * Morphology, data = data_fucus)
summary(m6)
m6_length = lmer(TDMC ~ Condition * Morphology + length_mm + (1 | transectlabel), data = data_fucus)
summary(m6_length)

#TDMC showed a dependent on length alone, but not condition or Morphology
emmeans(m6_length, pairwise~Condition*Morphology)
#TDMC did not differ between condition or morphology

#TDMC plot
box_tdmc = ggboxplot(x = "Condition", y = "TDMC", color = "Morphology", 
          xlab = "Condition", ylab = "TDMC",
          palette = c("#FF6666", "#0099FF"),
          data = data_fucus)

#SA:P
m7 = lm(SA_P ~ Condition * Morphology, data = data_fucus)
summary(m7)
m7_length = lmer(SA_P ~ Condition * Morphology + length_mm + (1 | transectlabel), data = data_fucus)
summary(m7_length)
#SA:P showed a dependent on length alone, but not condition or Morphology
emmeans(m7_length, pairwise~Condition*Morphology)
#SA:P showed difference between Morphology on sheltered shores and morphology on sheltered serrated vs exposed non serrated

#SAP plot
box_sap = ggboxplot(x = "Condition", y = "SA_P", color = "Morphology", 
          xlab = "Condition", ylab = "SA:P",
          palette = c("#FF6666", "#0099FF"),
          data = data_fucus)
          #SAP showed no difference between the Conditions and or Morphologies

#LW-ratio
m8 = lm(LW_ratio ~ Condition * Morphology, data = data_fucus)
summary(m8)
m8_length = lmer(LW_ratio ~ Condition * Morphology + length_mm + (1 | transectlabel), data = data_fucus)
summary(m8_length)
#LW-ratio showed a slight dependent on Morphology
emmeans(m8_length, pairwise~Condition*Morphology)
#LW-ratio showed difference between not serrated sheltered and serrated exposed

#LW-ratio plot
box_lw = ggboxplot(x = "Condition", y = "LW_ratio", color = "Morphology", 
          xlab = "Condition", ylab = "Length width ratio",
          palette = c("#FF6666", "#0099FF"),
          data = data_fucus)

#LP-ratio
m9 = lm(LP_ratio ~ Condition * Morphology, data = data_fucus)
summary(m9)
m9_length = lmer(LP_ratio ~ Condition * Morphology + length_mm + (1 | transectlabel), data = data_fucus)
summary(m9_length)
#LP-ratio showed dependent on morphology and length
emmeans(m9_length, pairwise~Condition*Morphology)
#LP-ratio showed difference between not serrated sheltered and serrated exposed germlings

#LP_ratio plot
box_lp = ggboxplot(x = "Condition", y = "LP_ratio", color = "Morphology", 
          xlab = "Condition", ylab = "Length perimeter ratio",
          palette = c("#FF6666", "#0099FF"),
          data = data_fucus)

#STA
m10 = lm(STAmm2_g ~ Condition * Morphology, data = data_fucus)
summary(m10)
m10_length = lmer(STAmm2_g ~ Condition * Morphology + length_mm + (1 | transectlabel), data = data_fucus)
summary(m10_length)
#STA depend on morphology but also slightly on length
emmeans(m10_length, pairwise~Condition*Morphology)
#STA differ mostly from morphology but also a little depending on condition

#STA plot
box_sta = ggboxplot(x = "Condition", y = "STAmm2_g", color = "Morphology", 
          xlab = "Condition",
          palette = c("#FF6666", "#0099FF"),
          data = data_fucus) + 
  labs(y = bquote('STA' ~mm^2 ~g^-1))
          #STA showed great difference between Conditions, not Morphology

#Box plots arranged

box_plots = ggarrange(box_tdmc, box_sap, box_lw, box_lp, box_sta,
                       nrow =2, ncol = 3, common.legend = TRUE,
                       legend = "top",
                       labels= c("A","B","C","D","E")) +
  theme(text = element_text(size = 50))

ggsave("C:/Users/hugom/OneDrive/Skrivbord/examensarbete/Graphs/box_plots.png",
       plot = box_plots, width = 20, height = 14, dpi = 300)


#PCA not used - future project????
m11 = lmer(TDMC ~ scale(length_mm) + scale(depth_interpolated) + (1 | transectlabel), data = data_fucus)
summary(m11)
qqp(resid(m11))

#Select the data for the pca
data_pca = data_fucus %>% select(TDMC, SA_P, LW_ratio, LP_ratio, STAmm2_g, wetweight_g, depth_interpolated, Condition)
data_pca = na.omit(data_pca)
data_pca_result = prcomp(data_pca[c(-6, -7, -8)], scale = TRUE)

#autoplot for pca
data_pca_x = data.frame(data_pca_result$x)

ggscatter(x = "PC1", y = "PC2", color = "PC3", add = "reg.line", data = data_pca_x)

autoplot(data_pca_result, loadings = TRUE, loadings.label = TRUE, color = "Condition", size = "depth_interpolated", data = data_pca)

