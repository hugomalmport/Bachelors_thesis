
#Library
library(ggplot2)
library(ggpubr)
library(dplyr)
library(googlesheets4)
library(lme4)
library(lmerTest)
library(directlabels)
library(patchwork)


data_measurements = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "measurements")
data_measurements$transect_pos_cm = data_measurements$knot*10-10


#Add missing depths
data_depth = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "depth")


# here i recalculate the depth 
data_depth$depth_corrected = data_depth$`viva(cm)` + data_depth$`depth_measured (cm)`

# here i remove the depth corrected when there is no measurement
data_depth$depth_corrected[is.na(data_depth$`depth_measured (cm)`)]=NA
data_depth$transect_pos_cm = data_depth$knot*10-10


library(tidyr)
library(dplyr)


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

data_measurements = left_join(data_measurements,data_depth, by =c("transectlabel","transect_pos_cm"))
data_measurements$depth_interpolated = round(data_measurements$depth_interpolated, 2)


######


data_measurements$depth_interpolated = data_measurements$depth_interpolated*-1



data_fucus = data_measurements %>% filter(Morphology != "A. nodosum")

#####Do traits of fucus change on a depth gradient

#TDMC
m1 = lm(TDMC ~ depth_interpolated, data = data_fucus)
summary(m1)
m1_length = lm(TDMC ~ depth_interpolated + length_mm, data = data_fucus)
summary(m1_length)

#TDMC plot
plot_tdmc <- ggscatter(x = "depth_interpolated", y = "TDMC",
                       xlab = "Depth (cm)", ylab = "TDMC",
                       color = "Morphology",
                       palette = c("#3399FF", "#FF6666"),
                       data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m1$coefficients[1], slope = m1$coefficients[2], size = 1) +
  guides(color = guide_legend(title = "Morphology",
                              label.theme = element_text(size = 15)))

cor.test(data_fucus$TDMC , data_fucus$depth_interpolated, method = "pearson")


#SA:P
m2 = lm(SA_P ~ depth_interpolated, data = data_fucus)
summary(m2)
m2_length = lm(SA_P ~ depth_interpolated + length_mm, data = data_fucus)
summary(m2_length)

#SA:P plot
plot_sap <- ggscatter(x = "depth_interpolated", y = "SA_P",
                      xlab = "Depth (cm)", ylab = "SA:P",
                      color = "Morphology",
                      palette = c("#3399FF", "#FF6666"),
                      data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m2$coefficients[1], slope = m2$coefficients[2], size = 1) +
  guides(color = FALSE)

cor.test(data_fucus$SA_P , data_fucus$depth_interpolated, method = "pearson")


#LW ratio
m3 = lm(LW_ratio ~ depth_interpolated, data = data_fucus)
summary(m3)
m3_length = lm(LW_ratio ~ depth_interpolated + length_mm, data = data_fucus)
summary(m3_length)

#LW ratio plot
plot_lw <- ggscatter(x = "depth_interpolated", y = "LW_ratio",
                     xlab = "Depth (cm)", ylab = "Length width ratio",
                     color = "Morphology",
                     palette = c("#3399FF", "#FF6666"),
                     data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m3$coefficients[1], slope = m3$coefficients[2], size = 1) +
  guides(color = FALSE)

#1. Different species occur in different depths
#2. Serration is a suboptimal identifyer in small life stages are there other things changing with depth?
m3 = lm(LW_ratio ~ depth_interpolated, data = data_fucus)
summary(m3)
#3. Lenght width ratio also changes with depth

#4 is this just because the shallow germlings are smaller? - control for length
m3 = lm(LW_ratio ~ depth_interpolated + length_mm, data = data_fucus)
summary(m3)

#5 no, it still persists if we control for length - p = 0.00448

t.test(LW_ratio ~ Morphology, data = data_fucus)

t.test(LW_ratio ~ Morphology, data = data_fucus[data_fucus$length_mm > 30,])
ggboxplot(x = "Morphology", y = "LW_ratio",
          xlab = "Morphology", ylab = "Length width ratio",
          color = "Morphology", palette = c("#FF6666", "#3366CC" ),
          data = data_fucus[data_fucus$length_mm > 40,])

cor.test(data_fucus$LW_ratio, data_fucus$depth_interpolated, method = "pearson")

#LP ratio
m4 = lm(LP_ratio ~ depth_interpolated, data = data_fucus)
summary(m4)
m4_length = lm(LP_ratio ~ depth_interpolated + length_mm, data = data_fucus)
summary(m4_length)

#LP ratio plot
plot_lp <- ggscatter(x = "depth_interpolated", y = "LP_ratio",
                     xlab = "Depth (cm)", ylab = "Length perimeter ratio",
                     color = "Morphology",
                     palette = c("#3399FF", "#FF6666"),
                     data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m4$coefficients[1], slope = m4$coefficients[2], size = 1) +
  guides(color = FALSE)

cor.test(data_fucus$LP_ratio , data_fucus$depth_interpolated, method = "pearson")

#Something happens to the LP ratio when the germlings reaches 25 mm, maybe Serratus developes serration then?
ggscatter(x = "depth_interpolated", y = "LP_ratio",
          title = "",
          xlab = "Depth (cm)", ylab = "Length perimeter ratio",
          color = "Morphology", palette = c("darkgreen", "orange"),
          data = data_fucus[data_fucus$length_mm>25,], add = "reg.line")

ggscatter(x = "depth_interpolated", y = "LP_ratio",
          title = "",
          xlab = "Depth (cm)", ylab = "Length perimeter ratio",
          color = "Morphology", palette = c("darkgreen", "orange"),
          data = data_fucus, add = "reg.line")
# when they are bigger the lp ratio is differnt, when they are small, it is not


#STA(mm2/g-1)
m5 = lm(STAmm2_g ~ depth_interpolated, data = data_fucus)
summary(m5)
m5_length = lm(STAmm2_g ~ depth_interpolated + length_mm, data = data_fucus)
summary(m5_length)

#STA plot
plot_sta <- ggscatter(x = "depth_interpolated", y = "STAmm2_g",
                     xlab = "Depth (cm)", ylab = "STA (mm^2/g)",
                     color = "Morphology",
                     palette = c("#3399FF", "#FF6666"),
                     data = data_fucus, size = 4, alpha = 0.7) +
  geom_abline(intercept = m5$coefficients[1], slope = m5$coefficients[2], size = 1) +
  guides(color = FALSE)

cor.test(data_fucus$STAmm2_g , data_fucus$depth_interpolated, method = "pearson")

#Combined plots - arranged
combined_plots <- plot_tdmc + plot_sap + plot_lw + plot_lp + plot_sta +
  plot_layout(ncol = 3, nrow = 2) +
  plot_annotation(title = "Traits and depth correlations", theme = theme(plot.title = element_text(size = 15))) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15))

combined_plots


####Is there a difference between exposed and sheltered
m6 = lm(LW_ratio ~ Condition, data = data_fucus)
summary(m6)
m6 = lm(LW_ratio ~ Condition + length_mm, data = data_fucus)
summary(m6)

#LW_ratio
ggboxplot(x = "Condition", y = "LW_ratio", color = "Condition", 
          xlab = "Condition", ylab = "Length width ratio",
          palette = c("#FF6699", "#0099FF"),
          data = data_fucus)

#TDMC
m7 = lmer(TDMC ~ Condition + length_mm + (1 | transectlabel), data = data_fucus)
summary(m7)
ggscatter(x = "depth_interpolated", y = "TDMC", 
          color = "Condition", add = "reg.line",
          xlab = "Depth (cm)", ylab = "TDMC",
          data = data_fucus, 
          palette = c("#009900", "#0099FF"),
          ellipse = TRUE, mean.point = TRUE,
          star.plot = TRUE)

ggboxplot(x = "transectlabel", y = "TDMC", color = "Condition", data = data_fucus)

#Is it possible to identify germlings based on their traits
#Is serratus distinguishable by traits - most likely with area

ggboxplot(x = "Condition", y = "depth_interpolated", data = data_measurements)

#pca
library(car)
library(ggfortify)
library(MuMIn)

m8 = lmer(TDMC ~ scale(length_mm) + scale(depth_interpolated) + (1 | transectlabel), data = data_fucus)
summary(m8)
r.squaredGLMM(m8)
qqp(resid(m8))

data_pca = data_fucus %>% select(TDMC, SA_P, LW_ratio, LP_ratio, STAmm2_g, wetweight_g, depth_interpolated, Condition)
data_pca = na.omit(data_pca)
data_pca_result = prcomp(data_pca[c(-6, -7, -8)], scale = TRUE)

#PCA screeplot
eigenvalues <- data_pca_result$sdev^2
pcs <- 1:length(eigenvalues)
regression_line <- lm(eigenvalues ~ pcs)

screeplot(data_pca_result, main = "PCA", ylab = "Eigenvalues", xlab = "PCs")

variance_percentage <- eigenvalues / sum(eigenvalues) * 100
max_height <- max(data_pca_result$sdev)
label_position <- eigenvalues + 0.000005 * max_height

labels <- c("38%", "29%", "17%", "12%", "4%")
label_pcs <- c(1, 2, 3, 4, 5)

text(label_pcs, label_position - 0.05 * max_height, labels, pos = 3, adj = c(0.5, 0.5))

abline(regression_line, col = "black")


data_pca_x = data.frame(data_pca_result$x)

ggscatter(x = "PC1", y = "PC2", color = "PC3", add = "reg.line", data = data_pca_x)

autoplot(data_pca_result, loadings = TRUE, loadings.label = TRUE, color = "Condition", size = "depth_interpolated", data = data_pca)


#Does the distribution of germlings look different from the adult distribution, filter for sheltered

####Multi density chart
tra_dat = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "adult", col_types = "ccccccccnnncnncc")

# subset the data needed to interpolate the depth to the different points
# correct the depth by the water level from the nearby station
depth_data <- 
  tra_dat %>% 
  mutate(depth_correct = (water_level_cm + depth_cm) ) %>%
  select(date, transect_id, position, depth_correct) %>%
  distinct()

# problems with the start of transect 2, remove positions 0 to 4
# we do not have a starting depth
depth_data <- 
  depth_data %>%
  filter( !(transect_id == 2 & position %in% 0:4) )

# split into a list
depth_list <- split(depth_data, depth_data$transect_id)

# loop over all transects
depth_out <- vector("list", length = length(depth_list))
for(i in 1:length(depth_list)) {
  
  # initialise a data.frame to work with
  df <- depth_list[[i]]
  
  # get the dividers
  dividers <- which(!is.na(df$depth_correct) )
  
  # duplicate middles for which the data are needed for multiple calculations
  if (length(dividers) > 2) {
    
    dups <- dividers[-c(1, length(dividers))]
    
    df <- 
      df[c(1:nrow(df), dups), ] %>%
      arrange(date, transect_id, position)
    
  }
  
  # add an ID column
  x <- vector("list", length = (length(dividers)-1))
  for(j in 1:(length(dividers)-1) ) {
    
    x[[j]] <- rep(j, ( (dividers[j+1] - dividers[j])+1 ) )
    
  }
  
  # write the ID column into the data.frame
  df$position_id <- unlist(x)
  
  # split data by the position ID
  df_list <- split(df, df$position_id)
  
  # for each block, interpolate the depth between the points
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

# bind the list into a data.frame
depth_out <- bind_rows(depth_out)

# join this back to the full dataset
tra_a <- 
  full_join(tra_dat, depth_out, by = c("date", "transect_id", "position")) %>%
  select(date, transect_id, site_code, time, position, water_level_cm, 
         depth_cm, depth_correct, depth_interpolated, binomial_code, length_cm, circum_cm, 
         field_observer, notes)

# check the summary statistics
summary(tra_a)

# check for unique values, especially for the binomial codes
lapply(tra_a, function(x) unique(x))

# remove missing binomial codes and the missing values "", NA
tra_a <- 
  tra_a %>%
  filter( !(binomial_code %in% c("-9999", "") | is.na(binomial_code) | is.na(depth_interpolated) )  )

# work with the transect data
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

# load the raw allometric data
allo_dat = read_sheet("https://docs.google.com/spreadsheets/d/1fNnrs96uAcEXcz0evvIH6nLL-HG18aKn-hUdE6uLRuw/edit#gid=0", sheet = "adult2", col_types = "ccccdcccccdcddddccdddccdccc")


# remove the missing values from allo_dat
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

# plot out the depth distribution
allo_dat %>%
  ggplot(data = ., 
         mapping = aes(x = binomial_code, y = depth_correct)) +
  geom_point()

# get summary statistics
allo_summary <- 
  allo_dat %>%
  group_by(binomial_code) %>%
  summarise(mean_depth = mean(depth_correct, na.rm = TRUE),
            min_depth = min(depth_correct, na.rm = TRUE),
            max_depth = max(depth_correct, na.rm = TRUE),
            quant_20 = quantile(depth_correct, 0.20, na.rm = TRUE),
            quant_80 = quantile(depth_correct, 0.80, na.rm = TRUE))

# add an allometry identifier
allo_summary$data_id <- "allometry_data"

allo_summary <- 
  allo_summary %>%
  select(data_id, binomial_code, mean_depth:quant_80)

# we used these summary statistics to choose our experimental depths

# plot the depth data for each species
names(allo_dat)
names(tra_a)

# bind these data together
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

#plot of all 4 species
#plot of morphology of germlings and adults

# Filter the relevant data
adult_data <- all_depth %>% 
  filter(!is.na(binomial_code) & !is.na(depth_correct) & data_set == "transect")

# Create a named vector for custom legend labels
legend_labels <- c("fu_ve" = "F. vesiculosus",
                   "fu_se" = "F. serratus",
                   "fu_sp" = "F. spiralis",
                   "as_no" = "A. nodosum")

# Create the multi-density plot with positive x-axis values
chart <- ggplot(adult_data, aes(x = abs(depth_correct), fill = binomial_code)) +
  geom_density(alpha = 0.5) +
  labs(title = "Distribution of Adult Species",
       x = "Depth (cm)",
       y = "Density (%)",
       fill = "Species") +
  theme(legend.position = "bottom") +
  scale_fill_manual(values = c("fu_sp" = "#99FF33",
                               "fu_ve" = "#6699FF",
                               "as_no" = "#FF99CC",
                               "fu_se" = "#FF6666"),
                    labels = legend_labels) +
  scale_x_continuous(labels = abs) + 
  scale_y_percent()

print(chart)

#Adult morphology distribution - flip and reverse and stuff and add aschophyllum
ggplot(adult_data, aes(x = depth_correct, fill = serrated)) +
  geom_density(alpha = 0.7) +
  xlab("Depth (cm)") +
  ylab("Density (%)") +
  scale_fill_manual(values = c("#FF99CC", "#6699FF", "#FF6666"), drop = TRUE) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) + 
  guides(fill = guide_legend(title = "Morphology",
                             label.theme = element_text(size = 15))) + 
  ggtitle("Distribution of Adult Morphology")


#Exposed and Sheltered Germlings 

#Sheltered germlings
# Remove NA values from Morphology variable
data_sheltered = data_measurements %>% filter(Condition != "Exposed")
data_sheltered <- data_sheltered[!is.na(data_sheltered$Morphology), ]

ggplot(data_sheltered, aes(x = depth_interpolated, fill = Morphology)) +
  geom_density(alpha = 0.7) +
  xlab("Depth (cm)") +
  ylab("Density (%)") +
  scale_fill_manual(values = c("#FF99CC", "#6699FF", "#FF6666"), drop = TRUE) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) + 
  guides(fill = guide_legend(title = "Morphology",
                             label.theme = element_text(size = 15))) + 
  ggtitle("Distribution of Sheltered Germling Morphology")

#Exposed germlings - NA still there?

data_exposed = data_measurements %>% filter(Condition != "Sheltered")
data_exposed <- data_exposed[!is.na(data_exposed$Condition), ]

ggplot(data_exposed, aes(x = depth_interpolated, fill = Morphology)) +
  geom_density(alpha = 0.7) +
  xlab("Depth (cm)") +
  ylab("Density (%)") +
  scale_fill_manual(values = c("#6699FF", "#FF6666"), drop = TRUE) +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 15)) + 
  guides(fill = guide_legend(title = "Morphology",
                             label.theme = element_text(size = 15))) + 
  ggtitle("Distribution of Exposed Germling Morphology")
