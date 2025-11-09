#Instalando paquetes
#install.packages("ggplot2")  # If you haven't installed ggplot2
#install.packages("dplyr")
#install.packages("car")
#install.packages("tidyverse")
#install.packages(c("multcompView", "agricolae"))
library(multcompView)
library(agricolae)
library(tidyverse)
library(car)
library(dplyr)
library(ggplot2)
data <- read.csv("C:/Users/User/Desktop/master_thesis/Biochemical_traits_pseudomonas.csv")
biochemical <- na.omit(data)


#----------------------------------------

# Calculate the index of replicates for each bacteria group
biochemical <- biochemical %>% 
  mutate(index = (diameter_halo + diameter_colony)/(diameter_colony))

# Remove rows with NA values in column1 and column2
clean_biochem <- biochemical %>%
  filter(!is.na(index))



#protease ###########################################
protease <- clean_biochem %>%
  filter(biochemical_test == "Protease")


# Visualize original data
ggplot(protease, aes(sample = index)) + 
  stat_qq() + 
  stat_qq_line() +
  ggtitle("Q-Q Plot for Index")


# Shapiro-Wilk test for normality
# Prueba de Shapiro-Wilk para cada bacteria y variable
by(protease$index, protease$biochemical_test, shapiro.test)


# Levene's test for homogeneity of variances
leveneTest(index ~ strains, data = clean_biochem)


# Anova test for protease 
anova_protease <- aov(index ~ strains, data = protease)
summary(anova_protease)

# Perform Tukey's HSD post-hoc test 
tukey_protease <- TukeyHSD(anova_protease, "strains")

# Print the results
print(tukey_protease)

#ggplot
# Prepare data for plotting
tukey_ggplot_protease <- as.data.frame(tukey_protease$strains)
tukey_ggplot_protease$comparison <- rownames(tukey_ggplot_protease)

# Determine group letters
tukey_groups_phos <- multcompLetters4(anova_protease, tukey_protease)$strains


#phosphate_solubilization ###############################
phosphate_solubilization <- clean_biochem %>%
  filter(biochemical_test == "Phosphate_solubilization")
  
# Visualize original data
ggplot(phosphate_solubilization, aes(sample = index)) + 
  stat_qq() + 
  stat_qq_line() +
  ggtitle("Q-Q Plot for Index")


# Shapiro-Wilk test for normality
# Prueba de Shapiro-Wilk para cada bacteria y variable
by(phosphate_solubilization$index, 
   phosphate_solubilization$biochemical_test, 
   shapiro.test)


# Levene's test for homogeneity of variances
leveneTest(index ~ strains, data = clean_biochem)

# Anova test for phosphate
anova_phosphate <- aov(index ~ strains, data = phosphate_solubilization)
summary(anova_phosphate)


# Perform Tukey's HSD post-hoc test
tukey_phosphate <- TukeyHSD(anova_phosphate, "strains")

# Print the results
print(tukey_phosphate)
  
#ggplot
# Prepare data for plotting
tukey_ggplot_phosphate <- as.data.frame(tukey_phosphate$strains)
tukey_ggplot_phosphate$comparison <- rownames(tukey_ggplot_phosphate)

# Determine group letters
tukey_groups_phos <- multcompLetters4(anova_phosphate, tukey_phosphate)$strains


# Prepare data for plotting
mean_values_pro <- protease %>%
  group_by(strains) %>%
  summarize(
    mean_strains = mean(index),
    sem_strains = sd(index) / sqrt(n())
  )

# Prepare data for plotting
mean_values_pho <- phosphate_solubilization %>%
  group_by(strains) %>%
  summarize(
    mean_strains = mean(index),
    sem_strains = sd(index) / sqrt(n())
  )



