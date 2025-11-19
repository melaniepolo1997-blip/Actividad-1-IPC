# install.packages(c("multcomp", 'tidyverse', 
#                    'agricolae', 'ggplot2',
#                    'car', 'dplyr'))

#Libraries 
library(multcomp)
library(tidyverse)
library(agricolae)
library(ggplot2)
library(car)
library(dplyr)

#Data loading
setwd('C:/Users/User/Desktop/master_thesis')
fungal_table <- read.csv('fungal_test.csv')
bact_table <- read.csv("bacteria.csv")

#Data curation
fungal_subset <- fungal_table %>% select(strains, pathogen,IR)
bac_subset <- bact_table %>%  select(strains, pathogen,IR)

#Working table including only strains in thesis project
inhibition_rate <- rbind(fungal_subset, bac_subset) %>%
  filter(!strains %in% c("P_protegens", "P_fragi"))

remove(bac_subset, fungal_subset)

# Q-Q plot of data
ggplot(inhibition_rate, aes(sample = IR)) + 
  stat_qq() + 
  stat_qq_line() +
  ggtitle("Q-Q Plot for Inhibition Rate (IR)")

# Levene's test for homogeneity of variances
leveneTest(IR ~ pathogen*strains, data = inhibition_rate)

# Shapiro-Wilk test for normality - Fungi
# Prueba de Shapiro-Wilk for each fungi strain and variable
by(inhibition_rate$IR, inhibition_rate$strains, shapiro.test)
by(inhibition_rate$IR, inhibition_rate$pathogen, shapiro.test)

#ANOVA test
anova_pathogen <- aov(IR ~ strains*pathogen, 
                      data = inhibition_rate)
anova_summary <- summary(anova_pathogen)

#Results from ANOVA in a data frame
#Comparar las medias de los grupos para determinar diferencias 
#estadísticamente significativas entre ellos 
anova_table <- anova_summary[[1]]
anova_p <- as.data.frame(anova_table)
remove(anova_table, anova_summary)


#-----------------------------------------------------------
# AUTOMATIZACIÓN – Tukey + HSD para cada patógeno
#-----------------------------------------------------------

pathogens <- unique(inhibition_rate$pathogen)

tukey_results_list <- list()
hsd_groups_list <- list()

for(p in pathogens){
  
  # Subset de datos
  data_sub <- inhibition_rate %>%
    filter(pathogen == p)
  
  # ANOVA
  model <- aov(IR ~ strains, data = data_sub)
  
  # Tukey HSD
  tukey_res <- TukeyHSD(model)
  tukey_df <- as.data.frame(tukey_res$strains)
  tukey_df$comparison <- rownames(tukey_df)
  
  # HSD agricolae (letras)
  hsd_res <- HSD.test(model, trt = "strains", group = TRUE)
  hsd_df <- as.data.frame(hsd_res$groups)
  hsd_df$strain <- rownames(hsd_df)
  
  # Guardar resultados
  tukey_results_list[[p]] <- tukey_df
  hsd_groups_list[[p]] <- hsd_df
  
}

# Para inspeccionar:
tukey_results_list     # Tukey por patógeno
hsd_groups_list        # Letras de significancia por patógeno

#-----------------------------------------------------------

# Prepare data for plotting
mean_all <- inhibition_rate %>%
  group_by(strains, pathogen) %>%
  summarize(
    mean_strains = mean(IR),
    sem_strains = sd(IR) / sqrt(n()))

mean_all$letters <- c('bc', 'c', 'bc', 'b', 'a', 'bc', 'bc', 'bc', 'a', "bc")

# Plot the results with ggplot2
ggplot(mean_all, aes(x = pathogen, y = mean_strains, fill = strains)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_strains - sem_strains, ymax = mean_strains + sem_strains), 
                position = position_dodge(0.7), width = 0.25)  +
  geom_text(aes(label = letters , y = mean_strains + sem_strains + 1), 
            position = position_dodge(0.7), vjust = 0) +  # Adjust y position for letters
  theme_minimal() +
  labs(title = "Pathogen inhibition rate (IR)",
       x = "Phytopathogens",
       y = "Inhibition %") +
  scale_fill_brewer(palette = "Set1")+
  theme(legend.position = "left")

# Se propone estas modificaciones en el código para mejorar legibilidad y aporte visual
library(ggplot2)

ggplot(mean_all, aes(x = pathogen, y = mean_strains, fill = strains)) +
  geom_bar(stat = "identity", position = position_dodge(), width = 0.7) +
  geom_errorbar(aes(ymin = mean_strains - sem_strains, ymax = mean_strains + sem_strains),
                position = position_dodge(0.7), width = 0.25) +
  geom_text(aes(label = paste0(round(mean_strains, 1), "%"),
                y = mean_strains + sem_strains + 2),
            position = position_dodge(0.7), vjust = 0, size = 3.5, color = "black") +
  theme_minimal(base_size = 14) +
  labs(title = "Pathogen Inhibition Rate (IR)",
       subtitle = "Mean inhibition percentage with SEM error bars",
       x = "Phytopathogens",
       y = "Inhibition (%)",
       fill = "Strains") +
  scale_fill_brewer(palette = "Set1") +
  theme(legend.position = "top",
        plot.title = element_text(face = "bold", hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))

remove(bact_table, fungal_table, 
       inhibition_rate, mean_all, anova_pathogen)





