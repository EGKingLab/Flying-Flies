args = commandArgs(trailingOnly = TRUE)
target_chrom <- args[1]

message(paste("Start Chromosome:", target_chrom))

library(tidyverse)

FF_genomes_filtered <- readRDS("../../../Big_data/FF_genomes_filtered.RDS") %>% 
  filter(CHROM == target_chrom)

FF_matrix_raw <- FF_genomes_filtered %>% 
  select(starts_with("freq_")) %>% 
  as.matrix()
  
FF_depth <- FF_genomes_filtered %>% 
  select(starts_with("N_")) %>% 
  as.matrix()

# linear model
lm_out <- numeric(length = nrow(FF_matrix))
# model comparison of linear models
comp_out <- numeric(length = nrow(FF_matrix))

# t.test unequal variance (Welch) - default t.test
tue_out <- numeric(length = nrow(FF_matrix))
# t.test equal variance
te_out <- numeric(length = nrow(FF_matrix))

# quasibinomial without weights
qb_b_out <- numeric(length = nrow(FF_matrix))
# quasibinomial with weights
qb_w_out <- numeric(length = nrow(FF_matrix))


perc_tracker <- round(nrow(FF_matrix) * 0.05)

paste0("Progress: START")

for(ii in 1:nrow(FF_matrix)){
  FF_temp <- data.frame(
    "Treatment" = rep(c("CON","FP"), each = 4), 
    "Frequency" = FF_matrix[ii,],
    "Depth" = FF_depth[ii,])
  
  # Arcsine squareroot transform data before fitting linear models
  FF_temp <- FF_temp %>%
    mutate(Freq_asqrt = asin(sqrt(Frequency)))
  
  # linear model with anova
  lm_FF_temp <- lm(data = FF_temp, 
                    Freq_asqrt ~ Treatment)
  
  anova_lm_FF_temp <- anova(lm_FF_temp)
  
  lm_out[ii] <- anova_lm_FF_temp$`Pr(>F)`[1]
  
  # lm model comparison
  comp_lm_mean <- lm(data = FF_temp,
                     Freq_asqrt ~ 1)
  comp_anova <- anova(comp_lm_mean, lm_FF_temp)
  
  comp_out[ii] <- comp_anova$`Pr(>F)`[2]
  
  
  # t test Welch's (default, unequal variances)
  tue_FFtemp <- t.test(data = FF_temp, 
                         Freq_asqrt ~ Treatment)
  tue_out[ii] <- tue_FFtemp$p.value
  
  # t test equal var
  te_FFtemp <- t.test(data = FF_temp,
                       Freq_asqrt ~ Treatment, var.equal = TRUE)
  te_out[ii] <- te_FFtemp$p.value
  
  # Following GLMs do not use arcsine-squareroot transform
  # quasibinomial model - WITHOUT weights
  qb_b_FFtemp <- glm(data = FF_temp,
                     Frequency ~ Treatment,
                     family = quasibinomial(link = "logit")
                     )
  
  qb_b_out[ii] <- as.double(broom::tidy(qb_b_FFtemp)[2,5])
  
  # quasibinomial model - WITH weights
  qb_w_FFtemp <- glm(data = FF_temp,
                           Frequency ~ Treatment,
                           family = quasibinomial(link = "logit"),
                           weights = Depth)
  
  qb_w_out[ii] <- as.double(broom::tidy(qb_w_FFtemp)[2,5])
  
  
  
  if (ii %% perc_tracker == 0){
    current_per <- (ii/nrow(FF_matrix) * 100)
    cat("Progress: ", current_per, "%\n")
  }
}

FF_results <- FF_genomes_filtered[, c(1, 2)]

FF_results$lm_Pval <- lm_out
FF_results$comp_Pval <- comp_out
FF_results$tue_Pval <- tue_out
FF_results$te_Pval <- te_out
FF_results$qb_b_Pval <- qb_b_out
FF_results$qb_w_Pval <- qb_w_out

saveRDS(FF_results, 
        file = paste0("../../../Big_data/SNP_diff", target_chrom, ".RDS"))
