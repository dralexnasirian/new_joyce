library(readr)
library(ggplot2)
library(tibble)
library(dplyr)
library(readxl)
library(tidyr)
library(openxlsx)

# uploading CSV files ####

setwd("C:/Users/alexn/Downloads/")

# List all CSV files in the directory
csv_files <- list.files(pattern = "\\.csv$")

# Read all CSV files into separate data frames
list_of_dfs <- lapply(csv_files, read.csv)

# Merge all data frames into a single data frame
df <- bind_rows(list_of_dfs)
df <- as_tibble(df)
head(df)

dim(df)

# Wrangling ####

df$noMulti <- df$noSkill
# add the skill per time period
df$noCas_PerTime <- df$noCas_perSkill * df$J

head(df)

df <- df %>% 
  select(mean_zp, mean_zc, mean_zm, std_D, H, J, T, noPer, noMulti, noCas_perSkill, noCas_PerTime ) %>% 
  mutate(totSkill = df$noPer + df$noMulti ) %>%
  mutate(totSkill_2 = df$noPer + df$noMulti + df$noCas_PerTime ) %>%
  mutate(ratio = (H/J) ) %>% 
  mutate(availMulti = (H*J) - H  ) %>% 
  mutate(perPercent = noPer / H ) %>% 
  mutate(multiPercent = noMulti / availMulti  ) %>% 
  mutate(casPercent = noCas_PerTime / J)  %>%
  mutate(ratio_HJ = case_when(
    ratio <= 1.17 ~ ' H/J <= 1.4',
    TRUE ~ '1.4 < H/J'
  )) %>% 
  mutate(r1 = mean_zc / mean_zp ) %>% 
  mutate(r2 = mean_zp / mean_zm) %>% 
  filter(J >= 5)   #exclude J < 5
  
  

df %>% arrange(desc(T ))

colnames(df)

# Find the values for 33%, 66%, and 100% quantiles
q_33 <- quantile(df$ratio, 0.33)
q_66 <- quantile(df$ratio, 0.66)
q_100 <- quantile(df$ratio, 1)
cat("Value for which 33% of the data is less than:", q_33, "\n")
cat("Value for which 66% of the data is less than:", q_66, "\n")
cat("Value for which 100% of the data is less than:", q_100, "\n")

# --- -------------------- ---
# write to excel ####
# --- -------------------- ---
setwd("C:/Users/alexn/Downloads")
write.xlsx(df, "df.xlsx")

df

# --- -------------------- ---
# Stats ####
# --- -------------------- ---

hist(df$ratio)

hist(df$H)

table(df$H)
table(df$J)

ggplot(df, aes(x=J)) + geom_histogram()

# --- ------------------------------------------------ ---
# T ####
# --- ------------------------------------------------ ---

# df_filtered <- df

df_filtered <- df %>% 
  filter(H <= 10)

df_filtered$H <- ifelse(df_filtered$H == 10, paste("H =", df_filtered$H), paste(" H =", df_filtered$H))

# df_filtered$H <- paste("H =", df_filtered$H)

# Plot with customized facet labels
ggplot(df_filtered, aes(x = T, y = noMulti, color = ratio_HJ, shape = ratio_HJ)) +
  geom_jitter(alpha = 0.5) +
  facet_wrap(~H, ncol = 3) + # Use custom labeller
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) + 
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Time") +
  ylab("Number of secondary skills") +
  labs(color = NULL)+
  coord_cartesian(ylim = c(0, 30))

ggplot(df, aes(x = T, y = noPer, color = ratio_HJ, shape = ratio_HJ)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~H, ncol = 5) + # Use custom labeller
  geom_smooth(se = FALSE) + 
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Time") +
  ylab("Number of permanent workers") +
  labs(color = NULL)

ggplot(df, aes(x = T, y = noCas_perSkill , color = ratio_HJ, shape = ratio_HJ)) +
  geom_point(alpha = 0.5) +
  facet_wrap(~H, ncol = 5) + # Use custom labeller
  geom_smooth(se = FALSE) + 
  theme_bw() +
  theme(legend.position = "none") +
  xlab("Time") +
  ylab("Average number of casual workers per skill") +
  labs(color = NULL)

# --- -------------------- ---
# std_D ####
# --- -------------------- ---


df_filtered <- df %>% 
  filter(T > 120) #%>% 
  filter(ratio > 1.2)
  
df_filtered$J <- ifelse(df_filtered$J == 10, paste("J =", df_filtered$J), paste(" J =", df_filtered$J))  

hist(df_filtered$std_D)

df_filtered <- df_filtered %>%
  mutate(group = ifelse(std_D < 0.6, "std_D < 0.6", "std_D > 0.6"))
result <- df_filtered %>%
  group_by(group) %>%
  summarize(
    max_noMulti = max(noMulti),
    min_noMulti = min(noMulti),
    avg_noMulti = mean(noMulti)
  )

max(df_filtered$noCas_PerTime)
min(df_filtered$noCas_PerTime)
hist(df_filtered$noCas_PerTime)

ggplot(df_filtered, aes(x = std_D, y = noMulti, color = noCas_PerTime, alpha = 0.5)) + 
  geom_jitter(width = 0.2, height = .3) +
  facet_wrap(~J) +
  geom_smooth(method = "lm", formula = y ~ poly(x, 2), se = FALSE) +
  scale_color_gradient2(low = "lightblue", mid = "lightblue3", high = "darkblue", 
                        midpoint = 2.6, limits = c(0, 5.4)) +
  theme_bw() +
  xlab("Standard deviation of demand") +
  ylab("Number of secondary skills")  +
  guides(alpha = FALSE) +
  theme(legend.position = "top") +  
  labs(color = "Color code for the average number of casual workers per time") 


ggplot(df_filtered,aes(x=std_D, y=noPer)) + 
  geom_jitter()+
  facet_wrap(~J)

ggplot(df_filtered,aes(x=std_D, y=noCas_PerTime)) + geom_jitter()

# --- ------------------------------------------ ---
ggplot(df_filtered, aes(x=T, y=noPer, color = ratio_HJ)) +
  geom_point() +
  facet_wrap(~J)

ggplot(df_filtered, aes(x=T, y=totSkill, color = ratio_HJ)) +
  geom_point() +
  facet_wrap(~J)

# --- ------------------------------------------ ---
# no Permanent ####
# --- ------------------------------------------ ---

df_filtered <- df %>% 
  filter(J >= 5) %>% 
  filter(T>= 120)

df_filtered$J <- ifelse(df_filtered$J == 10, paste("J =", df_filtered$J), paste(" J =", df_filtered$J))

ggplot(df_filtered, aes(x = noPer, y = noMulti, color = noCas_perSkill)) + 
  geom_jitter(width = 0.4, height = 1) +
  facet_wrap(~J) +
  scale_color_gradient2(low = "lightblue", mid = "lightblue3", high = "darkblue", 
                        midpoint = 0.3, limits = c(0, 0.72)) +
  theme_bw() + 
  xlab("Number of available permanent workers") + 
  ylab("Number of secondary skills") +
  scale_x_continuous(labels = function(x) ifelse(x %% 2.5 == 0, x*2, x*2.5)) +
  labs(color = "Color code for the average number of casual workers per time per skill") +
  theme(legend.position = "top")


ggplot(df_filtered, aes(x=noPer, y= noMulti, color = noCas_perSkill)) + 
  geom_jitter()+
  facet_wrap(~J) +
  scale_color_gradient2(low = "lightblue", mid = "lightblue3", high = "darkblue", 
                        midpoint = 0.3, limits = c(0, 0.72)) +
  theme_bw() + 
  xlab("Number of available permanent workers") + 
  ylab("Number of secondary skills")




# --- ------------------------------------------ ---
# z^c, z^p & z^m ####
# --- ------------------------------------------ ---

df_filtered <- df %>% 
  filter(T > 120) %>% 
  filter(ratio > 1.3)

df_filtered$J <- ifelse(df_filtered$J == 10, paste("J =", df_filtered$J), paste(" J =", df_filtered$J))

ggplot(df_filtered, aes(x=r1, y=noMulti, color= noCas_perSkill, alpha = 0.5 )) + 
  geom_jitter(width = 0.3, height = 0.6) + 
  facet_wrap(~J) + 
  theme_bw() + 
  xlab("mean_zc / mean_zp")

ggplot(df_filtered, aes(x=r2, y=noMulti, color= noCas_perSkill, alpha = 0.5 )) + 
  geom_jitter(width = 0.3, height = 0.6) + 
  facet_wrap(~J)+ 
  theme_bw() + 
  xlab("mean_zp / mean_zm")





ggplot(df_filtered, aes(x=totSkill_2 )) + geom_histogram() 


# partitioning by H
ggplot(df, aes(x = std_D , y = noPer , color = ratio  )) + 
  geom_jitter(width = 0.1, height = 0.1) + 
  facet_wrap(~H) +
  scale_color_gradient(low = "lightblue", high = "darkblue")  # Reverse the color scale

ggplot(df, aes(x = std_D , y = noMulti , color = J  )) + 
  geom_jitter(width = 0.1, height = 0.1) + 
  facet_wrap(~H) +
  scale_color_gradient(low = "lightblue", high = "darkblue")

ggplot(df, aes(x = std_D , y = casPercent, color = J  )) + 
  geom_jitter(width = 0.1, height = 0.1) + 
  facet_wrap(~H) +
  scale_color_gradient(low = "lightblue", high = "darkblue")


hist(df$T)

# --- --------------------------------------------------------- ---
# draft
# --- --------------------------------------------------------- ---




