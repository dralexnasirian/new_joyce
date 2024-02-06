
library(readr)
library(ggplot2)
library(dplyr)
#df_3 <- read_csv("C:/Users/alexn/Downloads/df_3.csv")
df_4 <- read_csv("C:/Users/alexn/Downloads/df_4.csv")
df_5 <- read_csv("C:/Users/alexn/Downloads/df_5.csv")
df_6 <- read_csv("C:/Users/alexn/Downloads/df_6.csv")
df_7 <- read_csv("C:/Users/alexn/Downloads/df_7.csv")
df_8 <- read_csv("C:/Users/alexn/Downloads/df_8.csv")

df <- rbind(df_4, df_5, df_6, df_7, df_8)


df$H <- sprintf("H=J=%d", df$H)


df <- df %>%
  select( -O, -J, -noPer, -costPer, -costCas, -costTrain, -costTot)

df$r3 <- df$mean_zc / df$mean_zm

model <- lm(noSkill  ~ ., df)
model$coefficients
summary(model)

# jitter plot ###########
df

ggplot(df, aes(x = T, y = noSkill, color = mean_zc)) +
  geom_jitter(width = 0.2, height = 0.2, alpha = 0.3) +
  facet_wrap(~ H, ncol = 5) +
  theme_bw() +
  geom_smooth() +
  theme(legend.position = "none") 




# Assuming df is your dataset
rows_with_condition <- df[df$noPer > df$noSkill, ]

# Display the resulting subset
print(rows_with_condition)


ggplot(df, aes(x=T,color=noMulti)) + 
  geom_boxplot() +
  facet_wrap(~ H) +
  theme_classic()


# Assuming df is your dataset
ggplot(df, aes(x = mean_zc, y = mean_zm, size = noSkill)) +
  geom_jitter(width = 0.2, height = 0, color = "blue", alpha = 0.2) +
  scale_size_continuous(range = c(1, 5), breaks = seq(0, max(df$noSkill), by = 1)) +
  labs(title = "Jitter Plot", x = "mean_zc", y = "mean_zm")+
  facet_wrap(~ H)


result <- df %>%
  group_by(H) %>%
  summarize(max_noSkill = max(noSkill))
