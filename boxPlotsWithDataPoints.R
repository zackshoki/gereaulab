library(dplyr)
library(ggplot2)


all.data <- read.csv("combinedHALOdata.csv")
script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
new_folder <- file.path(script_dir, "Plots")  # name folder
dir.create(new_folder, showWarnings = FALSE)
destination_dir = paste0(script_dir, "/Plots")
# unpack
grik2.pvalb.positive <- all.data$X..Grik2.Pvalb..Positive.Cells
grik2.ntrk3.positive <- all.data$X..Grik2.Ntrk3..Positive.Cells
grik2.positive <- all.data$X..Grik2..Positive.Cells
positive.cells <- all.data$X......Positive.Cells
channel1 <- all.data$Avg.Channel.1.Copies.Per.Cell
channel2 <- all.data$Avg.Channel.2.Copies.Per.Cell
channel3 <- all.data$Avg.Channel.3.Copies.Per.Cell


# plot 1: Grik2+Pvalb+ positive cells
df = data.frame(
  Group = "Grik2+Pvalb+ positive cells",
  Value = grik2.pvalb.positive
)
mean <- mean(df$Value)
sem  <- sd(df$Value) / sqrt(nrow(df))
print(length(grik2.pvalb.positive))
print(summary(grik2.pvalb.positive))
plot1 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "bisque", width = 0.5) + 
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "Grik2+Pvalb+ positive cells", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )


# plot 2: Grik2+Ntrk3+ positive cells
df = data.frame(
  Group = "Grik2+Ntrk3+ positive cells",
  Value = grik2.ntrk3.positive
)

mean <- mean(df$Value)
sem  <- sd(df$Value) / sqrt(nrow(df))

plot2 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "skyblue", width = 0.5) + 
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "Grik2+Ntrk3+ positive cells", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )


# plot 3: Grik2+ positive cells
df = data.frame(
  Group = "Grik2+ positive cells",
  Value = grik2.positive
)

plot3 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "darkkhaki", width = 0.5) +  
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "Grik2+ positive cells", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )



# plot 4: +++Positive cells
df = data.frame(
  Group = "+++Positive cells",
  Value = positive.cells
)


plot4 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "darkorchid1", width = 0.5) + 
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "+++Positive cells", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )



# plot 5: Avg. Channel 1 copies per cell
df = data.frame(
  Group = "Avg. Channel 1 copies per cell",
  Value = channel1
)


plot5 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "chocolate", width = 0.5) + 
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "Avg. Channel 1 copies per cell", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )



# plot 6: Avg. Channel 2 copies per Cell
df = data.frame(
  Group = "Avg. Channel 2 copies per cell",
  Value = channel2
)


plot6 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "darkseagreen3", width = 0.5) + 
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "Avg. Channel 2 copies per cell", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )



# plot 7: Avg. Channel 3 copies per Cell
df = data.frame(
  Group = "Avg. Channel 3 copies per cell",
  Value = channel3
)


plot7 <- ggplot(df, aes(x = Group, y = Value)) +
  geom_bar(stat = "summary", fill = "darkgray", width = 0.5) +  
  geom_jitter(width = 0.1, size = 2, alpha = 0.7) +
  geom_point(aes(x = 1, y = mean), color = "red", size = 3) +
  geom_errorbar(aes(x = 1, ymin = mean - sem, ymax = mean + sem),
                width = 0.1, color = "red") +
  theme_classic() +
  labs(title = "Avg. Channel 3 copies per cell", x = NULL, y = "Number of Cells") +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()   
  )



ggsave(filename = file.path(destination_dir, "grik2_pvalb_positive.jpg"), height = 7, width = 5, plot = plot1, quality = 50)
ggsave(filename = file.path(destination_dir, "grik2_ntrk3_positive.jpg"), height = 7, width = 5, plot = plot2, quality = 50)
ggsave(filename = file.path(destination_dir, "grik2_positive.jpg"), height = 7, width = 5, plot = plot3, quality = 50)
ggsave(filename = file.path(destination_dir, "positive.jpg"), height = 7, width = 5, plot = plot4, quality = 50)
ggsave(filename = file.path(destination_dir, "channel1.jpg"), height = 7, width = 5, plot = plot5, quality = 50)
ggsave(filename = file.path(destination_dir, "channel2.jpg"), height = 7, width = 5, plot = plot6, quality = 50)
ggsave(filename = file.path(destination_dir, "channel3.jpg"), height = 7, width = 5, plot = plot7, quality = 50)
