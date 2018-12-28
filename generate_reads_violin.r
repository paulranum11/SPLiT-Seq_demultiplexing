library(readr)
numlines_data <- read_table2("linespercell.txt", 
                            col_names = FALSE)
numlines <- numlines_data$X1
numreads <- numlines/4
avg <- mean(as.numeric(numreads))
stdev <- sd(numreads)
numreads_data <- data.frame(numreads, avg, stdev)
colnames(numreads_data) <- c("reads", "avg", "stdev")

library(ggplot2)
# Basic violin plot
p <- ggplot(numreads_data, aes(x="Cells", y=reads)) + 
  geom_violin()
p + geom_jitter(shape=16, position=position_jitter(0.2))

ggsave("reads_violin.pdf", width=5, height=5)

# Basic violin plot
numreads_data["log_reads"] <- log(numreads_data$reads)
p <- ggplot(numreads_data, aes(x="Cells", y=log_reads)) + 
  geom_violin()
p + geom_jitter(shape=16, position=position_jitter(0.2))

ggsave("log_reads_violin.pdf", width=5, height=5)

