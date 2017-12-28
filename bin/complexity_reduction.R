library('Rtsne')
directory <- "../data/processed_raw_counts.tsv"
data <- read.table(directory, sep='\t', header=TRUE, row.names=1)
h <- sdata[,grepl("H",names(data))]
s <- sdata[,grepl("S",names(data))]
n <- sdata[,grepl("N",names(data))]
c <- sdata[,grepl("C",names(data))]

set.seed(13)

tsne_h <- Rtsne(as.matrix(unique(h)), theta=0.20)
tsne_s <- Rtsne(as.matrix(unique(s)), theta=0.20)
tsne_n <- Rtsne(as.matrix(unique(n)), theta=0.20)
tsne_c <- Rtsne(as.matrix(unique(c)), theta=0.20)


th <- ggplot(data.frame(tsne_h$Y), aes(X1, X2)) + 
  geom_point(alpha=0.55, size=1.25) +
  ggtitle ("Healthy")
ggsave("../pictures/tSNE_healthy.png", plot=last_plot(), scale=1.5,
       dpi=450, width=5, height=3)

ts <- ggplot(data.frame(tsne_s$Y), aes(X1, X2)) + 
  geom_point(alpha=0.55, size=1.25) +
  ggtitle ("Steatosis")
ggsave("../pictures/tSNE_steatosis.png", plot=last_plot(), scale=1.5,
       dpi=450, width=5, height=3)

tn <- ggplot(data.frame(tsne_n$Y), aes(X1, X2)) + 
  geom_point(alpha=0.55, size=1.25) +
  ggtitle ("NASH")
ggsave("../pictures/tSNE_NASH.png", plot=last_plot(), scale=1.5,
       dpi=450, width=5, height=3)

tc <- ggplot(data.frame(tsne_c$Y), aes(X1, X2)) + 
  geom_point(alpha=0.55, size=1.25) +
  ggtitle ("HCC")
ggsave("../pictures/tSNE_HCC.png", plot=last_plot(), scale=1.5,
       dpi=450, width=5, height=3)