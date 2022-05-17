library(reshape2)
library(ggrepel)
library(ggplot2)
library(ggthemes)

dataf = read.csv("filename.csv")
dataf=dataf[,-1] # Remove the column with row indices

