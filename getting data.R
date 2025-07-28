# pipelines for the GEO-UC data
library(tidyverse)
library(dplyr)
library(ggplot2)
df <- read.delim("data/GSE243625_RNAseq_colon.txt", check.names = FALSE)
df
row.names(df) <- df$`Gene symbol`
# the above will return error due to 1-mar 2-mar as they are are read by r as two similar values
rownames(df) <- make.unique(df$`Gene symbol`) # this will append .1 , .2
df$`Gene symbol` <- NULL

# create a metadata 
samples <- colnames(df)
condition <- ifelse(grepl("^N", samples), "Control",
                    ifelse(grepl("^I", samples), "Inactive_UC",
                           ifelse(grepl("^A", samples), "Active_UC", NA)))

# Filter out any columns not assigned
keep <- !is.na(condition)
df <- df[, keep]
condition <- condition[keep]

# Create metadata
coldata <- data.frame(row.names = colnames(df),
                      condition = factor(condition, levels = c("Control", "Inactive_UC", "Active_UC")))
# in the above code levels should be ordered the same as the order of ifelse loops in the condition vector
all(colnames(df) == rownames(coldata))  # should return TRUE to check if all the same 
# Subset samples
keep <- coldata$condition %in% c("Active_UC", "Control")
df_sub <- df[, keep]
coldata_sub <- coldata[keep, , drop = FALSE]












