# Code Appendix

# Step 1: Look at variables and models holistically, see what stands out
## Predictors: Look for correlation, outliers
## Classification: See best WCA performers, least variable
training <- read.csv("training.csv", stringsAsFactors = TRUE)
sort(abs(cor(training)["class", ]), decreasing = TRUE)[2:19]
training$class <- factor(training$class)
levels(training$class) <- c("NG", "OG", "TSG")
dim(training)
names(training)[c(1, 99)]
barplot(table(training$class))
table(training$class) / nrow(training)
any(is.na(training))
library(ggplot2)
scatter <- function(var) {
  ggplot(training, aes_string(var, 1)) +
    geom_jitter(width = 0.05, height = 0.1, size = 0.1,
                colour = rgb(0, 0, 0, alpha = 1 / 3))
}
scat_plot <- lapply(names(training)[-99], scatter)
library(gridExtra)
grid.arrange(grobs = scat_plot[1:20], ncol = 4)
grid.arrange(grobs = scat_plot[21:40], ncol = 4)
grid.arrange(grobs = scat_plot[41:60], ncol = 4)
grid.arrange(grobs = scat_plot[61:80], ncol = 4)
grid.arrange(grobs = scat_plot[81:98], ncol = 4)
library(ggplot2)
scatter <- function(var) {
  ggplot(training, aes_string(var, "class")) +
    geom_jitter(width = 0.05, height = 0.1, size = 0.1,
                colour = rgb(0, 0, 0, alpha = 1 / 3))
}
scat_plot <- lapply(names(training)[-99], scatter)
library(gridExtra)
grid.arrange(grobs = scat_plot[1:20], ncol = 4)
grid.arrange(grobs = scat_plot[21:40], ncol = 4)
grid.arrange(grobs = scat_plot[41:60], ncol = 4)
grid.arrange(grobs = scat_plot[61:80], ncol = 4)
grid.arrange(grobs = scat_plot[81:98], ncol = 4)
sig <- logical(98)
names(sig) <- names(training)[-99]
k <- 1
diffs <- logical(98)
for (var in names(training)[-99]) {
  model <- aov(training[[var]] ~ factor(training$class))
  sig[k] <- summary(model)[[1]][1, 5]
  diffs[k] <- all(TukeyHSD(model)$`factor(training$class)`[, 4] < 0.05)
  k <- k + 1
}
sort(sig[diffs])
score <- function (conf_mat) {
  print(sum(diag(conf_mat) * c(1, 20, 20)))
  print(sum(diag(conf_mat) * c(1, 20, 20)) / sum(apply(conf_mat, 2, sum) * c(1, 20, 20)))
}

classify <- function(probs) {
  if (any(probs[2:3] > 0.05)) {
    subset <- probs[2:3]
    output <- which(subset == max(subset))
    if (length(output) > 1) {
      output <- sample(1:2, 1)
    }
  } else {
    output <- 0
  }
  output
}
library(dplyr)
vars <- training %>% select(Broad_H3K9ac_percentage, N_LOF, pLOF_Zscore,
                            N_Splice, LOF_TO_Total_Ratio, VEST_score, BioGRID_log_degree, Broad_H3K79me2_percentage,
                            RVIS_percentile,
                            Polyphen2, Broad_H3K36me3_percentage, class)
vars$class <- factor(vars$class)
levels(vars$class) <- c("NG", "OG", "TSG")

# Visual display of correlation matrix for predictors
cor_mtx <- round(cor(vars[, names(vars) != "class"]), 2)
library(reshape2)
melted_cor_mtx <- melt(cor_mtx)
cor_heatmap <- ggplot(data = melted_cor_mtx, aes(x = Var1, y = Var2, fill = value)) + geom_tile()
cor_heatmap <- cor_heatmap +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                       midpoint = 0, limit = c(-1, 1), space = "Lab", name="Correlation") + 
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))
cor_heatmap

# Run various fitted models to see which perform best with respect to variability and WCA
scores_mat <- matrix(nrow = 10, ncol = 3)
k <- 1
set.seed(11092020)
for (i in runif(10, min = 1, max = 10^7)) {
  set.seed(i)
  vars <- training %>% select(Broad_H3K9ac_percentage, N_LOF, pLOF_Zscore,
                              N_Splice, LOF_TO_Total_Ratio, VEST_score, BioGRID_log_degree, Broad_H3K79me2_percentage,
                              RVIS_percentile,
                              Polyphen2, Broad_H3K36me3_percentage, class)
  library(caret)
  vars_test <- createDataPartition(vars$class, p = 0.7, 
                                   list = FALSE)
  vars_train <- vars[vars_test, ]
  vars_test <- vars[-vars_test, ]
  
  train_cont <- trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = TRUE)
  
  
  knn_ft <- train(class ~ ., data = vars_train, method = "knn", preProc = c("center", "scale"),
                  trControl = train_cont, tuneGrid = expand.grid(k = 5))
  preds <- predict(knn_ft, newdata = vars_test, type = "prob")
  knn_mod <- table("pred"=unlist(apply(preds, 1, classify)), "obs" = vars_test$class)
  scores_mat[k, 1] <- score(knn_mod)
  train_cont <- trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = TRUE)
  qda_ft <- train(class ~ ., data = vars_train, method = "qda", preProc = c("center", "scale"),
                  trControl = train_cont)
  preds <- predict(qda_ft, newdata = vars_test, type = "prob")
  
  qda_mod <- table("pred"=apply(preds, 1, classify), "obs" = vars_test$class)
  scores_mat[k, 2] <- score(qda_mod)
  train_cont <- trainControl(method = "cv", number = 5, classProbs = TRUE, savePredictions = TRUE)
  lda_ft <- train(class ~ ., data = vars_train, method = "lda", preProc = c("center", "scale"),
                  trControl = train_cont)
  preds <- predict(lda_ft, newdata = vars_test, type = "prob")
  
  lda_mod <- table("pred"=apply(preds, 1, classify), "obs" = vars_test$class)
  scores_mat[k, 3] <- score(lda_mod)
  k <- k + 1
}
colnames(scores_mat) <- c("KNN", "QDA", "LDA")
data.frame(scores_mat)
apply(scores_mat, 2, mean)
apply(scores_mat, 2, sd)

## Step 2: Consider removing observations and correlated predictors before 