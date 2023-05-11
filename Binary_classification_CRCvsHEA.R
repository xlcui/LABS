library(pROC)
library(randomForest)
library(caret)

##############################
###### Prepare data ##########
##############################

cohort <- read.csv("Cohort.csv")
cohort_CRC <- cohort[cohort$Type %in% c("CRC", "HEA"), ]
cohort_CRC$group <- cohort_CRC$Type
IDs <- cohort_CRC$Sample
set.seed(123)
train_idx <- createDataPartition(cohort_CRC$group, p = 0.75, list = FALSE)


###### methylation TSS -> PCs ######
methyl_TSS <- read.csv("CRC.HEA.gencode.v19.genes.tss.V2.csv")
gene_names <- methyl_TSS[, 2]
TSS_region <- methyl_TSS[, 1]
methyl_TSS <- methyl_TSS[, -c(1:2)]
methyl_sub <- t(methyl_TSS[, colnames(methyl_TSS)%in%IDs])
colnames(methyl_sub) <- as.factor(paste0("F", 1:ncol(methyl_sub)))

# Split the dataset into training and testing sets

train_data <- methyl_sub[train_idx, ]
test_data <- methyl_sub[-train_idx, ]

# Preprocess data: standardize/normalize features
preprocess_params <- preProcess(train_data, method = c("center", "scale"))
train_data_standardized <- predict(preprocess_params, train_data)
test_data_standardized <- predict(preprocess_params, test_data)

# Perform PCA
pca_train <- prcomp(train_data_standardized, center = TRUE, scale. = TRUE)

# Determine the number of principal components to retain
# (e.g., retain components that explain 95% of the variance)
explained_variance <- cumsum(pca_train$sdev^2) / sum(pca_train$sdev^2)
num_components <- which(explained_variance >= 0.95)[1]

# Get the reduced training data
train_data_pca <- pca_train$x[, 1:num_components]
colnames(train_data_pca) <- paste0("Methyl.", colnames(train_data_pca))

# Apply PCA transformation to the test set
test_data_pca <- predict(pca_train, test_data_standardized)[, 1:num_components]
colnames(test_data_pca) <- paste0("Methyl.", colnames(test_data_pca))


###### copy number -> PCs ######

copy_number <- read.csv("CRC.HEA.CN.csv")
bin_names <- copy_number[, 1]
copy_number <- copy_number[, -1]

flag <- apply(copy_number, 1, function(x){any(is.na(x) | is.infinite(x))})
copy_number <- copy_number[!flag, ]
copy_number_sub <- t(copy_number[, colnames(copy_number)%in%IDs])
colnames(copy_number_sub) <- as.factor(paste0("F", 1:ncol(copy_number_sub)))

# Split the dataset into training and testing sets
set.seed(123)
train_idx <- createDataPartition(cohort_CRC$group, p = 0.75, list = FALSE)
train_data <- copy_number_sub[train_idx, ]
test_data <- copy_number_sub[-train_idx, ]

# Preprocess data: standardize/normalize features
preprocess_params <- preProcess(train_data, method = c("center", "scale"))
train_data_standardized <- predict(preprocess_params, train_data)
test_data_standardized <- predict(preprocess_params, test_data)

# Perform PCA
pca_train <- prcomp(train_data_standardized, center = TRUE, scale. = TRUE)

# Determine the number of principal components to retain
# (e.g., retain components that explain 95% of the variance)
explained_variance <- cumsum(pca_train$sdev^2) / sum(pca_train$sdev^2)
num_components <- which(explained_variance >= 0.95)[1]

# Get the reduced training data
train_data_pca_copy_number <- pca_train$x[, 1:num_components]
colnames(train_data_pca_copy_number) <- paste0("CN.", colnames(train_data_pca_copy_number))
# Transpose the standardized test data

# Apply PCA transformation to the test set
test_data_pca_copy_number <- predict(pca_train, test_data_standardized)[, 1:num_components]
colnames(test_data_pca_copy_number) <- paste0("CN.", colnames(test_data_pca_copy_number))


###### immune prop ######
immune_prop <- read.csv("cell.type.decon.csv")
immune_prop_data <- immune_prop[immune_prop[,1]%in%IDs, -c(1, 7)]

train_prop_data <- immune_prop_data[train_idx, ]
test_prop_data <- immune_prop_data[-train_idx, ]

# Preprocess data: standardize/normalize features
preprocess_params_prop <- preProcess(train_prop_data, method = c("center", "scale"))
train_prop_data_standardized <- predict(preprocess_params_prop, train_prop_data)
test_prop_data_standardized <- predict(preprocess_params_prop, test_prop_data)




##############################
###### Random Forest #########
##############################

set.seed(123)
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)
                        

###### immune prop only
model_prop <- train(Group ~ .,
                      data = data.frame(train_prop_data_standardized, Group = as.factor(cohort_CRC$group[train_idx])),
                      method = "rf",
                      trControl = control,
                      metric = "ROC",
                      tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_prop, test_prop_data_standardized, type = "prob")
roc_test_prop <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_prop <- as.numeric(auc(roc_test_prop))
cat("ROC AUC on test set:", auc_test_prop, "\n")
plot(roc_test_prop, main = "ROC curve on test set")


###### methylation only
model <- train(Group ~ .,
               data = data.frame(train_data_pca, Group = as.factor(cohort_CRC$group[train_idx])),
               method = "rf",
               trControl = control,
               metric = "ROC",
               tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model, test_data_pca, type = "prob")
roc_test <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test <- as.numeric(auc(roc_test))
# ROC AUC on test set: 0.9583333 
cat("ROC AUC on test set:", auc_test, "\n")
plot(roc_test, main = "ROC curve on test set")


###### copy_number only

set.seed(123456)
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

model_copy_number <- train(Group ~ .,
                           data = data.frame(train_data_pca_copy_number, Group = as.factor(cohort_CRC$group[train_idx])),
                           method = "rf",
                           trControl = control,
                           metric = "ROC",
                           tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_copy_number, test_data_pca_copy_number, type = "prob")
roc_test_copy_number <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_copy_number <- as.numeric(auc(roc_test_copy_number))
# ROC AUC on test set: 0.9583333 
cat("ROC AUC on test set:", auc_test_copy_number, "\n")
plot(roc_test_copy_number, main = "ROC curve on test set")


###### methylation + copy_number 

train_data <- data.frame(train_data_pca, train_data_pca_copy_number)
test_data <- data.frame(test_data_pca, test_data_pca_copy_number)

model_tss_cn <- train(Group ~ .,
                      data = data.frame(train_data, Group = as.factor(cohort_CRC$group[train_idx])),
                      method = "rf",
                      trControl = control,
                      metric = "ROC",
                      tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_tss_cn, test_data, type = "prob")
roc_test_tss_cn <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_tss_cn <- as.numeric(auc(roc_test_tss_cn))
cat("ROC AUC on test set:", auc_test_tss_cn, "\n")
plot(roc_test_tss_cn, main = "ROC curve on test set")


###### methylation + copy_number + prop

train_data <- data.frame(train_data_pca, train_data_pca_copy_number, train_prop_data_standardized)
test_data <- data.frame(test_data_pca, test_data_pca_copy_number, test_prop_data_standardized)

model_all <- train(Group ~ .,
                      data = data.frame(train_data, Group = as.factor(cohort_CRC$group[train_idx])),
                      method = "rf",
                      trControl = control,
                      metric = "ROC",
                      tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_all, test_data, type = "prob")
roc_test_all <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_all <- as.numeric(auc(roc_test_all))
cat("ROC AUC on test set:", auc_test_all, "\n")
plot(roc_test_all, main = "ROC curve on test set")


importance <- varImp(model_all)
# Plot variable importance
plot(importance, main = "Variable Importance Plot")

pdf("importance.pdf", width = 6, height = 12)
# Plot variable importance
plot(importance, main = "Variable Importance Plot")
# Close the pdf file
dev.off()


plot(roc_test, col = "blue", main = "ROC Curve: CRC vs Health (Random Forest)")
plot(roc_test_tss_cn, col = "red", add = TRUE)
plot(roc_test_all, col = "black", add = TRUE)
legend("bottomright", legend = c("TSS Methyl Alone", "Methyl + Copy Number", "All"), col = c("blue", "red", "black"), lty = 1)






####################
###### SVM #########
####################


set.seed(123)
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)


###### immune prop only
model_prop <- train(Group ~ .,
                    data = data.frame(train_prop_data_standardized, Group = as.factor(cohort_CRC$group[train_idx])),
                    method = "svmRadial",
                    trControl = control,
                    metric = "ROC",
                    tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_prop, test_prop_data_standardized, type = "prob")
roc_test_prop <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_prop <- as.numeric(auc(roc_test_prop))
cat("ROC AUC on test set:", auc_test_prop, "\n")
plot(roc_test_prop, main = "ROC curve on test set")


###### methylation only
model <- train(Group ~ .,
               data = data.frame(train_data_pca, Group = as.factor(cohort_CRC$group[train_idx])),
               method = "svmRadial",
               trControl = control,
               metric = "ROC",
               tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model, test_data_pca, type = "prob")
roc_test <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test <- as.numeric(auc(roc_test))
# ROC AUC on test set: 0.9583333 
cat("ROC AUC on test set:", auc_test, "\n")
plot(roc_test, main = "ROC curve on test set")


###### copy_number only

set.seed(123456)
control <- trainControl(method = "cv", number = 5, classProbs = TRUE, summaryFunction = twoClassSummary)

model_copy_number <- train(Group ~ .,
                           data = data.frame(train_data_pca_copy_number, Group = as.factor(cohort_CRC$group[train_idx])),
                           method = "svmRadial",
                           trControl = control,
                           metric = "ROC",
                           tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_copy_number, test_data_pca_copy_number, type = "prob")
roc_test_copy_number <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_copy_number <- as.numeric(auc(roc_test_copy_number))
# ROC AUC on test set: 0.9583333 
cat("ROC AUC on test set:", auc_test_copy_number, "\n")
plot(roc_test_copy_number, main = "ROC curve on test set")


###### methylation + copy_number 

train_data <- data.frame(train_data_pca, train_data_pca_copy_number)
test_data <- data.frame(test_data_pca, test_data_pca_copy_number)

model_tss_cn <- train(Group ~ .,
                      data = data.frame(train_data, Group = as.factor(cohort_CRC$group[train_idx])),
                      method = "svmRadial",
                      trControl = control,
                      metric = "ROC",
                      tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_tss_cn, test_data, type = "prob")
roc_test_tss_cn <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_tss_cn <- as.numeric(auc(roc_test_tss_cn))
cat("ROC AUC on test set:", auc_test_tss_cn, "\n")
plot(roc_test_tss_cn, main = "ROC curve on test set")


###### methylation + copy_number + prop

train_data <- data.frame(train_data_pca, train_data_pca_copy_number, train_prop_data_standardized)
test_data <- data.frame(test_data_pca, test_data_pca_copy_number, test_prop_data_standardized)

model_all <- train(Group ~ .,
                   data = data.frame(train_data, Group = as.factor(cohort_CRC$group[train_idx])),
                   method = "svmRadial",
                   trControl = control,
                   metric = "ROC",
                   tuneLength = 5)

# Model performance on the test set
test_preds <- predict(model_all, test_data, type = "prob")
roc_test_all <- roc(cohort_CRC$group[-train_idx], test_preds[, "HEA"])
auc_test_all <- as.numeric(auc(roc_test_all))
cat("ROC AUC on test set:", auc_test_all, "\n")
plot(roc_test_all, main = "ROC curve on test set")



library(kernlab)

# Perform Recursive Feature Elimination (RFE) for SVM
set.seed(123)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
svm_rfe <- rfe(x = train_data, y = as.factor(cohort_CRC$group[train_idx]), sizes = 1:ncol(train_data),
               rfeControl = control, method = "svmRadial", metric = "ROC")
svm_rfe
pdf("SVM_importance.pdf", width = 12, height = 6)
# Visualize feature importance for SVM (RFE method)
plot(svm_rfe$results$Variables, svm_rfe$results$Accuracy, type = "b",
     xlab = "Number of Features", ylab = "Mean Accuracy", main = "Recursive Feature Elimination with SVM")
top_features <- svm_rfe$optVariables[1:10]
top_feature_indices <- match(top_features, colnames(train_data))
text(top_feature_indices, svm_rfe$results$Accuracy[top_feature_indices], top_features, pos = 1, cex = 0.8, col = "red")
dev.off()


plot(roc_test, col = "blue", main = "ROC Curve: CRC vs Health (SVM)")
plot(roc_test_tss_cn, col = "red", add = TRUE)
plot(roc_test_all, col = "black", add = TRUE)
legend("bottomright", legend = c("TSS Methyl Alone", "Methyl + Copy Number", "All"), col = c("blue", "red", "black"), lty = 1)

