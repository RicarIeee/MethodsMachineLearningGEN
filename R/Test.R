library(tidyverse)
library(caret)
library(class)
library(gmodels)
library(psych)
library(DynamicCancerDriverKM)
library(rpart)
library(randomForest)
library(e1071)
library(kernlab)


normaldata <-(DynamicCancerDriverKM::BRCA_normal)
dataPt <-(DynamicCancerDriverKM::BRCA_PT)
final_data <- bind_rows(normaldata, dataPt)


rate_min_10 <- final_data %>%
  summarise_all(~ mean(. <400, na.rm = TRUE))


colmn_delate <- names(rate_min_10[, rate_min_10 >= 0.8])


data_filter <- final_data %>%
  select(-one_of(colmn_delate))

data_filter2 <- data_filter

data_PPI<-(DynamicCancerDriverKM::PPI)

data_PPIn <- data_PPI %>%
  pivot_longer(cols = c(`Input-node Gene Symbol`, `Output-node Gene Symbol`), names_to = "variable", values_to = "gen") %>%
  group_by(gen, variable) %>%
  summarise(frecuencia = n()) %>%
  pivot_wider(names_from = variable, values_from = frecuencia, values_fill = 0)

data_PPInR <- data_PPIn %>%
  mutate(total_mode = `Input-node Gene Symbol` + `Output-node Gene Symbol`) %>%
  select(total_mode) %>%
  arrange(desc(total_mode))


print(data_PPInR)

data_filter_x<-colnames(data_filter)[ 8:ncol(data_filter)]
aux2 <- AMCBGeneUtils::changeGeneId(data_filter_x, from = "Ensembl.ID")

names(data_filter)[8:12631] <- aux2$HGNC.symbol


final_data_gen <- colnames(data_filter)


data_PPInR_filtrado <- data_PPInR %>%
  filter(gen %in% final_data_gen)

# k-NN model
Predictors <- as.vector(head(data_PPInR_filtrado[, 1], 100))
Predictors <- as.character(unlist(Predictors))

colnames(data_filter)[is.na(colnames(data_filter))] <- paste0("xs", seq_along(colnames(data_filter) == ""))
set.seed(1)

data_filter <- data_filter %>%
  group_by(sample_type) %>%
  sample_n(123, replace = TRUE) %>%
  ungroup()

sample.index <- sample(1:nrow(data_filter), nrow(data_filter) * 0.6, replace = FALSE)

train.data <- data_filter[sample.index, c(Predictors, "sample_type"), drop = FALSE]
test.data <- data_filter[-sample.index, c(Predictors, "sample_type"), drop = FALSE]

train.data$sample_type <- factor(train.data$sample_type)
test.data$sample_type <- factor(test.data$sample_type)

# Train the k-NN model
ctrl <- trainControl(method = "cv", p = 0.6)
knnFit <- train(sample_type ~ .,
                data = train.data,
                method = "knn",
                trControl = ctrl,
                preProcess = c("range"),  # c("center", "scale") for z-score
                tuneLength = 50)

# Plot k-NN model
plot(knnFit)

knnPredict <- predict(knnFit, newdata = test.data)

# Create the confusion matrix for k-NN
confusionMatrix(data = knnPredict, reference = test.data$sample_type)

# Linear regression

data_filter <- data_filter %>%
  mutate(sample_type = ifelse(sample_type == "Solid Tissue Normal", 1, 0))

train.data <- data_filter[sample.index, c(Predictors, "sample_type"), drop = FALSE]
test.data <- data_filter[-sample.index, c(Predictors, "sample_type"), drop = FALSE]

# Fit linear regression model
ins_model <- lm(sample_type ~ ., data = train.data)

# Summary of linear regression model
summary(ins_model)

# Train the linear regression model
train.control <- trainControl(method = "cv", number = 10)
model <- train(sample_type ~ .,
               data = train.data,
               method = "lm",
               trControl = train.control)

# Summarize the results of linear regression model
print(model)

fit <- rpart(sample_type ~ .,
             method = "anova",
             data = data_filter[, c(Predictors, "sample_type")],
             control = rpart.control(xval = 10))

# Print the decision tree
print(fit)

# Plot the decision tree
rpart.plot::rpart.plot(fit)

## First Random Forest

fit.rf <- randomForest(sample_type ~ .,
                       data = data_filter[, c(Predictors, "sample_type")])
prediction.rf <- predict(fit.rf, test.data)
table(test.data$sample_type, prediction.rf)

## Second Random Forest
fit.rf <- randomForest(sample_type ~ .,
                       data = data_filter[, c(Predictors, "sample_type")])


prediction.rf <- predict(fit.rf, test.data)
output <- data.frame(Actual = test.data$sample_type, Predicted = prediction.rf)
RMSE = sqrt(sum((output$Actual - output$Predicted)^2) / nrow(output))

print(head(output))

######vector
# Convierte la variable de respuesta a factor si no lo está
data_filter$sample_type <- as.factor(data_filter$sample_type)

# Divide los datos en conjuntos de entrenamiento y prueba

set.seed(13)
sample.index <- sample(1:nrow(data_filter), nrow(data_filter) * 0.7, replace = FALSE)
train.data <- data_filter[sample.index, c(Predictors, "sample_type"), drop = FALSE]
test.data <- data_filter[-sample.index, c(Predictors, "sample_type"), drop = FALSE]

# Realiza la búsqueda de hiperparámetros con e1071
tune_out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "linear",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

# Extrae el mejor modelo
best_model <- tune_out$best.model

# Configura el modelo SVM con un kernel lineal utilizando los mejores hiperparámetros
svm_model <- svm(sample_type ~ ., data = train.data, kernel = "linear", cost = best_model[["cost"]])

# Realiza predicciones en el conjunto de prueba
svm_predict <- predict(svm_model, newdata = test.data)

# Evalúa el rendimiento del modelo

confusionMatrix(data = svm_predict, reference = test.data$sample_type)


# Realiza la búsqueda de hiperparámetros con e1071
tune_out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "radial",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

best_model <- tune_out$best.model

svm_model <- svm(sample_type ~ ., data = train.data, kernel = "radial", cost = best_model[["cost"]])
# Realiza predicciones en el conjunto de prueba
svm_predict <- predict(svm_model, newdata = test.data)

# Evalúa el rendimiento del modelo
library(caret)
confusionMatrix(data = svm_predict, reference = test.data$sample_type)

############Segunda Parte


folder<-dirname(rstudioapi::getSourceEditorContext()$path)
parentFolder <-dirname(folder)

DataBulk <- file.path(parentFolder, "data/ExperimentsBulk.rdata")
load(DataBulk)
ls()

gen_scores <- results[["ENSG00000145675"]][["gen_scores"]]
View(gen_scores)

gen_scores2<- gen_scores%>%arrange(desc(score))
score_column <- gen_scores2$features

# Obtén la columna "features" de gen_scores2
features_column <- gen_scores2$features

# Aplica la función changeGeneId a los valores de la columna "features"
gen_scores2$features <- AMCBGeneUtils::changeGeneId(features_column, from = "Ensembl.ID")$HGNC.symbol

gen_scores2_filter <- gen_scores2 %>%
  filter(features %in% final_data_gen) #2

Predictors_2 <- head(gen_scores2_filter$features, 100)

# Convierte a caracteres si es necesario
Predictors_2 <- as.character(Predictors_2)

data_filter2 <- data_filter %>%
  group_by(sample_type) %>%
  sample_n(123, replace = TRUE) %>%
  ungroup()

sample.index <- sample(1:nrow(data_filter2), nrow(data_filter2) * 0.7, replace = FALSE)

train.data <- data_filter2[sample.index, c(Predictors_2, "sample_type"), drop = FALSE]
test.data <- data_filter2[-sample.index, c(Predictors_2, "sample_type"), drop = FALSE]

train.data$sample_type <- factor(train.data$sample_type)
test.data$sample_type <- factor(test.data$sample_type)

# Train the k-NN model
ctrl <- trainControl(method = "cv", p = 0.7)
knnFit <- train(sample_type ~ .,
                data = train.data,
                method = "knn",
                trControl = ctrl,
                preProcess = c("range"),  # c("center", "scale") for z-score
                tuneLength = 50)

# Plot k-NN model
plot(knnFit)

# Make predictions with k-NN
knnPredict <- predict(knnFit, newdata = test.data)

# Create the confusion matrix for k-NN
confusionMatrix(data = knnPredict, reference = test.data$sample_type)


# Linear regression

data_filter2 <- data_filter2  %>%
  mutate(sample_type = ifelse(sample_type == "Solid Tissue Normal", 1, 0))

train.data <- data_filter2[sample.index, c(Predictors, "sample_type"), drop = FALSE]
test.data <- data_filter2[-sample.index, c(Predictors, "sample_type"), drop = FALSE]

# Fit linear regression model
ins_model <- lm(sample_type ~ ., data = train.data)

# Summary of linear regression model
summary(ins_model)

# Train the linear regression model
train.control <- trainControl(method = "cv", number = 10)
model <- train(sample_type ~ .,
               data = train.data,
               method = "lm",
               trControl = train.control)

# Summarize the results of linear regression model
print(model)





##arboles de decision

fit <- rpart(sample_type ~ .,
             method = "anova",
             data = data_filter2[, c(Predictors, "sample_type")],
             control = rpart.control(xval = 10))

# Print the decision tree
print(fit)

# Plot the decision tree
rpart.plot::rpart.plot(fit)

###### Bosques aleatorios

fit.rf <- randomForest(sample_type ~ .,
                       data = data_filter2[, c(Predictors, "sample_type")])
prediction.rf <- predict(fit.rf, test.data)
table(test.data$sample_type, prediction.rf)


fit.rf <- randomForest(sample_type ~ .,
                       data = data_filter2[, c(Predictors, "sample_type")])


prediction.rf <- predict(fit.rf, test.data)
output <- data.frame(Actual = test.data$sample_type, Predicted = prediction.rf)
RMSE = sqrt(sum((output$Actual - output$Predicted)^2) / nrow(output))

print(head(output))

#######################################


data_filter2$sample_type <- as.factor(data_filter2$sample_type)


set.seed(123)
sample.index <- sample(1:nrow(data_filter2), nrow(data_filter2) * 0.7, replace = FALSE)
train.data <- data_filter2[sample.index, c(Predictors, "sample_type"), drop = FALSE]
test.data <- data_filter2[-sample.index, c(Predictors, "sample_type"), drop = FALSE]


tune_out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "linear",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))


best_model <- tune_out$best.model


svm_model <- svm(sample_type ~ ., data = train.data, kernel = "linear", cost = best_model[["cost"]])


svm_predict <- predict(svm_model, newdata = test.data)



confusionMatrix(data = svm_predict, reference = test.data$sample_type)

tune_out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "radial",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

best_model <- tune_out$best.model

svm_model <- svm(sample_type ~ ., data = train.data, kernel = "radial", cost = best_model[["cost"]])

svm_predict <- predict(svm_model, newdata = test.data)


confusionMatrix(data = svm_predict, reference = test.data$sample_type)

# Realiza la búsqueda de hiperparámetros con e1071
tune_out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "sigmoid",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

best_model <- tune_out$best.model

svm_model <- svm(sample_type ~ ., data = train.data, kernel = "sigmoid", cost = best_model[["cost"]])
# Realiza predicciones en el conjunto de prueba
svm_predict <- predict(svm_model, newdata = test.data)

# Evalúa el rendimiento del modelo
confusionMatrix(data = svm_predict, reference = test.data$sample_type)

