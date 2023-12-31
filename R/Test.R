library(conflicted)
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


datanormal <-(DynamicCancerDriverKM::BRCA_normal)
dataPt <-(DynamicCancerDriverKM::BRCA_PT)
final_data <- bind_rows(datanormal, dataPt)


porcentaje_menor_10 <- final_data %>%
  summarise_all(~ mean(. <400, na.rm = TRUE))


columnas_a_eliminar <- names(porcentaje_menor_10[, porcentaje_menor_10 >= 0.8])


final_data_filtrado <- final_data %>%
  select(-one_of(columnas_a_eliminar))

final_data_filtrado2 <- final_data_filtrado

data_pii<-(DynamicCancerDriverKM::PPI)

data_piin <- data_pii %>%
  pivot_longer(cols = c(`Input-node Gene Symbol`, `Output-node Gene Symbol`), names_to = "variable", values_to = "gen") %>%
  group_by(gen, variable) %>%
  summarise(frecuencia = n()) %>%
  pivot_wider(names_from = variable, values_from = frecuencia, values_fill = 0)

data_piinR <- data_piin %>%
  mutate(total_mode = `Input-node Gene Symbol` + `Output-node Gene Symbol`) %>%
  select(total_mode) %>%
  arrange(desc(total_mode))


print(data_piinR)

final_data_filtradox<-colnames(final_data_filtrado)[ 8:ncol(final_data_filtrado)]
aux2 <- AMCBGeneUtils::changeGeneId(final_data_filtradox, from = "Ensembl.ID")

names(final_data_filtrado)[8:12631] <- aux2$HGNC.symbol


genes_en_final_data <- colnames(final_data_filtrado)

<<<<<<< HEAD
conflicts_prefer(dplyr::filter)
data_PPInR_filtrado <- data_PPInR %>%
  filter(gen %in% final_data_gen)
=======

data_piinR_filtrado <- data_piinR %>%
  filter(gen %in% genes_en_final_data)
>>>>>>> 94b89f9010ea4d337d4f45b6b69746919962aa85

# k-NN model
Predictores <- as.vector(head(data_piinR_filtrado[, 1], 100))
Predictores <- as.character(unlist(Predictores))

<<<<<<< HEAD
colnames(data_filter)[is.na(colnames(data_filter))] <- paste0("xs", seq_along(colnames(data_filter) == ""))
set.seed(13)
=======
colnames(final_data_filtrado)[is.na(colnames(final_data_filtrado))] <- paste0("xs", seq_along(colnames(final_data_filtrado) == ""))
set.seed(1)
>>>>>>> 94b89f9010ea4d337d4f45b6b69746919962aa85

final_data_filtradoe <- final_data_filtrado %>%
  group_by(sample_type) %>%
  sample_n(123, replace = TRUE) %>%
  ungroup()

sample.index <- sample(1:nrow(final_data_filtradoe), nrow(final_data_filtradoe) * 0.6, replace = FALSE)

train.data <- final_data_filtradoe[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe[-sample.index, c(Predictores, "sample_type"), drop = FALSE]

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

final_data_filtradoe <- final_data_filtradoe %>%
  mutate(sample_type = ifelse(sample_type == "Solid Tissue Normal", 1, 0))

train.data <- final_data_filtradoe[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe[-sample.index, c(Predictores, "sample_type"), drop = FALSE]

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
             data = final_data_filtradoe[, c(Predictores, "sample_type")],
             control = rpart.control(xval = 10))

# Print the decision tree
print(fit)

# Plot the decision tree
rpart.plot::rpart.plot(fit)

##primer bosque

fit.rf <- randomForest(sample_type ~ .,
                       data = final_data_filtradoe[, c(Predictores, "sample_type")])
prediction.rf <- predict(fit.rf, test.data)
table(test.data$sample_type, prediction.rf)

## segundo bosque
fit.rf <- randomForest(sample_type ~ .,
                       data = final_data_filtradoe[, c(Predictores, "sample_type")])


prediction.rf <- predict(fit.rf, test.data)
output <- data.frame(Actual = test.data$sample_type, Predicted = prediction.rf)
RMSE = sqrt(sum((output$Actual - output$Predicted)^2) / nrow(output))

print(head(output))
######vector

# Convierte la variable de respuesta a factor si no lo está
final_data_filtradoe$sample_type <- as.factor(final_data_filtradoe$sample_type)

# Divide los datos en conjuntos de entrenamiento y prueba
set.seed(1)
sample.index <- sample(1:nrow(final_data_filtradoe), nrow(final_data_filtradoe) * 0.7, replace = FALSE)
train.data <- final_data_filtradoe[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe[-sample.index, c(Predictores, "sample_type"), drop = FALSE]

# Realiza la búsqueda de hiperparámetros con e1071
tune.out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "linear",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

# Extrae el mejor modelo
bestmod <- tune.out$best.model

# Configura el modelo SVM con un kernel lineal utilizando los mejores hiperparámetros
svm_model <- svm(sample_type ~ ., data = train.data, kernel = "linear", cost = bestmod[["cost"]])

# Realiza predicciones en el conjunto de prueba
svm_predict <- predict(svm_model, newdata = test.data)

# Evalúa el rendimiento del modelo

confusionMatrix(data = svm_predict, reference = test.data$sample_type)


# Realiza la búsqueda de hiperparámetros con e1071
tune.out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "radial",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

bestmod <- tune.out$best.model

svm_model <- svm(sample_type ~ ., data = train.data, kernel = "radial", cost = bestmod[["cost"]])
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

<<<<<<< HEAD
gen_scores <- results[["ENSG00000145675"]][["geneScore"]]
View(gen_scores)

gen_scores2<- gen_scores%>% arrange(desc(score))
score_column <- gen_scores2$features
=======
geneScore <- results[["ENSG00000145675"]][["geneScore"]]

View(geneScore)
>>>>>>> 94b89f9010ea4d337d4f45b6b69746919962aa85

geneScore2<- geneScore%>%arrange(desc(score))

score_column <- geneScore2$features

# Obtén la columna "features" de geneScore2
features_column <- geneScore2$features

# Aplica la función changeGeneId a los valores de la columna "features"
geneScore2$features <- AMCBGeneUtils::changeGeneId(features_column, from = "Ensembl.ID")$HGNC.symbol

geneScore2_filtrado <- geneScore2 %>%
  filter(features %in% genes_en_final_data) #2

Predictores_2 <- head(geneScore2_filtrado$features, 100)

# Convierte a caracteres si es necesario
Predictores_2 <- as.character(Predictores_2)

final_data_filtradoe2 <- final_data_filtrado %>%
  group_by(sample_type) %>%
  sample_n(123, replace = TRUE) %>%
  ungroup()

sample.index <- sample(1:nrow(final_data_filtradoe2), nrow(final_data_filtradoe2) * 0.7, replace = FALSE)

train.data <- final_data_filtradoe2[sample.index, c(Predictores_2, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe2[-sample.index, c(Predictores_2, "sample_type"), drop = FALSE]

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

final_data_filtradoe2 <- final_data_filtradoe2  %>%
  mutate(sample_type = ifelse(sample_type == "Solid Tissue Normal", 1, 0))

train.data <- final_data_filtradoe2[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe2[-sample.index, c(Predictores, "sample_type"), drop = FALSE]

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
             data = final_data_filtradoe2[, c(Predictores, "sample_type")],
             control = rpart.control(xval = 10))

# Print the decision tree
print(fit)

# Plot the decision tree
rpart.plot::rpart.plot(fit)

###### Bosques aleatorios

fit.rf <- randomForest(sample_type ~ .,
                       data = final_data_filtradoe2[, c(Predictores, "sample_type")])
prediction.rf <- predict(fit.rf, test.data)
table(test.data$sample_type, prediction.rf)


fit.rf <- randomForest(sample_type ~ .,
                       data = final_data_filtradoe2[, c(Predictores, "sample_type")])


prediction.rf <- predict(fit.rf, test.data)
output <- data.frame(Actual = test.data$sample_type, Predicted = prediction.rf)
RMSE = sqrt(sum((output$Actual - output$Predicted)^2) / nrow(output))

print(head(output))

#######################################


final_data_filtradoe2$sample_type <- as.factor(final_data_filtradoe2$sample_type)


set.seed(123)
sample.index <- sample(1:nrow(final_data_filtradoe2), nrow(final_data_filtradoe2) * 0.7, replace = FALSE)
train.data <- final_data_filtradoe2[sample.index, c(Predictores, "sample_type"), drop = FALSE]
test.data <- final_data_filtradoe2[-sample.index, c(Predictores, "sample_type"), drop = FALSE]


tune.out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "linear",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))


bestmod <- tune.out$best.model


svm_model <- svm(sample_type ~ ., data = train.data, kernel = "linear", cost = bestmod[["cost"]])


svm_predict <- predict(svm_model, newdata = test.data)



confusionMatrix(data = svm_predict, reference = test.data$sample_type)

tune.out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "radial",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

bestmod <- tune.out$best.model

svm_model <- svm(sample_type ~ ., data = train.data, kernel = "radial", cost = bestmod[["cost"]])

svm_predict <- predict(svm_model, newdata = test.data)


confusionMatrix(data = svm_predict, reference = test.data$sample_type)

# Realiza la búsqueda de hiperparámetros con e1071
tune.out <- tune(svm,
                 sample_type ~ .,
                 data = train.data,
                 kernel = "sigmoid",
                 ranges = list(cost = c(0.001, 0.01, 0.1, 1, 5, 10, 100)))

bestmod <- tune.out$best.model

svm_model <- svm(sample_type ~ ., data = train.data, kernel = "sigmoid", cost = bestmod[["cost"]])
# Realiza predicciones en el conjunto de prueba
svm_predict <- predict(svm_model, newdata = test.data)

# Evalúa el rendimiento del modelo
confusionMatrix(data = svm_predict, reference = test.data$sample_type)

