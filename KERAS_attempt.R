####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###########
###### TEST TENSOR FLOW FOR MODELLING ################
####~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~###########
### https://tensorflow.rstudio.com/install/
### https://rstudio.github.io/cheatsheets/html/keras.html?_gl=1*1gw0y35*_ga*NzExNjcxNTg5LjE2ODI5Mjk5MTE.*_ga_2C0WZ1JHG0*MTY4NjA2MzYxOC42LjAuMTY4NjA2MzYxOC4wLjAuMA..
### https://rstudio.github.io/reticulate/articles/python_packages.html?_gl=1*1y93qbc*_ga*ODcyMzQ5NjAxLjE2ODM1NTU2MDM.*_ga_2C0WZ1JHG0*MTY4NjExNjQ2NS4xMi4xLjE2ODYxMTc4NDMuMC4wLjA.
### https://tensorflow.rstudio.com/guides/keras/training_with_built_in_methods




### attempt to improve nest detection with TensorFlow
### installation was a pain - virtual environment did not work, had to install miniconda and all python libraries manually (not sure whether that needs to be repeated)
## 7 June 2023 steffen.oppel@vogelwarte.ch

# install.packages("devtools", dependencies = TRUE) 
# library(devtools)
# devtools::install_github("steffenoppel/NestTool", dependencies=TRUE) # development version - add argument 'build_vignettes = FALSE' to 
# install.packages("tensorflow")
# install.packages("keras")
library(reticulate)
# path_to_python <- install_python()
# path_to_python <- "C:/Program Files/Python/python.exe"
# virtualenv_create("r-reticulate", python = path_to_python)
# py_install("pip")
# py_install("scipy")
# py_install("requests")
# py_install("Pillow")
# py_install("pandas")
# py_install("pydot")
# py_install("h5py")
# py_install("tensorflow")
# py_install("tensorflow-hub")
# py_install("tensorflow-datasets")
library(tensorflow)
# install_tensorflow(envname = "r-reticulate")

library(keras)
#install_keras()
library(tidymodels)
library(recipes)


# ########## COPY THE BASIC EXAMPLE ###########
# 
# # input layer: use MNIST images
# mnist <- dataset_mnist()
# x_train <- mnist$train$x
# y_train <- mnist$train$y 
# x_test <- mnist$test$x
# y_test <- mnist$test$y
# 
# # reshape and rescale
# x_train <- array_reshape(x_train, c(nrow(x_train), 784)) 
# x_test <- array_reshape(x_test, c(nrow(x_test), 784)) 
# x_train <- x_train / 255
# x_test <- x_test / 255
# 
# y_train <- to_categorical(y_train, 10) 
# y_test <- to_categorical(y_test, 10)
# 
# # defining the model and layers
# model <- keras_model_sequential() 
# model %>%
#   layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>%
#   layer_dropout(rate = 0.4) %>% 
#   layer_dense(units = 128, activation = 'relu') %>% 
#   layer_dense(units = 10, activation = 'softmax')
# 
# # compile (define loss and optimizer)
# model %>%
#   compile(
#     loss = 'categorical_crossentropy', 
#     optimizer = optimizer_rmsprop(), 
#     metrics = c('accuracy')
#   )
# 
# # train (fit)
# model %>% fit(
#   x_train, y_train,
#   epochs = 30, batch_size = 128, 
#   validation_split = 0.2
# )
# 
# model %>% evaluate(x_test, y_test) 
# model %>% predict_classes(x_test)




# LOAD PREPARED TRACKING DATA
## set root folder for project
setwd("C:/Users/sop/OneDrive - Vogelwarte/REKI/Analysis/Feeding")
track_sf<-fread("data/REKI_annotated_feeding.csv")
track_nofor_day<-track_sf %>% filter(FOREST==0) %>% filter(tod_=="day")
DATA <- track_nofor_day %>%
  mutate(YDAY=yday(t_), hour=hour(t_), month=month(t_)) %>%
  filter(!is.na(step_length)) %>%
  filter(!is.na(turning_angle)) %>%
  filter(!is.na(speed)) %>%
  filter(!is.na(mean_speed)) %>%
  filter(!is.na(mean_angle)) %>%
  select(-tod_,-FOREST,-geometry,-forest_size,-build_id)
head(DATA)


DATA_TRAIN<- DATA %>% filter(bird_id %in% selectids)
DATA_TEST<- DATA %>% filter(!(bird_id %in% selectids))

table(DATA_TRAIN$FEEDER)
table(DATA_TEST$FEEDER)


############################################## FOLLOW A SIMPLE EXAMPLE ###################################
#### adapted from https://tensorflow.rstudio.com/tutorials/keras/regression



# CREATE DUMMIES FOR FACTORS ###########

DATA <- recipe(mpg ~ ., DATA) %>%
  step_num2factor(sex, levels = c("m", "f")) %>%
  step_dummy(origin, one_hot = TRUE) %>%
  prep() %>%
  bake(new_data = NULL)



split <- initial_split(DATA, 0.75)
train_dataset <- training(split)
test_dataset <- testing(split)



train_dataset %>%
  select(mpg, cylinders, displacement, weight) %>%
  GGally::ggpairs()



train_features <- train_dataset %>% select(-mpg)
test_features <- test_dataset %>% select(-mpg)

train_labels <- train_dataset %>% select(mpg)
test_labels <- test_dataset %>% select(mpg)


my_skim <- skimr::skim_with(numeric = skimr::sfl(mean, sd))
train_dataset %>%
  select(where(~is.numeric(.x))) %>%
  pivot_longer(
    cols = everything(), names_to = "variable", values_to = "values") %>%
  group_by(variable) %>%
  summarise(mean = mean(values), sd = sd(values))

normalizer <- layer_normalization(axis = -1L)

normalizer %>% adapt(as.matrix(train_features))

print(normalizer$mean)

first <- as.matrix(train_features[1,])

cat('First example:', first)

cat('Normalized:', as.matrix(normalizer(first)))


#### FIT A SIMPLE UNIVARIATE REGRESSION MODEL ####

horsepower <- matrix(train_features$horsepower)
horsepower_normalizer <- layer_normalization(input_shape = shape(1), axis = NULL)
horsepower_normalizer %>% adapt(horsepower)

horsepower_model <- keras_model_sequential() %>%
  horsepower_normalizer() %>%
  layer_dense(units = 1)

summary(horsepower_model)


predict(horsepower_model, horsepower[1:10,])

horsepower_model %>% compile(
  optimizer = optimizer_adam(learning_rate = 0.1),
  loss = 'mean_absolute_error'
)


history <- horsepower_model %>% fit(
  as.matrix(train_features$horsepower),
  as.matrix(train_labels),
  epochs = 100,
  # Suppress logging.
  verbose = 0,
  # Calculate validation results on 20% of the training data.
  validation_split = 0.2
)


