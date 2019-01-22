library(keras)
library(tidyverse)

times <- read_csv("speed_est.csv") %>%
  na.omit()

# Split Datasets
train <- times %>%
  sample_frac(0.9)

val <- times %>%
  anti_join(train)

X_train <- train[, 1:3] %>%
  as.matrix()

X_val <- val[, 1:3] %>%
  as.matrix()

Y_train <- train[, 4, drop = F] %>%
  as.matrix()
Y_val <- val[,4, drop = F] %>%
  as.matrix()

model <- keras_model_sequential()
model %>% 
  layer_dense(units = 64, input_shape = c(3)) %>% 
  layer_activation('relu') %>% 
  layer_dense(units = 32, input_shape = c(3)) %>% 
  layer_activation('relu') %>% 
  layer_dense(units = 1) %>% 
  layer_activation('linear')

model %>% compile(
  optimizer = optimizer_adam(),
  loss = 'mse'
)

callbacks_list <- list(
  callback_reduce_lr_on_plateau(
    monitor = "val_loss",
    factor = 0.1,
    patience = 10
  ),
  callback_early_stopping(
    monitor = "val_loss",
    patience = 15
  ),
  callback_model_checkpoint(
    filepath = "time_model.h5",
    monitor = "val_loss",
    save_best_only = TRUE
  )
)

history <- model %>% fit(X_train,
  Y_train,
  epochs = 1000,
  batch_size = 2^6,
  validation_data = list(X_val, Y_val),
  callbacks = callbacks_list)

