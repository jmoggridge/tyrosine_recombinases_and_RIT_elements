# compress models
library(readr)
library(tidymodels)

knn <- read_rds('models/knn_classifier.rds')
write_rds(knn, 'models/knn.rds.gz', compress = 'gz')

rf <- read_rds('models/rf_classifier.rds')
write_rds(rf, 'models/rf.rds.gz', compress = 'gz')

glmnet <- read_rds('models/glmnet_classifier.rds')
write_rds(glmnet, 'models/glmnet.rds.gz', compress = 'gz')

install.packages('butcher')
library(butcher)

butcher::weigh(knn) |> print(n = 100)

knn_axe_env <- butcher::axe_env(knn)
knn_axe_data <- butcher::axe_data(knn_axe_env)
object.size(knn)
object.size(knn_axe_env)
object.size(knn_axe_data)


library(rlang)
env_print(knn$terms)


## Test predict method









