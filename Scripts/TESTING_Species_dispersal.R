


test <- read.csv("../Data/Species_data/TOM_TRAVERS_animal_dispersal_(acknowledgement).csv")


na.omit(test$mean_movement) %>%
  hist(.)
