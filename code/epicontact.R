#epicontacts
library(outbreaks)
library(projections)
library(ggplot2)
library(lubridate)
library(dplyr)
library(tidyverse)
library(epicontacts)

setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/metadata")
linelist <- read.csv('epicontact.csv')
linelist$Gender <- ifelse(linelist$Gender == "Male", "m", "f")
linelist$Collection.date <- as.Date(linelist$Collection.date, format= "%m/%d/%Y")

contacts <- linelist %>%
  transmute(
    infector = infector,
    id = id,
    type = type,
  ) %>%
  drop_na(infector)

contacts <- filter(contacts, infector != "")

## generate epicontacts object
epic <- make_epicontacts(
  linelist = linelist,
  contacts = contacts,
  id = "id",
  from = "infector",
  to = "id",
  directed = TRUE
)

epic
## plot epicontacts object
plot(
  epic,
  width = 700,
  height = 700
)


plot(
  epic, 
  node_color = "type",
  node_shape = "Gender",
  node_size = "Patient.age",
  col_pal = c(travel = "firebrick", `non-travel` = "green"),
  shapes = c(f = "female", m = "male"),
  size_range = c(40, 60),
  height = 700,
  width = 700
)

#time-series
plot(
  epic,
  x_axis = "Collection.date",
  node_color = "type",
  col_pal = c(travel = "firebrick", `non-travel` = "green"),
  arrow_size = 0.03,
  node_size = 5,
  label = FALSE,
  height = 700,
  width = 700
)

#flow plot
flow <- filter(epic$linelist, type == "travel", infector != "") %>%
  group_by(city, Infector.city) %>%
  summarize(counts = n()) %>%
  ungroup() %>%
  arrange(desc(counts))

library(parsetR) 

parset(flow, dimensions = c('Infector.city', 'city'), 
       value = htmlwidgets::JS("function(d){return d.counts}"), 
       tension = 0.7)

