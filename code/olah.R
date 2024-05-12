setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/metadata")
library(tidyverse)


x0 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(2).tsv', sep = '\t', header = TRUE)
x0 <- x0[,-c(17)]
x0$Patient.age <- as.numeric(x0$Patient.age)
head(x0)

x <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(3).tsv', sep = '\t', header = TRUE)
x <- x[,-c(17)]
head(x)

x2 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(4).tsv', sep = '\t', header = TRUE)
x2 <- x2[,-c(17)]
head(x2)

x3 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(5).tsv', sep = '\t', header = TRUE)
x3 <- x3[,-c(17)]
head(x3)

x4 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(6).tsv', sep = '\t', header = TRUE)
x4 <- x4[,-c(17)]
x4$Patient.age <- as.numeric(x4$Patient.age)
head(x4)

x5 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(7).tsv', sep = '\t', header = TRUE)
x5 <- x5[,-c(17)]
head(x5)

x6 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(8).tsv', sep = '\t', header = TRUE)
x6 <- x6[,-c(17)]
head(x6)

x7 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(9).tsv', sep = '\t', header = TRUE)
x7 <- x7[,-c(17)]
head(x7)

x8 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(10).tsv', sep = '\t', header = TRUE)
x8 <- x8[,-c(17)]
head(x8)
x8$Patient.age <- ifelse(x8$Patient.age == 'unknown', NA, x8$Patient.age)

x9 <- read.table(file = 'gisaid_hcov-19_2024_04_23_23(11).tsv', sep = '\t', header = TRUE)
x9 <- x9[,-c(17)]
head(x9)

bind <- bind_rows(x0, x, x2, x3, x4, x5, x6, x7, x8, x9)
genom <- bind

genom$id <- paste(genom$Virus.name, genom$Accession.ID, genom$Collection.date, sep = "|")
genom <- cbind(id = genom$id, genom[, -ncol(genom)])
genom <- genom[, !names(genom) %in% c("Virus.name", "Accession.ID")]

library(stringr)
x <- str_split_fixed(genom$Location, " / ", 4)
genom$country <- x[,2]
genom$region <- x[,3]
genom$city <- x[,4]
genom <- subset(genom, select = -Location)
genom <- genom %>% distinct(id, .keep_all = TRUE)
genom <- genom[, c('id','Collection.date', 'country', 'region', 'city', 'Patient.age', 'Gender', 'Patient.status','Last.vaccinated','Lineage')]

# Extract country from each string
country <- gsub("^hCoV-19/([^/]+)/.*", "\\1", genom$id)
# Generate sequential numbers
sequential_number <- ave(seq_along(country), country, FUN = seq_along)

# Concatenate country and sequential number
country <- paste(country, sequential_number, sep = "")
date <- sub(".+\\|", "", genom$id) %>%
  gsub("-", "_", .)
genom$id <- paste(country, date, sep = "_")

#add china
new_row <- c("China1_2019_12_12", "2019-12-12", "China", "Wuhan", NA, NA, "Male", "Hospitalized", NA, "Original")
genom <- rbind(genom, new_row)

write.csv(genom,'metadata_covid.csv', row.names = FALSE)



#sankey plot
library(ggsankey)
setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine")
metadata <- read.csv("metadata_covid_plus.csv")

sinkey <- data.frame(Gender = metadata$Gender, Age = metadata$age.group, Country = metadata$country.group)
df <- sinkey %>% make_long(Gender, Age, Country)

# 2sankey
reagg <- df%>%
  dplyr::group_by(node)%>%  # Here we are grouping the data by node and then we are taking the frequency of it 
  tally()

df2 <- merge(df, 
             reagg, 
             by.x = 'node', 
             by.y = 'node', 
             all.x = TRUE)

pl <- ggplot(df2, aes(x = x,                        
                      next_x = next_x,                                     
                      node = node,
                      next_node = next_node,        
                      fill = factor(node),
                      label = paste0(node, " : ", n)))             # This Creates a label for each node

pl <- pl +geom_sankey(flow.alpha = 0.3,          #This Creates the transparency of your node 
                      node.color = "black",     # This is your node color        
                      show.legend = TRUE)        # This determines if you want your legend to show

pl <- pl + geom_sankey_label(size = 5,
                             color = "black", 
                             fill = "white") # This specifies the Label format for each node 


pl <- pl + theme_minimal()
pl <- pl + theme(legend.position = 'none')
pl <- pl + theme(axis.title = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks = element_blank(),
                 panel.grid = element_blank())


pl <- pl + scale_fill_manual(values = c('National'    = "red",
                                        '19-50'    = "red",
                                        'Female'    = "red")) + theme(axis.text.x = element_text(size = 10))
pl <- pl + labs(title = "Sankey Diagram")
pl <- pl + labs(subtitle = "Distribution of Demographic Variables")
pl <- pl + labs(fill = 'Nodes')
pl


