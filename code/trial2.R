setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/fasta")

library(ape)
mydna <- read.dna('sample2.fasta', as.character = 'TRUE', format = 'fasta')
ape::as.DNAbin(mydna)

# Assuming 'mydna' contains multiple sequences
# Use square brackets with the index of the sequence to subset
mydna2 <- mydna[1]
combined_dna <- c(mydna2, mydna)

#muscle
mydna_ss <- sapply(as.character(mydna), paste0, collapse="")
class(mydna_ss); length(mydna_ss);names(mydna_ss)
mydna_ss[1]

mydna_ss <- Biostrings::DNAStringSet(mydna_ss)
class(mydna_ss)
mydna_ss
mydna_ss_aligned <- muscle::muscle(mydna_ss)

load('muscle_output.RData')
mydna_ss_aligned
class(mydna_ss_aligned)


#viz
library(Biostrings)
library(ggmsa)
library(gridExtra)
devtools::install_github("YuLab-SMU/ggmsa")
mydna_ss_aligned <- readDNAStringSet("sample_aligment.fasta")
mydna_ss_aligned <- mydna_ss_aligned[c(1:7,229:235),]
DNAStringSet(mydna_ss_aligned)
#Biostrings::writeXStringSet(mydna_ss_aligned, 'my_covid_aligment.fasta')
ggmsa(mydna_ss_aligned, start = 100, end = 130, char_width = 0.5, seq_name = T,
      color='Chemistry_NT') + 
  geom_msaBar()
ggmsa(mydna_ss_aligned, start = 170, end = 200, char_width = 0.5, seq_name = T,
           color='Chemistry_NT') + 
  geom_msaBar()
a
b
ggmsa::ggmsa(mydna_ss_aligned, start=0, end=40, color='Chemistry_NT')
ggmsa::ggmsa(mydna_ss_aligned,
             start=length(mydna_ss_aligned[[1]])-40,
             end=length(mydna_ss_aligned[[1]]), color='Chemistry_NT')


mydna_ss_aligned <- ape::as.DNAbin(mydna_ss_aligned) 

#snp
library(reshape2)
mydna_mx_aligned <- as.matrix(mydna_ss_aligned)
adegenet::snpposi.plot(mydna_mx_aligned, codon=FALSE) + theme_minimal()

mysnps <- adegenet::DNAbin2genind(mydna_mx_aligned)
class(mysnps)
mysnps
snpsdf <- as.matrix(mysnps$tab)
longData <- melt(snpsdf)
ggplot(longData, aes(x=Var2, y=Var1)) + geom_raster(aes(fill=value))+
  labs(x='Loci', y='Sequence') + theme_minimal() + 
  theme(axis.text.x=element_text(size=2, angle=90, vjust = 0.3))

#phylo
setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/fasta")
mydna <- read.dna('sample_aligment.fasta', as.character = 'TRUE', format = 'fasta')
ape::as.DNAbin(mydna)
mydna <- read.FASTA('sample_aligment.fasta')
class(mydna)

names(mydna) <- sub("_(\\d{4}_\\d{2}_\\d{2})", "", names(mydna))

#gen distance
D.n <- ape::dist.dna(mydna, model = "N")
class(D.n)
D.n

D_df <- as.data.frame(as.matrix(D.n))

library(pheatmap)
pheatmap(as.matrix(D.n), fontsize_row = 1,fontsize_col = 1)
pheatmap(as.matrix(D.n), display_numbers = TRUE,
         number_color = "black", 
         fontsize_number = 2)

colnames(D_df[2,])[which(D_df[2,] == 2)]

D.raw <- ape::dist.dna(mydna, model = "raw")
D.raw
D.raw_sim <- (1-(D.raw))*100
D_df_raw <- as.matrix(D.raw_sim)
mat_df <- melt(D_df_raw)
mat_df$value <- round(mat_df$value, 3)

# Plot with ggplot
ggplot(mat_df, aes(x = Var2, y = Var1, fill = value, label = value)) +
  geom_tile(color = "white") +
  geom_text(color = "gray", size = 1.5) +
  scale_fill_gradient(low = "white", high = "blue") +
  theme_minimal() +
  labs(x = " ", y = " ", fill = "Similarities") +
  ggtitle("Matrix Similarities") + theme(axis.text.x=element_text(size=6, angle=90),
                                          axis.text.y=element_text(size=6))

library(phangorn)
mymodels <- phangorn::modelTest(phangorn::as.phyDat(mydna))
mymodels[which.min(mymodels$BIC),]


D_TN93 <- ape::dist.dna(mydna ,model = "TN93")
D_TN93


#Inferring phylogenetic trees

tre <- ape::nj(D_TN93)
tre$edge.length[tre$edge.length<0] <- 0
class(tre)
tre
plot(tre, cex=0.4); title("NJ tree of sarscov2.fas")

#tree rooted
ape::is.rooted(tre)
tre$tip.label
myroot <- which(tre$tip.label=="China01_2019_12_12")
tre_rooted <- phytools::reroot(tre, myroot)

#ape::is.rooted(tre_rooted)
plot(tre_rooted, cex=0.7)

ape::write.tree(tre_rooted, 'tre_rooted.tre')
write.nexus(tre_rooted, file = "tre_rooted.nexus")

locations <- gsub("([^A-Z])+", "", x = tre_rooted$tip.label)

split_vector <- strsplit(tre_rooted$tip.label, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)
locations <- sapply(split_vector, `[`, 1)
mypalette <- grDevices::colorRampPalette(c("turquoise", "yellow", "red", "pink", 'grey'))

plot(tre_rooted, show.tip = FALSE)
ape::tiplabels(locations, bg = adegenet::fac2col(locations, col.pal =
                                                   mypalette), cex = 0.6, fg="transparent")
title("NJ tree of xbb1.5")

library(ggtree)


ggtree(tre_rooted, colour="#99d8c9") + geom_treescale() + geom_tiplab()

ggtree(tre_rooted, colour="#99d8c9") + geom_tiplab() + theme_tree2()
ggtree(tre_rooted, colour="#99d8c9") + geom_tiplab() +
  geom_tippoint(colour="#2ca25f") + geom_nodepoint(colour="#99d8c9") +
  theme_tree2()

locations.df <- data.frame(label=tre_rooted$tip.label, locations)

# Extracting numbers between country and year
extracted_numbers <- gsub("\\D+", ",", locations.df$label) # Replace non-digits with commas
extracted_numbers <- strsplit(extracted_numbers, ",") # Split at commas

# Extract the second element from each split
extracted_numbers <- sapply(extracted_numbers, function(x) x[2])

# Convert to numeric
extracted_numbers <- as.numeric(extracted_numbers)

# Add to the data frame
locations.df$id <- extracted_numbers

ggtree(tre_rooted) %<+% locations.df +
  geom_tiplab(aes(label=id)) +
  geom_tippoint(aes(colour=locations))

ggtree(tre_rooted, layout="circular") %<+% locations.df +
  geom_tiplab(aes(label=id)) +
  geom_tippoint(aes(colour=locations))

#bootstrapping
setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/fasta")
mybootsdata <- treeio::read.iqtree("sample_aligment.fasta.contree")
myroot <- which(mybootsdata@phylo$tip.label=="China1_2019_12_12")
mybootsdata@phylo <- phytools::reroot(mybootsdata@phylo , myroot)
mybootsdata
class(mybootsdata)
str(mybootsdata@phylo)

#summarize
mean(as.numeric(mybootsdata@phylo$node.label), na.rm = T)


mybootsdata@phylo
mybootstree <- mybootsdata@phylo
locations <- gsub("([^A-Z])+", "", x = mybootstree$tip.label)
locations.df <- data.frame(label=mybootstree$tip.label ,locations)
ggtree(mybootstree) + geom_tiplab(size = 0.8)+
  geom_text2(aes(subset=!isTip, label=label),
             size = 1,
             color = "#0063B1",
             hjust = 1,
             vjust = -1.5)
ggtree(mybootstree, branch.length='none') + geom_tiplab(size = 0.8)+
  geom_text2(aes(subset=!isTip, label=label),
             size = 1.2,
             color = "#0063B1",
             hjust = 1,
             vjust = -1.5) + ggplot2::xlim(0, 50)




#tree time
setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/fasta")
mytree <- treeio::read.iqtree('sample_aligment.fasta.contree')
myroot <- which(mytree@phylo$tip.label=="China1_2019_12_12")
mytree@phylo <- phytools::reroot(mytree@phylo, myroot)
mytree

#add metadata
setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/metadata")
genom <- read.csv('metadata_covid.csv')
node <- read.csv('metadata_node.csv')
genom$kode <- sub("_.*", "", genom$id)
genom$Collection.date <- as.Date(genom$Collection.date, format = "%Y-%m-%d")
sorted_df <- arrange(genom, Collection.date, country, region)
genom$region <- ifelse(genom$region == "", genom$country, genom$region)

#alt time scale
sts <- genom$Collection.date
sts <- lubridate::as_date(sts)
sts <- lubridate::decimal_date(sts)
sts
names(sts) <- genom$id
sts

#time scales2
setwd("C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/fasta")
mydna <- ape::read.FASTA('sample_aligment.fasta')
mylength <- sapply(mydna, length)
time.data <- treedater::dater(mytree@phylo, sts=sts, s = mylength, clock ="uncorrelated",  searchRoot=TRUE)
class(time.data)
time.data

time.data$timeOfMRCA
lubridate::date_decimal(time.data$timeOfMRCA)

time.tree <- time.data$intree
time.tree$edge.length <- time.data$edge.length

split_vector <- strsplit(mytree@phylo$tip.label, "(?<=[A-Za-z])(?=[0-9])|(?<=[0-9])(?=[A-Za-z])", perl = TRUE)
locations <- sapply(split_vector, `[`, 1)


locations.df <- data.frame(label=time.tree$tip.label, locations)
# Extracting numbers between country and year
extracted_numbers <- gsub("\\D+", ",", mytree@phylo$tip.label) # Replace non-digits with commas
extracted_numbers <- strsplit(extracted_numbers, ",") # Split at commas

# Extract the second element from each split
extracted_numbers <- sapply(extracted_numbers, function(x) x[2])

# Convert to numeric
extracted_numbers <- as.numeric(extracted_numbers)

# Add to the data frame
locations.df$id <- extracted_numbers
genom <- left_join(locations.df, genom, by = c("label" = "id"))

most_recent_date <- max(lubridate::as_date(sub(".*?_", "", time.tree$tip.label)))

gg <- ggtree(time.tree, mrsd=most_recent_date, size = 0.25, color = 'black') %<+% locations.df +
  geom_tiplab(aes(label=id), size = 0.7) +
  geom_tippoint(aes(colour=locations), size = 0.7) +
  theme_tree2()
gg2 <- gg + geom_text2(aes(subset=!isTip, label=label),
                size = 0.5,
                color = "#0063B1",
                hjust = 1,
                vjust = -0.5) + 
  theme(legend.position = c(0.2, 0.7))
gg2
gg + geom_text(aes(label=node), size = 0.4, hjust = 1,
               vjust = -0.5, color = 'darkgreen')


ggsave('test.tiff',dpi=900)

#extract molecular time
mol_tim <- data.frame(id = gg$data$label, molecular_time = gg$data$x)
mol_tim$molecular_time <- as.numeric(mol_tim$molecular_time) 
genom <- left_join(genom, mol_tim, by = c("label" = "id"))
genom$molecular_time <- date_decimal(genom$molecular_time)
genom <- arrange(genom, molecular_time)


#animate

library(gganimate)
animation <- gg +
  transition_states(x, wrap = FALSE, transition_length = 3, state_length = 3) +
  shadow_mark()
animate(animation, nframes = 100, duration = 10)

ggtree(mybootstree, branch.length='none') + geom_tiplab(size = 2)+
  geom_text2(aes(subset=!isTip, label=label),
             size = 3,
             color = "#0063B1",
             hjust = 1,
             vjust = -1.5) + ggplot2::xlim(0, 50)

#MRCA(tree, tip=c("Indonesia1_2022_11_17"))

#branch
branch <- c(30, 33, 43, 44, 96, 97, 106, 145, 157, 158, 173, 202, 203, 205, 244, 247, 248, 251,
             268, 305, 311, 312, 322, 332, 337, 339, 340, 344, 345, 347, 349, 357, 359, 371, 373,
             383, 385, 392, 412, 413, 416, 417, 420, 422, 428, 440, 449, 453, 454, 455, 462)


ggtree(time.tree, mrsd=most_recent_date) %<+% genom  + 
  geom_hilight(node = branch, fill="purple") +
  theme_tree2()

gg2 + geom_hilight(node = branch, fill="purple")


#plotly

id <- time.tree$tip.label
p1 <- ggtree(time.tree, mrsd=most_recent_date)
metat <- p1$data %>%
  dplyr::inner_join(genom, c('label' = 'label')) 
metat <- metat %>% inner_join(node, c('label' = 'id'))
p2 <- p1 +
  geom_point(data = metat,
             aes(
               colour = country,
               label = id))
p3 <- p1 +
  geom_point(data = metat,
             aes(colour=country, text=paste('label:', label, 
                                           "<br>date:", Collection.date,
                                           '<br>country:', country,
                                           '<br>region:', region,
                                           '<br>city:', city,
                                           '<br>gender:', Gender,
                                           '<br>age:', Patient.age,
                                           '<br>cluster:', cluster))) +
  geom_text(data = metat,
            aes(label=id), colour="black", size=3, vjust = 0, nudge_y = 0, hjust=0, nudge_x = 0.000009) +
  labs(title = "XBB 1.5 Indonesia")
saveRDS(p3, "C:/Users/dhihr/OneDrive - London School of Hygiene and Tropical Medicine/covid/plot/p3.rds")
plotly::ggplotly(p3, tooltip = "text")




