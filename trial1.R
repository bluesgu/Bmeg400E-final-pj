BiocManager::install("Rbowtie2")
BiocManager::install("msa")
BiocManager::install("DECIPHER")
BiocManager::install("seqinr")
install.packages("seqinr")
install.packages("adegenet")
install.packages("viridis")
install.packages("ape")
install.packages("remotes")
remotes::install_github("nvelden/NGLVieweR")
install.packages("bio3d", dependencies=TRUE)
library(bio3d)
library(NGLVieweR)
library(remotes)
library(viridis)
library(msa)
library(ggmsa)
library(Biostrings)
library(DECIPHER)
library(seqinr)
library(ape)
library(ade4)
library(ggmsa)

fasta_seq <- readDNAStringSet("summary.fasta",format = "fasta")
#Read the raw data of nulcited sequence from five betacoronavirus samples in fasta format. 
aligned <- AlignSeqs(fasta_seq)
#align those 5 sequence together. 

BrowseSeqs(aligned,htmlFile="E:/BMEG 400E",highlight = 0)
#view the alignment through html page

writeXStringSet(aligned,file="Cornavirus_aligned.fasta")
#save the aligned fasta file

seq <- read.alignment("Cornavirus_aligned.fasta",format = "fasta")
#read the aligned file and store it into a 

seq$nam <- c("Bat_RaTH13","Covid-19 ref seq" , "Bat Btrs", "Bat Sars like", "Sars Urbani")
dist <- dist.alignment(seq,matrix = "similarity")
#calculate the similarity between each sequence

t <- as.data.frame(as.matrix(dist))

table.paint(t)
#create a heat shows the similarity

tree <- nj(dist)
tree<-ladderize(tree) 
plot(tree)
#Generate the evolutionary tree

#align the receptor binding region 
spike_aa <- readAAStringSet("spike.fasta", format = "fasta")
# read the amino acid sequences of the spike proteins on five virus
aligned_spike <- AlignSeqs(spike_aa)
#do alignment for amino acid sequence
ggmsa(aligned_spike,start =450, end = 500,seq_name = TRUE)
#visualize the alignment 

#visualize the spike protein of SARS-coV-2 Reference sequence and  the closest related bat coronavirus. 
NGLVieweR("7BNM") %>%
  addRepresentation("spacefill")
#visualize the spike protein on SARS-coV-2
NGLVieweR("6ZGF") %>%
  addRepresentation("spacefill")
#visualize the spike protein on bat_RaTH13
