library(tidyverse)

####
#Create SNP Matrix
#####
SNP<-read_tsv("GATKOut/SNP_matrix.012", col_names=F)%>%
  select(-X1)

SNP[SNP==-1]<-0

sampleIDs<-read_tsv("GATKOut/SNP_matrix.012.indv", col_names=F)%>%
  pull(X1)
SNPIDs<-read_tsv("GATKOut/SNP_matrix.012.pos", col_names=F)%>%
  unite(col = "SNPID", X1, X2, sep="_")#%>%

#Some SNPs have duplicates for some reason. I assume this is because
# some SNPs have multiple alternative alleles.
dupSNPs<-group_by(SNPIDs, SNPID)%>%
  summarize(Count=n())%>%
  filter(Count>1)%>%
  pull(SNPID)
SNPIDs<-pull(SNPIDs, SNPID)
colnames(SNP)<-SNPIDs

SNP<-select(SNP, -dupSNPs) #Remove duplicated SNPs

#Because of the grouping structure of the data I don't think a
# simple presence type filter is approprite. Instead I will just remove
# SNPs that are all alternative alleles, and any SNP with more than one alternative
# allel.
SNP_Mat<-mutate(SNP,SampleID=sampleIDs)%>%
  select(SampleID, colnames(SNP))%>%
  column_to_rownames("SampleID")%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("SNPID")%>%
  mutate(rowSum=rowSums(select(., -SNPID)))%>%
  filter( rowSum < 56)%>% #Removes rows that are all alternative alleles as these will likely be errors in the assembly.
  select(-rowSum)

#Remove SNPs with more than one alternative allele. 
SNP_Mat_notwos<-SNP_Mat[apply(SNP_Mat, 1, function(x) all(x!=2)),]
  
#Calculate MAF of all SNPS:
SNP_MAF<-SNP_Mat_notwos%>%
  mutate(MAF=rowSums(select(.,-SNPID))*100/56,
         MAF = ifelse(MAF > 50, 100-MAF, MAF))

SNP_MAF.filtered<-filter(SNP_MAF, MAF > 5)%>%
  select(-MAF)
SNP_Mat_notwos<-SNP_MAF.filtered
write_tsv(SNP_Mat_notwos, "SNP_Matrix.txt")


######
#Create SNP Location File
#####
vcf<-read_tsv("GATKOut/Merged.vcf", comment="#", 
              col_names=c("CHROM",	"POS",	"ID",	"REF",	"ALT",	"QUAL",	
                          "FILTER",	"INFO",	"FORMAT",	"SRR13275195",	
                          "SRR13275196",	"SRR13275197",	"SRR13275198",	
                          "SRR13275199",	"SRR13275200",	"SRR13275201",	
                          "SRR13275202",	"SRR13275203",	"SRR13275204",	
                          "SRR13275205",	"SRR13275206",	"SRR13275207",	
                          "SRR13275208",	"SRR13275209",	"SRR13275210",	
                          "SRR13275211",	"SRR13275212",	"SRR13275213",	
                          "SRR13275214",	"SRR13275215",	"SRR13275216",	
                          "SRR13275217",	"SRR13275218",	"SRR13275219",	
                          "SRR13275220",	"SRR13275221",	"SRR13275222",	
                          "SRR13275223",	"SRR13275224",	"SRR13275225",	
                          "SRR13275226",	"SRR13275227",	"SRR13275228",	
                          "SRR13275229",	"SRR13275230",	"SRR13275231",	
                          "SRR13275232",	"SRR13275233",	"SRR13275234",	
                          "SRR13275235",	"SRR13275236",	"SRR13275237",	
                          "SRR13275238",	"SRR13275239",	"SRR13275240",	
                          "SRR13275241",	"SRR13275242",	"SRR13275243",	
                          "SRR13275244",	"SRR13275245",	"SRR13275246",	
                          "SRR13275247", "SRR13275248", "SRR13275249", 
                          "SRR13275250"))%>%
  select(CHROM, POS)%>%
  mutate(SNPID=paste(CHROM, POS, sep="_"))%>%
  select(SNPID, CHROM, POS)%>%
  filter(SNPID %in% pull(SNP_Mat_notwos, SNPID))%>%
  distinct()

write_tsv(vcf, "SNP_Posistions.txt")


#####
#Create Gene Location File
#####

genes<-read_tsv("gffID_to_GeneID.txt")%>%
  select(GeneID, Start, End)%>%
  mutate(Chromosome="NZ_CP018664.1")%>%
  select(GeneID, Chromosome, Start, End)

write_tsv(genes, "Gene_Posistions.txt")

######
# Create a covariates file. 
# Samples as columns covariates as rows
######

sampleData<-read_tsv("SampleTable.txt")

covar<-sampleData%>%
  select(-`Sample Name`, -FileName)%>%
  mutate(Antibiotic=case_when(
    Antibiotic == "NDC" ~ 0,
    Antibiotic == "MOX" ~ 1,
    Antibiotic == "CIP" ~ 2,
    Antibiotic == "COL" ~ 3,
    Antibiotic == "TOB" ~ 4,
    Antibiotic == "IMI" ~ 5,
    Antibiotic == "MER" ~ 6,
    TRUE ~ -1
  ))%>%
  column_to_rownames("Run")%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column("CoVariate")%>%
  write_tsv("CovariateMatrix.txt")
