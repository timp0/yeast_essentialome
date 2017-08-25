library(tidyverse)
library(GenomicAlignments)

workdir="/atium/Data/NGS/Aligned/170823_gale"
rawdir="/mithril/Data/NGS/Raw/170707_ANG_QiSEQ_Rend_PE_CusPrimer2Trial12pmlibrary"

##Make ref
ref="/mithril/Data/NGS/Reference/cglabrata/CGD138"
if (FALSE) {
    system(paste0("bowtie2-build --threads 8 ",
                  "/mithril/Data/NGS/Reference/cglabrata/C_glabrata_CBS138_current_chromosomes.fasta ",
                  ref))
}


bammy=file.path(workdir, "170707.sorted.bam")
##Align
if (TRUE) {
    
    read1=list.files(path=rawdir, pattern=".*R1.*fastq.gz", full.names=T)
    
    system(paste0("atropos -T 6 --discard-untrimmed -g ^TTGTTGAAGTTCTCTG -o /tmp/gale.trim.fq.gz ",
                  "-q 20 -se ", read1,
                  " >", file.path(workdir, "170707_atropos.log")))

    system(paste0("bowtie2 -p 10 -t ",
                  "-x ", ref,
                  " -U /tmp/gale.trim.fq.gz",
                  " 2> ", file.path(workdir, "170707_bowtie2.log"),
                  " | awk '$19 !~ /MD:Z:0/' | ", ##Filter the MD:Z:0 tag with awk in stream
                  "samtools view -bS - | samtools sort - -o ", bammy))
                                   



}


if (TRUE) {
    features=read_tsv(file.path(workdir, "All.gff3.txt"), comment="##",
                      col_names=c("chr", "source", "type", "start", "end", "score", "strand", "phase", "attributes")) %>%
        filter(type=="ORF") %>%
        separate(attributes, into=c("ID", "Note", "orf_classification"), sep=";", extra="drop") %>%
        mutate(ID=gsub("ID=Sequence:", "", ID))

    features.gr=GRanges(seqnames=features$chr, ranges=IRanges(start=features$start, end=features$end),
                        strand=features$strand, names=features$ID)

    ##Trim 10% off end
    perc=.1
    features.gr=resize(features.gr, width=width(features.gr)*(1-perc))

}


if (TRUE) {
    ##ok - filter bam to remove MD:Z:0
    raw=readGAlignments(bammy)

    features$num.hits=countOverlaps(features.gr, raw, type="any")

    
}



            #5'
            if orientation == "-":
                end=end-firstxth
            if orientation == "+":
                beginning=beginning+firstxth                
            if beginning <= insert and insert <= end: #checks if it isbetween that gene
                i=str(i)
                identifier=name
                if identifier in ChromcountingUnique:
                    ChromcountingUnique[identifier]+=1 
                else:
                    ChromcountingUnique.update({identifier:01}) #each gene how many times there was one in it
                Chromcountinglist.append(identifier)
    if chrom == "ChrL_C_glabrata_CBS138":  #for this one if it is in chrom I
        for i in featureschromL: #cycles through the features
            countchromLfeaturekey+=1
            insert=int(startlocation)
            beginning=int(featureschromL[i][0])
            end=int(featureschromL[i][1])
            name=featureschromL[i][3]
            orientation=featureschromL[i][2]
            lastxth=threeprimetrimarea*abs(beginning-end)
            firstxth=fiveprimetrimarea*abs(beginning-end)            
            #3'
            if orientation == "-":
                beginning=beginning+lastxth
            if orientation == "+":
                end=end-lastxth            
            #5'
            if orientation == "-":
                end=end-firstxth
            if orientation == "+":
                beginning=beginning+firstxth                
            if beginning <= insert and insert <= end: #checks if it isbetween that gene
                i=str(i)
                identifier=name
                if identifier in ChromcountingUnique:
                    ChromcountingUnique[identifier]+=1 
                else:
                    ChromcountingUnique.update({identifier:01}) #each gene how many times there was one in it
                Chromcountinglist.append(identifier)
    if chrom == "ChrM_C_glabrata_CBS138":  #for this one if it is in chrom I
        for i in featureschromM: #cycles through the features
            countchromMfeaturekey+=1
            insert=int(startlocation)
            beginning=int(featureschromM[i][0])
            end=int(featureschromM[i][1])
            name=featureschromM[i][3]
            orientation=featureschromM[i][2]
            lastxth=threeprimetrimarea*abs(beginning-end)
            firstxth=fiveprimetrimarea*abs(beginning-end)            
            #3'
            if orientation == "-":
                beginning=beginning+lastxth
            if orientation == "+":
                end=end-lastxth            
            #5'
            if orientation == "-":
                end=end-firstxth
            if orientation == "+":
                beginning=beginning+firstxth                
            if beginning <= insert and insert <= end: #checks if it isbetween that gene
                i=str(i)
                identifier=name
                if identifier in ChromcountingUnique:
                    ChromcountingUnique[identifier]+=1 
                else:
                    ChromcountingUnique.update({identifier:01}) #each gene how many times there was one in it
                Chromcountinglist.append(identifier)#




sorted_ChromcountingUnique = sorted(ChromcountingUnique.items(), key=operator.itemgetter(1))

print sorted_ChromcountingUnique





"""outputlist = open("outputoflist.txt", "w") 
outputlist.write(sorted_ChromcountingUnique)
outputlist.close
"""
count=0
for line in ChromcountingUnique:
    count+=1
print count

