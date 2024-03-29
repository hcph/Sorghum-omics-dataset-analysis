#####MNase-seq analysis
#raw reads filtering and clean reads mappping are same to ChIP-seq
#produce input bed file
cd /public/home/chaohe/sorghum/mnase/align
cut -f 2,3 ../meta/macs2.meta | while read i j; do
bedtools bamtobed -i ../align/"$j".final.bam > "$j".final.bed
done
#usng iNPS to identify the postion of nucleosome
conda activate py3
mkdir ../nucleosome 
cut -f 3 ../meta/macs2.meta | while read j; do
python3 /public/home/chaohe/sorghum/mnase/iNPS_V1.2.2.py -i "$j".final.bed -o ../nucleosome/"$j" --s_p p
done
#merge the bed files
cd /public/home/chaohe/sorghum/mnase/nucleosome
cat CRMN1_*.like_bed >CRMN1.like.bed
#draw metaplot
#CRMN
cd /public/home/chaohe/sorghum/mnase/nucleosome
cat CRMN_reo0*.like_wig | awk -F'[\t]' '{print $1"\t"$6}' | sed '/track/d;/Coordinate/d' >CRMN.like.wig
perl -p -i -e 's/chromosome/chrom/g' CRMN.like.wig
wigToBigWig  -clip CRMN.like.wig /public/home/chaohe/db/chromInfo.bed CRMN.like.bw
bsub  -J 242-1 -n 8 -o %J.3.out -e %J.3.err  -q smp -R "rusage[mem=80GB]" "computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R CR_Q1.txt CR_Q2.txt CR_Q3.txt CR_Q4.txt --binSize 20  -S CRMN.like.bw --skipZeros --averageTypeBins mean -o CRMN_nucleosome.gz  --outFileSortedRegions CRMN_MN_genes.bed"
plotProfile --dpi 720 -m CRMN_nucleosome.gz -out CRMN_nucleosome.pdf --plotFileFormat pdf
#CLMN
cd /public/home/chaohe/sorghum/mnase/nucleosome
cat CLMN_reo0*.like_wig | awk -F'[\t]' '{print $1"\t"$6}' | sed '/track/d;/Coordinate/d' >CLMN.like.wig
perl -p -i -e 's/chromosome/chrom/g' CLMN.like.wig
wigToBigWig  -clip CLMN.like.wig /public/home/chaohe/db/chromInfo.bed CLMN.like.bw
bsub  -J 242-1 -n 8 -o %J.3.out -e %J.3.err  -q smp -R "rusage[mem=80GB]" "computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R CR_Q1.txt CR_Q2.txt CR_Q3.txt CR_Q4.txt --binSize 20  -S CLMN.like.bw --skipZeros --averageTypeBins mean -o CLMN_nucleosome.gz  --outFileSortedRegions CLMN_MN_genes.bed"
plotProfile --dpi 720 -m CLMN_nucleosome.gz -out CLMN_nucleosome.pdf --plotFileFormat pdf
#PRMN
cd /public/home/chaohe/sorghum/mnase/nucleosome
cat PRMN_reo0*.like_wig | awk -F'[\t]' '{print $1"\t"$6}' | sed '/track/d;/Coordinate/d' >PRMN.like.wig
perl -p -i -e 's/chromosome/chrom/g' PRMN.like.wig
wigToBigWig  -clip PRMN.like.wig /public/home/chaohe/db/chromInfo.bed PRMN.like.bw
bsub  -J 242-1 -n 8 -o %J.3.out -e %J.3.err  -q smp -R "rusage[mem=80GB]" "computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R CR_Q1.txt CR_Q2.txt CR_Q3.txt CR_Q4.txt --binSize 20  -S PRMN.like.bw --skipZeros --averageTypeBins mean -o PRMN_nucleosome.gz  --outFileSortedRegions PRMN_MN_genes.bed"
plotProfile --dpi 720 -m PRMN_nucleosome.gz -out PRMN_nucleosome.pdf --plotFileFormat pdf
#PLMN
cd /public/home/chaohe/sorghum/mnase/nucleosome
cat PLMN_reo0*.like_wig | awk -F'[\t]' '{print $1"\t"$6}' | sed '/track/d;/Coordinate/d' >PLMN.like.wig
perl -p -i -e 's/chromosome/chrom/g' PLMN.like.wig
wigToBigWig  -clip PLMN.like.wig /public/home/chaohe/db/chromInfo.bed PLMN.like.bw
bsub  -J 242-1 -n 8 -o %J.3.out -e %J.3.err  -q smp -R "rusage[mem=80GB]" "computeMatrix reference-point --referencePoint TSS -b 1000 -a 1000 -R CR_Q1.txt CR_Q2.txt CR_Q3.txt CR_Q4.txt --binSize 20  -S PLMN.like.bw --skipZeros --averageTypeBins mean -o PLMN_nucleosome.gz  --outFileSortedRegions PLMN_MN_genes.bed"
plotProfile --dpi 720 -m PLMN_nucleosome.gz -out PLMN_nucleosome.pdf --plotFileFormat pdf

