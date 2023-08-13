conda activate medaka
minimap2 -x map-ont -t 50 -a assembly.fasta ../assembled_by_repeates/1/all.fastq.gz | samtools sort -o al.bam --write-index -

micromamba activate graphmb
jgi_summarize_bam_contig_depths --outputDepth depth.txt al.bam

micromamba activate checkm2 

mkdir edges
cd edges

micromamba activate utils

seqkit split ../assembly.fasta -s 999

for files in ../assembly.fasta.split/*; do 
	cat $files | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'
done

cd ..
find edges/ -name "* *" -type f | rename 's/ /_/g'

# evaluate edges
checkm taxonomy_wf -t 20 -x fa domain Bacteria edges/ checkm_edges/
checkm qa -t 20  checkm_edges/Bacteria.ms checkm_edges/ -f checkm_edges_polished_results.txt --tab_table -o 2 

graphmb --assembly . --outdir results  --cuda --tsne --post writebins --numcores 50



