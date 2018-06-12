echo "Starting blastx search at:"
date
wait

blastx -query ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta \
-db ../../../data-reference/protein-databases/$1/$1_ProteinDatabase \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp" \
-out ../../data-raw/020_blastx/$1/Nagy_transcriptome_ProteinBlast.csv \
-num_threads 12 \
-max_target_seqs 1

wait
echo "All done at:"
date
