echo "Starting blastx search at:"
date
wait

blastn -query ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta \
-db ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_transcriptome_database/Nagy_Transcriptome_db \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp" \
-out ../../data-raw/010_self-blast/Nagy_transcriptome_SelfBlast.csv \
-num_threads 12 \
-max_target_seqs 2

wait
echo "All done at:"
date
