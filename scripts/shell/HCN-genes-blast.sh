echo "Starting blastn search at:"
date
wait

blastn -query ../../data-raw/011_HCN-genes-blast/T-repens_HCN_genes.fasta \
-db ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_transcriptome_database/Nagy_Transcriptome_db \
-outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovs qcovhsp" \
-out ../../data-clean/011_HCN-genes-blast/HCN-genes_NagyTranscriptBlast.txt \
-num_threads 12 \
-max_target_seqs 1

wait
echo "All done at:"
date
