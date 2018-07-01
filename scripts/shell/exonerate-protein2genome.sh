echo "Starting exonerate runs at:"
date
wait

# $1 is the full path to the file containing the protein sequences for the species of interest.
# For T. subterraneum and M. truncatula, these are contained in data-raw/021_retrieve-proteins whereas
# for T. pratense these are contained in data-clean/021_retrieve-proteins since the T. pratense
# sequences needed to be cleaned.

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 1 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk1.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 2 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk2.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 3 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk3.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 4 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk4.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 5 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk5.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 6 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk6.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 7 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk7.output &

exonerate --query $1 --target ../../../data-reference/reference-genomes/$2/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 8 --targetchunktotal 8 > ../../data-raw/030_exonerate_p2g/$2/Nagy_transcriptome_p2g_chunk8.output &

wait
echo "All done at:"
date
