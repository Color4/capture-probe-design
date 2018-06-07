echo "Starting exonerate runs at:"
date
wait

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 1 --targetchunktotal 8 > $2/$2_protein2genome_chunk1.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 2 --targetchunktotal 8 > $2/$2_protein2genome_chunk2.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 3 --targetchunktotal 8 > $2/$2_protein2genome_chunk3.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 4 --targetchunktotal 8 > $2/$2_protein2genome_chunk4.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 5 --targetchunktotal 8 > $2/$2_protein2genome_chunk5.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 6 --targetchunktotal 8 > $2/$2_protein2genome_chunk6.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 7 --targetchunktotal 8 > $2/$2_protein2genome_chunk7.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 8 --targetchunktotal 8 > $2/$2_protein2genome_chunk8.output &

exonerate --query $1 --target ../../../../$2/Genome/*.fasta -m protein2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --geneseed 250 --percent 95 --bestn 1 --softmasktarget yes --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 1 --targetchunktotal 8 > $2/$2_protein2genome_chunk4.output &

wait
echo "All done at:"
date
