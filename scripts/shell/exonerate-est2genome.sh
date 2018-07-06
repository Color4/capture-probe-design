echo "Starting exonerate runs at:"
date
wait

# $1 is the species for which the coding sequences should be mapped to the T. repens transcriptome

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 1 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk1.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 2 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk2.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 3 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk3.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 4 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk4.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 5 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk5.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 6 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk6.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95 --exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 7 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk7.output &

exonerate --query ../../data-clean/031_coding-sequences/$1/*.fasta --target ../../../data-reference/transcriptomes/Trifolium-repens/Nagy/Nagy_Transcriptome.fasta -m est2genome --showtargetgff yes --showalignment no --fsmmemory 1024 --percent 95--exhaustive --bestn 1 --ryo ">%qi length=%ql alnlen=%qal\n>%ti length=%tl alnlen=%tal\n" --targetchunkid 8 --targetchunktotal 8 > ../../data-raw/040_exonerate_e2g/$1/Nagy_transcriptome_e2g_chunk8.output &

wait
echo "All done at:"
date
