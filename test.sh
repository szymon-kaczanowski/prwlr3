$1 prwlr.py species_list  list_apoptotic_genes test    
$1 prwlr_ko.py species_list  list_apoptotic_genes test1   
diff test.tsv test/test.tsv | wc | awk '$1==0 {print "prwlr.py properly intalled"}'
diff test.tsv test/test1.tsv | wc | awk '$1==0 {print "prwlr_ko.py properly intalled"}'
