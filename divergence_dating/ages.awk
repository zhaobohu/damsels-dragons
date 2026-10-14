# shared: posterior node ages from an MCMCtree out.txt, tab separated
/^t_n/ { gsub(/[(),]/,""); print d"\t"$1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7 }
