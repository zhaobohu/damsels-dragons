#!/bin/bash
# t_n15 ages of the dS subsets and their controls, with permutation p = (1 + k) / (n + 1) per fraction
set -eo pipefail     # no -x: a trace over ~3000 run directories is unreadable

W=/path/to/project/divergence_dating
cd $W
mkdir -p modelsel
MAN=controls_manifest.tsv
[ -s "$MAN" ] || { echo "missing $MAN" >&2; exit 1; }

# node ages of the observed and control runs
printf 'run\tnode\tmean\tet_lo\tet_hi\thpd_lo\thpd_hi\thpd_width\n' \
  > modelsel/ages_controls.tsv
for d in full_chain1 ds_top75 ds_top50 ds_top25 rand??_??? prior_only; do
  [ -s "$d/out.txt" ] || continue
  awk -v d="$d" -f ages.awk "$d/out.txt"
done >> modelsel/ages_controls.tsv
echo "rows in modelsel/ages_controls.tsv: $(( $(wc -l < modelsel/ages_controls.tsv) - 1 ))"

# per fraction: control mean, SD, max, k (controls >= observed), p; warn if controls are missing
awk -F'\t' 'NR>1 && $2=="t_n15" {print $1"\t"$3*100}' modelsel/ages_controls.tsv |
awk -F'\t' -v man="$MAN" '
  BEGIN {                                   # controls expected per fraction
    while ((getline line < man) > 0) {
      split(line, c, "\t")
      if (c[1] ~ /^rand/) { split(c[1], a, "_"); expct[substr(a[1], 5)]++ }
    }
  }
  $1 ~ /^ds_top/ { obs[substr($1,7)] = $2; next }
  $1 ~ /^rand/   { split($1,a,"_"); f=substr(a[1],5)
                   n[f]++; ctl[f,n[f]]=$2; next }
                 { other[$1]=$2 }
  END {
    printf "\nfull dataset: %.2f Ma\n\n", other["full_chain1"]
    printf "%-8s %7s %6s %9s %8s %8s %8s %8s %8s\n",
           "subset","obs","n","ctl_mean","ctl_sd","ctl_max","k>=obs","perm_p","dS_eff"
    split("75 50 25", fr, " ")
    short = 0
    for (j = 1; j <= 3; j++) {
      f = fr[j]
      if (!(f in obs)) { printf "top%-5s no observed run\n", f; continue }
      if (!n[f])       { printf "top%-5s no controls read\n", f; short = 1; continue }
      m = 0; for (i=1;i<=n[f];i++) m += ctl[f,i]; m /= n[f]
      v = 0; for (i=1;i<=n[f];i++) v += (ctl[f,i]-m)^2
      sd = (n[f] > 1) ? sqrt(v/(n[f]-1)) : 0
      k = 0; mx = 0
      for (i=1;i<=n[f];i++) { if (ctl[f,i] >= obs[f]) k++
                              if (ctl[f,i] > mx) mx = ctl[f,i] }
      printf "top%-5s %7.2f %6d %9.2f %8.2f %8.2f %8d %8.4f %8.2f\n",
             f, obs[f], n[f], m, sd, mx, k, (1+k)/(n[f]+1), obs[f]-m
      if (n[f] < expct[f]) short = 1
    }
    print "\ndS_eff = observed minus control mean, i.e. the elevation not"
    print "explained by dataset size."
    if (short) {
      print "\nWARNING: fewer controls read than the manifest lists:"
      for (j = 1; j <= 3; j++) { f = fr[j]
        if (n[f] < expct[f]) printf "  top%s: %d of %d\n", f, n[f], expct[f] }
      print "Resubmit 10.controls.sbatch (finished controls are skipped), then rerun this."
    }
  }'
