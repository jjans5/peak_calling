ls -lah fragment_files/enterocytes/bonobo_enterocytes.fragments.tsv.gz

mkdir subsampled_fragments/enterocytes

cat -n 100 fragment_files/enterocytes/bonobo_enterocytes.fragments.tsv.gz > subsampled_fragments/enterocytes/bonobo_enterocytes.fragments.tsv.gz