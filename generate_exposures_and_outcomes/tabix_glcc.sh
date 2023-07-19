# use tabix

while read p; do
filename=$(basename "$p" .gz)


tabix  "$p" 11:6137923-6751280 | awk -v fname="$filename"  '{print $0 "\t" fname "\t" region}' >> glcc_output
tabix  "$p" 2:113875548-113891591 | awk -v fname="$filename"  '{print $0 "\t" fname "\t" region}' >> glcc_output
tabix  "$p" 1:154077819-154741926 | awk -v fname="$filename"  '{print $0 "\t" fname "\t" region}' >> glcc_output

rm *.tbi
done < glcc_files