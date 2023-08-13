# Collect unique filenames
file_names=$(ls *_GO.txt | sed 's/_GO\.txt//')
file_names_header=$(ls *_GO.txt | sed 's/_GO\.txt//' | sed 's/.faa.tsv//' | tr "\n" "\t")

# Create header line
echo -e "Pattern\t$file_names_header" > occurence_table1.tsv

# Process each pattern and count occurrences in each file
for pattern in $(grep -Eo 'GO:[0-9]+' --no-filename *_GO.txt | sort -u); do
    row="$pattern"
    for file in $file_names; do
        count=$(grep -Eo "$pattern" "${file}_GO.txt" | wc -l)
        row="$row\t$count"
    done
    echo -e "$row" >> occurence_table1.tsv
done
