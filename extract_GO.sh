for file in *.faa.tsv; do
    if [[ -f "$file" ]]; then
        output_file="${file%.faa}_GO.txt"
        grep -Eo 'GO:[0-9]+' "$file" > "$output_file"
    fi
done
