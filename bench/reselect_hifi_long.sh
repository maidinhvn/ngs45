#!/bin/bash
# Re-select a LONG-READ HiFi run per species (avgLength >= MINLEN) so reads can
# span a full 45S unit -> required by easy45. Picks the smallest such run for a
# fast download. Writes manifest_hifi_long.tsv: order species srr avgLen size_MB
export PATH="/path/to/.conda/envs/edirect_env/bin:$PATH"
OUT=/path/to/ngs245/bench/manifest_hifi_long.tsv
MINLEN=7000
: > "$OUT"

species=(
 "Poales|Oryza sativa"                 "Zingiberales|Musa acuminata"
 "Solanales|Solanum lycopersicum"      "Asterales|Helianthus annuus"
 "Lamiales|Sesamum indicum"            "Caryophyllales|Beta vulgaris"
 "Vitales|Vitis vinifera"              "Fabales|Glycine max"
 "Malvales|Gossypium hirsutum"         "Sapindales|Citrus sinensis"
)

pick_long_hifi() {  # $1 species -> "SRR  avgLen  MB"  (smallest size with avgLen>=MINLEN)
  timeout 110 bash -c '
    esearch -db sra -query "\"'"$1"'\"[Organism] AND \"PacBio SMRT\"[Platform] AND (HiFi[All Fields] OR CCS[All Fields])" </dev/null 2>/dev/null \
     | efetch -format runinfo -stop 80 2>/dev/null' \
   | awk -F, -v m='"$MINLEN"' 'NR>1 && $1~/^[SED]RR/ && ($20 ~ /Sequel|Revio/) && $7+0>=m && $8+0<=20000 {print $8"\t"$1"\t"$7}' \
   | sort -n | head -1 | awk -F'\t' '{print $2"\t"$3"\t"$1}'
}

for e in "${species[@]}"; do
  order="${e%%|*}"; sp="${e##*|}"
  r=$(pick_long_hifi "$sp")
  srr="${r%%$'\t'*}"; rest="${r#*$'\t'}"; avg="${rest%%$'\t'*}"; mb="${rest##*$'\t'}"
  printf '%s\t%s\t%s\t%s\t%s\n' "$order" "$sp" "${srr:-NA}" "${avg:-NA}" "${mb:-NA}" | tee -a "$OUT"
done
echo "DONE_RESELECT" >> "$OUT"
