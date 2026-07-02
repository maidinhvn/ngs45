#!/bin/bash
# Resolve the smallest suitable HiFi + Illumina-paired-WGS run per species.
# Writes bench/manifest.tsv:  order  species  hifi_srr  hifi_MB  illu_srr  illu_MB
export PATH="/path/to/.conda/envs/edirect_env/bin:$PATH"
OUT=/path/to/ngs245/bench/manifest.tsv
: > "$OUT"

species=(
 "Poales|Oryza sativa"
 "Zingiberales|Musa acuminata"
 "Solanales|Solanum lycopersicum"
 "Asterales|Helianthus annuus"
 "Lamiales|Sesamum indicum"
 "Caryophyllales|Beta vulgaris"
 "Vitales|Vitis vinifera"
 "Fabales|Glycine max"
 "Malvales|Gossypium hirsutum"
 "Sapindales|Citrus sinensis"
)

# Only scan the first ~60 runinfo records (fast even for species with 10^5 runs);
# we just need *a* small suitable run, not the global minimum. 100s cap per query.
pick_hifi() {  # $1 species -> "SRR<TAB>MB"
  timeout 100 bash -c '
    esearch -db sra -query "\"'"$1"'\"[Organism] AND \"PacBio SMRT\"[Platform] AND (HiFi[All Fields] OR CCS[All Fields])" </dev/null 2>/dev/null \
     | efetch -format runinfo -stop 60 2>/dev/null' \
   | awk -F, 'NR>1 && $1~/^[SED]RR/ && ($20 ~ /Sequel|Revio/) && $8+0>=150 && $8+0<=12000 {print $8"\t"$1}' \
   | sort -n | head -1 | awk -F'\t' '{print $2"\t"$1}'
}
pick_illu() {  # $1 species -> "SRR<TAB>MB"
  timeout 100 bash -c '
    esearch -db sra -query "\"'"$1"'\"[Organism] AND ILLUMINA[Platform] AND WGS[Strategy] AND biomol_dna[Properties] AND paired[Layout]" </dev/null 2>/dev/null \
     | efetch -format runinfo -stop 60 2>/dev/null' \
   | awk -F, 'NR>1 && $1~/^[SED]RR/ && $16=="PAIRED" && $8+0>=150 && $8+0<=900 {print $8"\t"$1}' \
   | sort -n | head -1 | awk -F'\t' '{print $2"\t"$1}'
}

for e in "${species[@]}"; do
  order="${e%%|*}"; sp="${e##*|}"
  h=$(pick_hifi "$sp"); i=$(pick_illu "$sp")
  hs="${h%%$'\t'*}"; hm="${h##*$'\t'}"; is="${i%%$'\t'*}"; im="${i##*$'\t'}"
  printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$order" "$sp" "${hs:-NA}" "${hm:-NA}" "${is:-NA}" "${im:-NA}" | tee -a "$OUT"
done
echo "DONE_SELECT" >> "$OUT"
