# Kiểm tra kết quả buổi sáng — ngs45 overnight benchmark

Job chạy nền (nice -19 + ionice idle → chỉ dùng CPU nhàn rỗi, nhường người khác).

## Cách kiểm tra nhanh

```bash
cd /path/to/ngs245/bench

# 1. Xong chưa?  (file này xuất hiện = đã hoàn tất toàn bộ)
ls OVERNIGHT_DONE 2>/dev/null && echo "HOÀN TẤT" || echo "còn đang chạy"

# 2. Tiến độ / log realtime
tail -40 overnight.log

# 3. Bảng kết quả (điền dần theo từng loài)
column -t -s$'\t' benchmark_results.tsv

# 4. Bản tóm tắt đẹp (ghi khi kết thúc)
cat SUMMARY.txt
```

## Bảng `benchmark_results.tsv` — ý nghĩa cột

| Cột | Nghĩa |
|---|---|
| ngs45_unit_len / its_len | độ dài đơn vị 45S & ITS từ **Illumina (ngs45)** |
| ngs45_ribo_sites | số vị trí ribotype dị hợp (tín hiệu lai/allopolyploid) |
| easy45_unit_len / its_len | từ **HiFi (easy45)** |
| **unitPID/mm** | %identity đơn vị 45S: **ngs45(Illumina) ↔ easy45(HiFi)** /số mismatch |
| ITS_ngsVSez | %id ITS: ngs45 ↔ easy45 |
| **ITS_ngsVSgb** | %id ITS: ngs45 ↔ **GenBank** (validation độc lập) |
| ITS_ezVSgb | %id ITS: easy45 ↔ GenBank |

10 loài đa dạng bộ: Poales, Zingiberales (2 monocot) + Solanales, Asterales, Lamiales,
Caryophyllales, Vitales, Fabales, Malvales, Sapindales (8 eudicot).
(Sesamum: không có HiFi → chỉ ngs45 + GenBank.)
Cộng Apiales đã làm trước: **Panax ginseng** (ngs45↔easy45 = 99.97%) & **Polyscias filicifolia** (100% vs ref).

## Nếu một loài lỗi
Không sao — driver bỏ qua và chạy tiếp. Xem log riêng:
`results/<Ten_loai>/ngs45.log`, `results/<Ten_loai>/easy45.log`, và `data/<Ten_loai>/*.log` (download).
Chạy lại chỉ cần: `bash run_overnight.sh` (tự bỏ qua bước đã xong — resumable).

## Dữ liệu
- Download: `data/<loài>/illu_1.fastq, illu_2.fastq, hifi.fastq`
- Kết quả: `results/<loài>/ngs45_out/, easy45_out/, its_ref.fasta`
