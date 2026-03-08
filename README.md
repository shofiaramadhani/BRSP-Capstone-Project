**Topik: Analisis Transkriptomik Jalur Interferon Tipe I pada Psoriasis**

Tools: GEO, GEO2R, R (Rstudio)

Proyek ini mengeksplorasi disregulasi genetik pada pasien Psoriasis menggunakan dataset GSE14905, berfokus pada peran Interferon Tipe I dan kaskade inflamasi.

**Catatan:**

Healthy = Kontrol sehat;
Disease = Pasien Psoriasis

**Langkah-Langkah Analisis (GEO2R & RStudio)**
1. Identifikasi & Pengelompokan: Akses dataset GSE14905 di GEO dan bagi sampel menjadi dua kategori: 21 sampel "Healthy" dan 61 sampel "Disease".
2. Ekstraksi R Script: Gunakan fitur R script pada GEO2R untuk mendapatkan kerangka kode analisis diferensial ekspresi gen (DEG).
3. Pemrosesan di RStudio: Jalankan script tersebut di RStudio untuk melakukan pengunduhan data, normalisasi dengan paket limma, serta transformasi data.
4. Analisis Statistik: Hitung nilai log2 fold change dan adjusted p-value (< 0,05) menggunakan metode Benjamini-Hochberg untuk mengidentifikasi gen signifikan.
5. Visualisasi Data: Hasilkan plot distribusi gen seperti MD Plot, Volcano Plot, Venn Diagram, dan Heatmap untuk melihat perbedaan profil ekspresi gen yang kontras.

**Hasil**

Diagram Mean Difference Plot
<img width="707" height="423" alt="MD plot" src="https://github.com/user-attachments/assets/4ba6463c-0c78-4bdf-b55c-988467da29c7" />

Diagram Volcano Plot
<img width="732" height="424" alt="Volcano plot 2" src="https://github.com/user-attachments/assets/4b4e1632-6411-415c-9f66-2c5552ebbc57" />

Diagram Venn
<img width="707" height="423" alt="Venn diagram" src="https://github.com/user-attachments/assets/8d9f777a-3880-475c-b80f-7ac78b3cdd5a" />

1. Total Gen Signifikan: Ditemukan sebanyak 13.619 gen yang mengalami ekspresi diferensial signifikan antara kondisi sehat dan psoriasis.
2. Pola Ekspresi: Terjadi peningkatan drastis gen yang diinduksi interferon seperti STAT1 dan ISG15 pada jaringan lesi psoriasis.

**Analisis Enrichment**
1. Gene Ontology (GO): Menunjukkan dominasi proses "Response to virus" yang mengindikasikan aktivasi jalur interferon tipe I sebagai respons imun utama.
<img width="800" height="1200" alt="GO 3" src="https://github.com/user-attachments/assets/ad46d704-d269-4629-bfb1-6607b8c5dbb0" />

2. KEGG Pathway: Mengidentifikasi aktivasi jalur IL-17 signaling, TNF signaling, serta DNA Replication yang mendasari inflamasi kronis dan hiperproliferasi sel kulit.
<img width="800" height="500" alt="KEGG 2" src="https://github.com/user-attachments/assets/5f448b4c-3521-4099-a00a-59599ee2ad02" />
