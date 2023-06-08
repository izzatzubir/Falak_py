from KiraanWaktuSolat import Takwim

## Wujudkan sebuah kelas bagi Takwim()
Penang = Takwim()

#Tetapkan lokasi. Jika tidak ditetapkan, maka lokasi Pusat Falak Sheikh Tahir akan digunakan
Penang.latitude = 5.411
Penang.longitude = 100.2

#Tetapkan hari cerapan, bulan dan tahun
Penang.year = 2023
Penang.month = 6
Penang.day = 18

#Hasilkan takwim
takwim = Penang.efemeris_hilal()

#pindahkan ke dalam excel. pastikan format xlsx
nama_file = 'Efemeris_Hilal_ZulHijjah_1444.xlsx'

#takwim.to_excel(nama_file) 