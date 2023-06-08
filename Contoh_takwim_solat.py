from KiraanWaktuSolat import Takwim

# Wujudkan sebuah kelas bagi Takwim()
Penang = Takwim()

#Tetapkan lokasi. Jika tidak ditetapkan, maka lokasi Pusat Falak Sheikh Tahir akan digunakan
Penang.latitude = 5.411
Penang.longitude = 100.2

#Tetapkan bulan dan tahun
Penang.year = 2023
Penang.month = 7

#Hasilkan takwim
takwim = Penang.takwim_solat_bulanan()

#pindahkan ke dalam excel. pastikan format xlsx
nama_file = 'Takwim_Solat_Ogos_2023.xlsx' 

#takwim.to_excel(nama_file)
#print(takwim)
