from KiraanWaktuSolat import Takwim

# Wujudkan sebuah kelas bagi Takwim()
Penang = Takwim()

#Tetapkan lokasi. Jika tidak ditetapkan, maka lokasi Pusat Falak Sheikh Tahir akan digunakan
Penang.latitude = 5.411
Penang.longitude = 100.2

#Cari azimut
azimut = Penang.azimut_kiblat()


#print(azimut)