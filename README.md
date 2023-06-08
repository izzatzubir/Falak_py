# Falak_py

<p>EphemSahabatFalak ialah sebuah file python yang mengandungi sebuah python Class iaitu Takwim() yang mempunyai pelbagai method. </p>
<p> Program ini menggunakan data ephemeris dari JPL Horizon (NASA) iaitu de440s. Pengguna boleh untuk mengubah kepada ephemeris lain </p>
<p> Program ini menggunapakai pakej Skyfield </p>
<p> Program ini dihasilkan oleh Izzat Zubir, dengan sokongan dari keluarga, beserta maklum balas dari Dr. Abdul Halim Abdul Aziz, Panel Pakar Falak JAKIM </p>


<h1> Parameter </h1>
<p> Terdapat 13 parameter yang boleh diubah:</p>
<p> nilai dalam kurungan ialah nilai 'default' </p>
<p>1. latitud (5.41144) </p>
<p>2. longitud (100.19672) --> latitud dan longitud Pusat Falak Sheikh Tahir</p>
<p>3. elevation (40) --> dalam unit meter</p>
<p>4. year (tahun semasa)</p>
<p>5. month (bulan semasa)</p>
<p>6. day (hari semasa)</p>
<p>7. hour (jam semasa)</p>
<p>8. minute (minit semasa)</p>
<p>9. second (saat semasa) --> tetapan semasa diambil menggunakan method datetime.now()</p>
<p>10. zone ('Asia/Kuala_Lumpur')</p>
<p>11. temperature (27) --> unit darjah Celcius</p>
<p>12. pressure (1010) --> unit millibar</p>
<p>13. ephem ('de440s.bsp') --> ephemeris de440s adalah ephemeris terkini dari JPL Horizon. Rujuk: https://rhodesmill.org/skyfield/planets.html#choosing-an-ephemeris </p>
<hr>
<h2>Method dalam Takwim</h2>

<h2> Waktu Solat </h2>
<h3>waktu_subuh(time_format = 'default', altitude = 'default')</h3>
<p>Waktu subuh ialah ketika altitud matahari -18.0</p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu subuh dalam python datetime </p>
<p> time_format = 'string' -> format waktu subuh dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu subuh dalam skyfield.Time </p>

<h4>altitude</h4>
<p> altitude boleh dipilih antara -24 hingga -12. Nilai default ialah -18 </p>
<hr>

<h3>waktu_syuruk(time_format = 'default', altitude = 'default')</h3>
<p>Waktu syuruk ialah ketika altitud matahari -0.8333 di ufuk timur</p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu dalam python datetime </p>
<p> time_format = 'string' -> format waktu dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu dalam skyfield.Time </p>

<h4>altitude</h4>
<p> altitude boleh dipilih antara 0 hingga -4. Nilai default ialah -0.8333 </p>
<hr>

<h3>waktu_zawal(time_format = 'default')</h3>
<p>Waktu zawal ialah ketika pusat matahari melintasi meridian</p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu dalam python datetime </p>
<p> time_format = 'string' -> format waktu dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu dalam skyfield.Time </p>
<hr>
<h3>waktu_zohor(time_format = 'default')</h3>
<p>Waktu zohor ialah waktu zawal + 1 minit 5 saat</p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu dalam python datetime </p>
<p> time_format = 'string' -> format waktu dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu dalam skyfield.Time </p>
<hr>
<h3>waktu_asar(time_format = 'default')</h3>
<p>Waktu asar ialah ketika bayang matahari bersamaan bayang sesebuah objek ketika waktu zawal ditambah panjang bayang tersebut</p>
<p>Rumus bagi altitud matahari asar ialah acot(cot(altitud matahari ketika zawal)+1)<p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu dalam python datetime </p>
<p> time_format = 'string' -> format waktu dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu dalam skyfield.Time </p>
<hr>
<h3>waktu_maghrib(time_format = 'default', altitude = 'default')</h3>
<p>Waktu maghrib ialah ketika altitud matahari -0.8333 di ufuk barat</p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu dalam python datetime </p>
<p> time_format = 'string' -> format waktu dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu dalam skyfield.Time </p>

<h4>altitude</h4>
<p> altitude boleh dipilih antara 0 hingga -4. Nilai default ialah -0.8333 </p>
<hr>
<h3>waktu_isyak(time_format = 'default', altitude = 'default')</h3>
<p>Waktu isyak ialah ketika altitud matahari -18.0</p>
<h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu subuh dalam python datetime </p>
<p> time_format = 'string' -> format waktu subuh dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu subuh dalam skyfield.Time </p>

<h4>altitude</h4>
<p> altitude boleh dipilih antara -24 hingga -12. Nilai default ialah -18 </p>
<hr>
<h3>bayang_searah_kiblat(time_format = 'default')</h3>
<p> Metod ini memberikan dua waktu, iaitu waktu mula dan waktu akhir ketika azimut matahari berada setentang dengan arah kiblat </p>
<p> Selisih yang digunapakai ialah 0.3 darjah atau 18 arka minit, ke kiri dan ke kanan azimut kiblat sebenar (36 arka minit)</p>
 <h3>parameter</h3>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu subuh dalam python datetime </p>
<p> time_format = 'string' -> format waktu subuh dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu subuh dalam skyfield.Time </p>
<hr>

## Efemeris matahari dan bulan

<h3>moon_altitude(t = current_time, angle_format = 'skylib', temperature = None, pressure = None)</h3>
  <p> Metod ini memberikan altitud bulan pada waktu yang ditetapkan </p>
  <p> Suhu dan Tekanan adalah berdasarkan tetapan yang dibuat </p>
 <h4>angle_format</h4>
 <p> angle_format = 'skylib' -> format sudut dalam skylib.Angle.radians
 <p> angle_format = 'degree' -> format sudut dalam skylib.Angle.degrees
 <p> angle_format = 'string' -> format sudut dalam string, darjah° minit' saat"

<hr>
<h3>moon_azimuth(t = None, angle_format = 'skylib')</h3>
<p> Metod ini memberikan azimut bulan pada waktu yang ditetapkan </p>
<p> Azimut bulan tidak mengambil kira pembiasan atmosfera (secara tradisi) </p>
 <h4>angle_format</h4>
 <p> angle_format = 'skylib' -> format sudut dalam skylib.Angle.radians
 <p> angle_format = 'degree' -> format sudut dalam skylib.Angle.degrees
 <p> angle_format = 'string' -> format sudut dalam string, darjah° minit' saat"
<hr>
<h3>moon_distance(t = None, angle_format = 'skylib')</h3>
  <p> Metod ini memberikan jarak bulan dalam unit au </p>
  <hr>
<h3>sun_altitude(t = current_time, angle_format = 'skylib', temperature = None, pressure = None)</h3>
  <p> Metod ini memberikan altitud matahari pada waktu yang ditetapkan </p>
  <p> Suhu dan Tekanan adalah berdasarkan tetapan yang dibuat </p>
 <h4>angle_format</h4>
 <p> angle_format = 'skylib' -> format sudut dalam skylib.Angle.radians
 <p> angle_format = 'degree' -> format sudut dalam skylib.Angle.degrees
 <p> angle_format = 'string' -> format sudut dalam string, darjah° minit' saat"

<hr>
<h3>sun_azimuth(t = None, angle_format = 'skylib')</h3>
<p> Metod ini memberikan azimut matahari pada waktu yang ditetapkan </p>
<p> Azimut matahari tidak mengambil kira pembiasan atmosfera (secara tradisi) </p>
 <h4>angle_format</h4>
 <p> angle_format = 'skylib' -> format sudut dalam skylib.Angle.radians
 <p> angle_format = 'degree' -> format sudut dalam skylib.Angle.degrees
 <p> angle_format = 'string' -> format sudut dalam string, darjah° minit' saat"
<hr>
<h3>sun_distance(t = None, angle_format = 'skylib')</h3>
  <p> Metod ini memberikan jarak sun dalam unit au </p>
  <hr>
<h3>elongation_moon_sun(t = None, topo = 'topo', angle_format = 'skylib')</h3>
<p> Metod ini memberikan elongasi antara bulan dan matahari, berdasarkan tetapan yang dibuat </p>
 <h4>topo</h4>
 <p> topo = 'geo' atau 'geocentric' -> nilai elongasi berdasarkan elongasi geocentric. </p>
 <p> topo - 'topo' atau selainnya -> nilai elongasi berdasarkan elongasi topocentric </p>
  
<h3>moon_set(time_format = 'default')</h3>
<p> Metod ini memberikan waktu terbenam bulan, iaitu ketika altitud bulan di bawah -0.8333 darjah </p>
<p> Jika tiada waktu terbenam, maka ia akan memberikan "Moon does not set on tarikh semasa"</p>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu subuh dalam python datetime </p>
<p> time_format = 'string' -> format waktu subuh dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu subuh dalam skyfield.Time </p>
<hr>
 <h3>moon_rise(time_format = 'default')</h3>
<p> Metod ini memberikan waktu terbit bulan, iaitu ketika altitud bulan di melebihi -0.8333 darjah </p>
<p> Jika tiada waktu terbit, maka ia akan memberikan "Moon does not rise on tarikh semasa"</p>
<h4>time_format</h4>
<p> time_format = 'datetime' -> format waktu subuh dalam python datetime </p>
<p> time_format = 'string' -> format waktu subuh dalam string (hh:mm:ss) </p>
<p> time_format = 'default' -> format waktu subuh dalam skyfield.Time </p>
<hr>

<h3>azimut_kiblat()</h3>
  <p> Metod ini memberikan azimut kiblat bagi kawasan yang ditetapkan (berdasarkan latitud dan longitud) </p>
  <p> unit azimut kiblat ialah dalam skyfield.Angle.degrees </p>
  <p> latitud dan longitud kaabah yang ditetapkan ialah 21.422487 Utara dan 39.826206 Timur. </p>
  <p> jejari Bumi yang digunakan ialah 6371000 meter </p>
  <p> formula yang digunakan ialah berdasarkan https://www.omnicalculator.com/other/azimuth</p>
  <hr>
 <h3>jarak_kaabah()</h3>
  <p> Metod ini memberikan jarak ke kaabah bagi kawasan yang ditetapkan (berdasarkan latitud dan longitud) </p>
  <p> unit jarak adalah dalam meter (integer) </p>
  <p> latitud dan longitud kaabah yang ditetapkan ialah 21.422487 Utara dan 39.826206 Timur. </p>
  <p> jejari Bumi yang digunakan ialah 6371000 meter </p>
  <p> formula yang digunakan ialah berdasarkan https://www.omnicalculator.com/other/azimuth iaitu haversine formula</p>
  <hr>

  <h2>Jadual</h2>
  <h3>efemeris_hilal(topo = 'topo')</h3>
  <p> Metod ini memberikan sebuah jadual efemeris dalam bentuk pandas.Dataframe </p>
  <p> Jadual ini mengandungi waktu, elongasi bulan-matahari, altitud bulan, azimut bulan, altitud matahari dan azimut matahari</p>
  <p> Efemeris ini bermula satu jam sebelum waktu maghrib dan berakhir sejam selepas waktu maghrib </p>
  <h4>topo</h4>
  <p> tetapan topo boleh dirujuk pada method elongasi_moon_sun()</p>
  <hr>
  <h3>takwim_solat_bulanan(altitud_subuh ='default', altitud_syuruk ='default', altitud_maghrib ='default', altitud_isyak ='default', saat = 'tidak')</h3>
  <p> Metod ini memberikan sebuah takwim solat dalam bentuk pandas.Dataframe</p>
  <p> Jadual ini mengandungi tarikh, mula bayang searah kiblat, tamat bayang searah kiblat, subuh, syuruk, zohor, asar, maghrib, isyak</p>
  <p> Jadual ini bermula pada 1 haribulan bagi bulan dan tahun yang ditetapkan, dan berakhir pada hari akhir bulan </p>
  <h4>tetapan altitud dan saat</h4>
  <p> Tetapan altitud subuh hingga isyak boleh dirujuk pada method waktu solat masing-masing </p>
  <p> Tetapan saat ialah bagi memilih sama ada mahu bundarkan kepada minit bawah (syuruk dan bayang tamat) atau minit atas (selainnya).</p>
  <p> Tetapan asal bagi saat ialah bundarkan, namun jika mahu nilai saat, ubah kepada saat = 'ya' </p>
  
  
  
