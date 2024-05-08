from flask import (Flask, render_template, request)
from kiraan_waktu_solat import Takwim
from negeri import (Pulau_Pinang)
from datetime import datetime


app = Flask(__name__)


@app.route("/", methods=["POST", "GET"])
def home():
    if request.method == "POST":
        print("Method is post")
        pemerhati = Takwim()
        latitude = 5.41144
        longitude = 100.19672
        elevation = 40.0
        year = datetime.now().year
        month = datetime.now().month
        day = datetime.now().day
        hour = datetime.now().hour
        minute = datetime.now().minute
        second = datetime.now().second
        temperature = 27.0
        pressure = None

        if request.form["latitud"] != "":
            latitude = float(request.form["latitud"])

        if request.form["longitud"] != "":
            longitude = float(request.form["longitud"])

        if request.form["elevation"] != "":
            elevation = float(request.form["elevation"])

        if request.form["temperature"] != "":
            temperature = float(request.form["temperature"])

        if request.form["pressure"] != "":
            pressure = float(request.form["pressure"])

        if request.form['datetime'] != "":
            parsed_date = datetime.strptime(request.form['datetime'],
                                            "%Y-%m-%dT%H:%M:%S")
            year = parsed_date.year
            month = parsed_date.month
            day = parsed_date.day
            hour = parsed_date.hour
            minute = parsed_date.minute
            second = parsed_date.second

        pemerhati = Takwim(
            latitude=latitude, longitude=longitude, elevation=elevation,
            temperature=temperature, pressure=pressure, year=year,
            month=month, day=day, hour=hour, minute=minute, second=second
        )
        if request.form["pilihan"] == "waktu-solat":
            subuh = pemerhati.waktu_subuh(time_format='string')
            syuruk = pemerhati.waktu_syuruk(time_format='string')
            zohor = pemerhati.waktu_zohor(time_format='string')
            asar = pemerhati.waktu_asar(time_format='string')
            maghrib = pemerhati.waktu_maghrib(time_format='string')
            isyak = pemerhati.waktu_isyak(time_format='string')
            return render_template(
                "index.html", subuh=subuh, syuruk=syuruk, zohor=zohor,
                asar=asar, maghrib=maghrib, isyak=isyak, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure)
        elif request.form["pilihan"] == "Kiblat":
            azimut_degree = float(pemerhati.azimut_kiblat())
            azimut = "{:.2f}".format(azimut_degree)
            jarak_float = float(pemerhati.jarak_kaabah())
            jarak = "{:.2f}".format(jarak_float)
            waktu_bayang = pemerhati.bayang_searah_kiblat(time_format="string")
            return render_template(
                "index.html", azimut=azimut, jarak=jarak,
                waktu_bayang=waktu_bayang, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure
            )

        elif request.form["pilihan"] == "efemerisKiblat":
            efemeris_kiblat = pemerhati.efemeris_kiblat(directory="web")
            return render_template("index.html",
                                   efemeris_kiblat=efemeris_kiblat)

        elif request.form["pilihan"] == "waktuSolatBulanan":
            takwim_bulanan = pemerhati.takwim_solat_bulanan(
                saat="ya", directory="web"
            )
            return render_template(
                "index.html", takwim_bulanan=takwim_bulanan)

    else:
        return render_template("index.html")


@app.route("/sahabatfalakpro", methods=["POST", "GET"])
def sahabatfalakpro():
    if request.method == "POST":
        year = 2024
        lokasi_dic = {}
        if request.form["tahun"] != "":
            year = int(request.form["tahun"])

        if request.form["negeri"] != "Tiada":
            if request.form["negeri"] == "Penang":
                pemerhati = Pulau_Pinang()
                pemerhati.year = year
                multipoint = pemerhati.takwim_solat_tahunan_multipoint_Penang(
                    directory="web"
                )
                lokasi_dic = {
                    "A": (5.282574, 100.177892, 40),
                    "B": (5.284480, 100.17720, 40),
                    "C": (5.142472, 100.502861, 40),
                    "D": (5.551423, 100.528866, 40),
                    "E": (5.472839, 100.180822, 40)}
                return render_template(
                    "sahabatfalakpro.html",
                    lokasi_dic=lokasi_dic,
                    takwim_multipoint=multipoint[0],
                    lokasi_pilihan=multipoint[1],
                    perbandingan_julat=multipoint[2],
                    analisis_subuh=multipoint[3],
                    analisis_syuruk=multipoint[4],
                    analisis_zohor=multipoint[5],
                    analisis_asar=multipoint[6],
                    analisis_maghrib=multipoint[7],
                    analisis_isyak=multipoint[8],
                    analisis_perbandingan_julat_subuh=multipoint[9],
                    analisis_perbandingan_julat_syuruk=multipoint[10],
                    analisis_perbandingan_julat_zohor=multipoint[11],
                    analisis_perbandingan_julat_asar=multipoint[12],
                    analisis_perbandingan_julat_maghrib=multipoint[13],
                    analisis_perbandingan_julat_isyak=multipoint[14]
                    )

        pemerhati = Takwim(year=year)
        if request.form["latitude0"] != "":
            for key, value in request.form.items():
                if key.startswith('titik'):
                    current_location = value
                    lokasi_dic[current_location] = ()
                elif key.startswith(('latitude', 'longitude', 'elevation')):
                    suffix = key.split(' ')[0]
                    lokasi_dic[current_location] += (
                        float(value) if suffix != "elevation" else int(value),
                    )

            multipoint = pemerhati.takwim_solat_tahunan_multipoint(
                directory="web", **lokasi_dic
            )
        return render_template(
            "sahabatfalakpro.html",
            lokasi_dic=lokasi_dic,
            takwim_multipoint=multipoint[0],
            lokasi_pilihan=multipoint[1],
            perbandingan_julat=multipoint[2],
            analisis_subuh=multipoint[3],
            analisis_syuruk=multipoint[4],
            analisis_zohor=multipoint[5],
            analisis_asar=multipoint[6],
            analisis_maghrib=multipoint[7],
            analisis_isyak=multipoint[8],
            analisis_perbandingan_julat_subuh=multipoint[9],
            analisis_perbandingan_julat_syuruk=multipoint[10],
            analisis_perbandingan_julat_zohor=multipoint[11],
            analisis_perbandingan_julat_asar=multipoint[12],
            analisis_perbandingan_julat_maghrib=multipoint[13],
            analisis_perbandingan_julat_isyak=multipoint[14]
            )
    return render_template("sahabatfalakpro.html")


@app.route("/tentangkami")
def tentangkami():
    return render_template("tentangkami.html")


if __name__ == "__main__":
    app.run(debug=True)
