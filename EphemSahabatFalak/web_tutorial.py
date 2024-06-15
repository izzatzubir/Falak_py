from flask import (Flask, render_template, request)
from kiraan_waktu_solat import Takwim
from negeri import (Pulau_Pinang)
from datetime import datetime
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import io
import base64


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
            original_date = request.form['datetime']
            try:
                parsed_date = datetime.strptime(request.form['datetime'],
                                                "%Y-%m-%dT%H:%M:%S")
                second = parsed_date.second
            except Exception:
                parsed_date = datetime.strptime(request.form['datetime'],
                                                "%Y-%m-%dT%H:%M")
                second = 0

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
                temperature=temperature, pressure=pressure,
                original_date=original_date)
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
                temperature=temperature, pressure=pressure,
                original_date=original_date
            )

        elif request.form["pilihan"] == "efemerisKiblat":
            efemeris_kiblat = pemerhati.efemeris_kiblat(directory="web")
            return render_template(
                "index.html", efemeris_kiblat=efemeris_kiblat,
                latitud=latitude, longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure,
                original_date=original_date)

        elif request.form["pilihan"] == "waktuSolatBulanan":
            takwim_bulanan = pemerhati.takwim_solat_bulanan(
                saat="ya", directory="web"
            )
            return render_template(
                "index.html", takwim_bulanan=takwim_bulanan, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure,
                original_date=original_date)

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
                    directory="web")
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


@app.route("/sahabatfalakplus", methods=["POST", "GET"])
def sahabatfalakplus():
    if request.method == "POST":
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
            original_date = request.form['datetime']
            try:
                parsed_date = datetime.strptime(request.form['datetime'],
                                                "%Y-%m-%dT%H:%M:%S")
                second = parsed_date.second
            except Exception:
                parsed_date = datetime.strptime(request.form['datetime'],
                                                "%Y-%m-%dT%H:%M")
                second = 0

        if request.form['timezone'] == "lain":
            zone = None
            timezone = 'lain'
        else:
            zone = 'Asia/Kuala_Lumpur'
            timezone = 'asia_KL'

            year = parsed_date.year
            month = parsed_date.month
            day = parsed_date.day
            hour = parsed_date.hour
            minute = parsed_date.minute

        pemerhati = Takwim(
            latitude=latitude, longitude=longitude, elevation=elevation,
            temperature=temperature, pressure=pressure, year=year,
            month=month, day=day, hour=hour, minute=minute, second=second,
            zone=zone
        )
        # Tambah pilihan zon waktu jika selain Asia/Kuala_Lumpur

        topo = request.form["topo_pilihan"]
        maghrib_0 = pemerhati.waktu_maghrib
        maghrib_1 = maghrib_0()
        maghrib = maghrib_0(time_format='string')
        pilihan = request.form["pilihan"]
        if request.form["pilihan"] == "dataHilalRingkas":
            moon_set = pemerhati.moon_set(time_format='string')
            ijtimak = pemerhati.moon_conjunction(time_format='string')
            moon_age = pemerhati.moon_age(topo=topo)
            elongation = pemerhati.elongation_moon_sun(
                t=maghrib_1, angle_format='degree', topo=topo)
            elongation = str(format(elongation, '.4f')) + "°"
            azimuth_difference = pemerhati.daz(
                t=maghrib_1, angle_format='degree')
            azimuth_difference = str(format(azimuth_difference, '.4f')) + "°"
            altitude = pemerhati.moon_altitude(
                t=maghrib_1, angle_format='degree', topo=topo)
            altitude = str(format(altitude, '.4f')) + "°"
            illumination = pemerhati.moon_illumination(
                t=maghrib_1, topo=topo)
            illumination = str(format(illumination, '.4f')) + "%"
            lag_time = pemerhati.lag_time()
            crescent_width = pemerhati.lunar_crescent_width(
                t=maghrib_1, topo=topo, angle_format='degree'
            )
            crescent_width = str(format(crescent_width*60, '.3f')) + \
                "' (arka minit)"
            if 'Moon' in lag_time:
                odeh, yallop = 'Mudah Kelihatan', 'Mudah Kelihatan'
            else:
                odeh = pemerhati.Odeh_criteria(value='description')
                yallop = pemerhati.Yallop_criteria(value='description')
            mabims2021 = pemerhati.Mabims_2021_criteria(value='description')
            malaysia2013 = pemerhati.Malaysia_2013_criteria(
                value='description')
            muhammadiyah = pemerhati.Muhammadiyah_wujudul_hilal_criteria(
                value='description'
            )
            istanbul78 = pemerhati.Istanbul_1978_criteria(value='description')
            gambar = pemerhati.gambar_hilal_mabims(save=False, topo=topo)
            pngImage = io.BytesIO()
            # gambar.savefig(pngImage, format='png')
            FigureCanvas(gambar).print_png(pngImage)
            pngImage.seek(0)
            img_base64 = base64.b64encode(pngImage.getvalue()).decode('utf-8')

            return render_template(
                "sahabatfalakplus.html", maghrib=maghrib, moon_set=moon_set,
                moon_age=moon_age, elongation=elongation, altitude=altitude,
                azimuth_difference=azimuth_difference, ijtimak=ijtimak,
                illumination=illumination,
                lag_time=lag_time, crescent_width=crescent_width,
                parsed_date=str(parsed_date)[:11], odeh=odeh, yallop=yallop,
                mabims2021=mabims2021, malaysia2013=malaysia2013,
                muhammadiyah=muhammadiyah, istanbul78=istanbul78, topo=topo,
                original_date=original_date, pilihan=pilihan,
                gambar_Hilal=img_base64, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure, timezone=timezone)

        elif request.form["pilihan"] == "PerbandinganElongasi":
            fasa = pemerhati.moon_phase(topo=topo, angle_format='degree')
            fasa = str(format(fasa, '.4f')) + "°"
            elong = pemerhati.elongation_moon_sun(topo=topo,
                                                  angle_format='degree')
            elong = str(format(elong, '.4f')) + "°"
            elong_min, sudut_min = pemerhati.moon_min_elongation(
                'string', topo)
            sudut_min = str(format(sudut_min, '.4f')) + "°"
            elong_max, sudut_max = pemerhati.moon_max_elongation(
                'string', topo)
            sudut_max = str(format(sudut_max, '.4f')) + "°"
            conj = pemerhati.moon_conjunction('string', topo)
            opp = pemerhati.moon_opposition('string', topo)

            return render_template(
                "sahabatfalakplus.html", topo=topo, fasa=fasa, elong=elong,
                parsed_date=str(parsed_date)[:19], original_date=original_date,
                pilihan=pilihan, elong_min=elong_min, elong_max=elong_max,
                conj=conj, opp=opp, sudut_min=sudut_min, sudut_max=sudut_max,
                latitud=latitude, longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure, timezone=timezone)

        elif request.form["pilihan"] == "dataHilalPenuh":
            efemeris_hilal = pemerhati.efemeris_hilal(topo, 'Yes')
            moon_set = pemerhati.moon_set(time_format='string')
            return render_template(
                "sahabatfalakplus.html", efemeris_hilal=efemeris_hilal,
                topo=topo, parsed_date=str(parsed_date)[:19], pilihan=pilihan,
                latitud=latitude, longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pressure,
                original_date=original_date, maghrib=maghrib,
                moon_set=moon_set, timezone=timezone)
    return render_template(
            "sahabatfalakplus.html")


@app.route("/tentangkami")
def tentangkami():
    return render_template("tentangkami.html")


if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8000, debug=True)
