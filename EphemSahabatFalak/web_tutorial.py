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
        latitude = 5.41144
        longitude = 100.19672
        elevation = 40.0
        now = datetime.now()
        year = now.year
        month = now.month
        day = now.day
        hour = now.hour
        minute = now.minute
        second = now.second
        original_date = now.strftime("%d-%m-%YT%H:%M:%S")
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

        try:

            if request.form['datetime'] != "":
                original_date = request.form['datetime']
                shortened_date = datetime.strptime(original_date,
                                                   "%Y-%m-%dT%H:%M:%S").strftime("%d-%m-%Y")
                parsed_date = datetime.strptime(original_date,
                                                "%Y-%m-%dT%H:%M:%S")
                second = parsed_date.second
            else:
                parsed_date = original_date
        except Exception as err:
            parsed_date = datetime.strptime(original_date,
                                            "%Y-%m-%dT%H:%M")
            shortened_date = datetime.strptime(original_date,
                                               "%Y-%m-%dT%H:%M").strftime("%d-%m-%Y")
            second = 0
            print(f"Error: {err}. Original Date: {original_date}")
        finally:

            year = parsed_date.year
            month = parsed_date.month
            day = parsed_date.day
            hour = parsed_date.hour
            minute = parsed_date.minute

        if request.form['timezone'] == "lain":
            zone = None
            timezone = 'lain'
        else:
            zone = 'Asia/Kuala_Lumpur'
            timezone = 'asia_KL'

        pemerhati = Takwim(
            latitude=latitude, longitude=longitude, elevation=elevation,
            temperature=temperature, pressure=pressure, year=year,
            month=month, day=day, hour=hour, minute=minute, second=second,
            zone=zone
        )
        if request.form["pilihan"] == "waktu-solat":
            pilihan = "waktu-solat"
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
                temperature=temperature, pressure=pemerhati.pressure,
                timezone=timezone,
                original_date=original_date, pilihan=pilihan)
        elif request.form["pilihan"] == "DataBulanMatahari":
            pilihan = "DataBulanMatahari"
            azimut_b = pemerhati.moon_azimuth(angle_format='string')
            altitud_b = pemerhati.moon_altitude(angle_format='string')
            fasa_bulan = pemerhati.moon_phase(angle_format='degree')
            if fasa_bulan <= 80:
                fasa_bulan = 'Bulan Sabit Muda'
            elif fasa_bulan > 80 and fasa_bulan <= 100:
                fasa_bulan = 'Bulan Separa Muda'
            elif fasa_bulan > 100 and fasa_bulan <= 170:
                fasa_bulan = 'Bulan Hampir Purnama Muda'
            elif fasa_bulan > 170 and fasa_bulan <= 190:
                fasa_bulan = 'Bulan Purnama'
            elif fasa_bulan > 190 and fasa_bulan <= 260:
                fasa_bulan = 'Bulan Hampir Purnama Tua'
            elif fasa_bulan > 260 and fasa_bulan <= 280:
                fasa_bulan = 'Bulan Separa Tua'
            else:
                fasa_bulan = 'Bulan Sabit Tua'

            jarak_lengkung = pemerhati.elongation_moon_sun(angle_format='string')
            azimut_m = pemerhati.sun_azimuth(angle_format='string')
            altitud_m = pemerhati.sun_altitude(angle_format='string')

            return render_template(
                "index.html", pilihan=pilihan, azimut_b=azimut_b, azimut_m=azimut_m,
                altitud_b=altitud_b, altitud_m=altitud_m, fasa_bulan=fasa_bulan,
                parsed_date=pemerhati.current_time('string'), jarak_lengkung=jarak_lengkung,
                original_date=original_date
            )

        elif request.form["pilihan"] == "Kiblat":
            pilihan = "Kiblat"
            azimut_degree = float(pemerhati.azimut_kiblat())
            azimut = "{:.2f}".format(azimut_degree)
            jarak_float = float(pemerhati.jarak_kaabah())
            jarak = "{:.2f}".format(jarak_float)
            waktu_bayang = pemerhati.bayang_searah_kiblat(time_format="string")
            return render_template(
                "index.html", azimut=azimut, jarak=jarak, pilihan=pilihan,
                waktu_bayang=waktu_bayang, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pemerhati.pressure,
                original_date=original_date, timezone=timezone
            )

        elif request.form["pilihan"] == "efemerisKiblat":
            pilihan = "efemerisKiblat"
            azimut_degree = float(pemerhati.azimut_kiblat())
            azimut = "{:.2f}".format(azimut_degree)
            efemeris_kiblat = pemerhati.efemeris_kiblat(directory="web")
            return render_template(
                "index.html", azimut=azimut, efemeris_kiblat=efemeris_kiblat,
                latitud=latitude, longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pemerhati.pressure,
                original_date=original_date, timezone=timezone, pilihan=pilihan)

        elif request.form["pilihan"] == "waktuSolatBulanan":
            pilihan = "waktuSolatBulanan"
            takwim_bulanan = pemerhati.takwim_solat_bulanan(
                saat="ya", directory="web"
            )
            return render_template(
                "index.html", takwim_bulanan=takwim_bulanan, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pemerhati.pressure,
                original_date=original_date, timezone=timezone, pilihan=pilihan)

        elif request.form["pilihan"] == "TukarKalendar":
            pilihan = "TukarKalendar"
            tukar_tarikh = pemerhati.tukar_ke_tarikh_hijri()
            return render_template(
                "index.html", tukar_tarikh=tukar_tarikh, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pemerhati.pressure,
                original_date=original_date, timezone=timezone, pilihan=pilihan,
                parsed_date=parsed_date, shortened_date=shortened_date)

    else:
        return render_template("index.html")


@app.route("/sahabatfalakpro", methods=["POST", "GET"])
def sahabatfalakpro():
    if request.method == "POST":
        year = 2024
        lokasi_dic = {}
        if request.form["tahun"] != "":
            year = int(request.form["tahun"])
        if request.form['timezone'] == "lain":
            zone = None
            timezone = 'lain'
        else:
            zone = 'Asia/Kuala_Lumpur'
            timezone = 'asia_KL'

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
                    analisis_perbandingan_julat_isyak=multipoint[14],
                    timezone=timezone
                    )

        pemerhati = Takwim(year=year, zone=zone)
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
            analisis_perbandingan_julat_isyak=multipoint[14],
            timezone=timezone
            )
    return render_template("sahabatfalakpro.html")


@app.route("/sahabatfalakplus", methods=["POST", "GET"])
def sahabatfalakplus():
    if request.method == "POST":
        latitude = 5.41144
        longitude = 100.19672
        elevation = 40.0
        now = datetime.now()
        year = now.year
        month = now.month
        day = now.day
        hour = now.hour
        minute = now.minute
        second = now.second
        temperature = 27.0
        pressure = None
        original_date = now.strftime("%d-%m-%YT%H:%M:%S")

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

        try:

            if request.form['datetime'] != "":
                original_date = request.form['datetime']
                parsed_date = datetime.strptime(original_date,
                                                "%Y-%m-%dT%H:%M:%S")
                second = parsed_date.second
            else:
                parsed_date = original_date
        except Exception as err:
            parsed_date = datetime.strptime(original_date,
                                            "%Y-%m-%dT%H:%M")
            second = 0
            print(f"Error: {err}. Original Date: {original_date}")
        finally:
            year = parsed_date.year
            month = parsed_date.month
            day = parsed_date.day
            hour = parsed_date.hour
            minute = parsed_date.minute

        if request.form['timezone'] == "lain":
            zone = None
            timezone = 'lain'
        else:
            zone = 'Asia/Kuala_Lumpur'
            timezone = 'asia_KL'

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
        maghrib_2 = maghrib_0(time_format='datetime')
        maghrib = maghrib_0(time_format='string')
        pilihan = request.form["pilihan"]
        if request.form["pilihan"] == "dataHilalRingkas":
            moon_set = pemerhati.moon_set(time_format='string')
            ijtimak = pemerhati.moon_conjunction(time_format='string', topo=topo)
            moon_age = pemerhati.moon_age(t=maghrib_1, topo=topo)
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
            moon_distance_ratio = pemerhati.moon_distance(t=maghrib_1, topo=topo, compare=True)
            moon_distance_ratio = str(format(moon_distance_ratio, '.4f'))
            illumination = str(format(illumination, '.4f')) + "%"
            lag_time = pemerhati.lag_time()
            moon_rise = pemerhati.moon_rise(time_format='string')
            sun_rise = pemerhati.waktu_syuruk(time_format='string')
            crescent_width = pemerhati.lunar_crescent_width(
                t=maghrib_1, topo=topo, angle_format='degree'
            )
            crescent_width = str(format(crescent_width*60, '.3f')) + \
                "' (arka minit)"
            if 'Moon' in lag_time:
                odeh, yallop = 'Mudah Kelihatan', 'Mudah Kelihatan'
            else:
                odeh = pemerhati.Odeh_criteria(value='description', maghrib_pre_calculated=maghrib_2)
                yallop = pemerhati.Yallop_criteria(value='description', maghrib_pre_calculated=maghrib_2)
            mabims2021 = pemerhati.Mabims_2021_criteria(value='description', maghrib_pre_calculated=maghrib_1)
            malaysia2013 = pemerhati.Malaysia_2013_criteria(
                value='description', maghrib_pre_calculated=maghrib_1)
            muhammadiyah = pemerhati.Muhammadiyah_wujudul_hilal_criteria(
                value='description', maghrib_pre_calculated=maghrib_1
            )
            istanbul78 = pemerhati.Istanbul_1978_criteria(value='description', maghrib_pre_calculated=maghrib_1)

            if "days" not in moon_age:
                gambar = pemerhati.gambar_hilal_mabims(save=False, topo=topo)
                pngImage = io.BytesIO()
                # gambar.savefig(pngImage, format='png')
                FigureCanvas(gambar).print_png(pngImage)
                pngImage.seek(0)
                img_base64 = base64.b64encode(pngImage.getvalue()).decode('utf-8')
            else:
                img_base64 = None

            return render_template(
                "sahabatfalakplus.html", maghrib=maghrib, moon_set=moon_set,
                moon_age=moon_age, elongation=elongation, altitude=altitude,
                azimuth_difference=azimuth_difference, ijtimak=ijtimak,
                illumination=illumination, moon_distance_ratio=moon_distance_ratio,
                lag_time=lag_time, crescent_width=crescent_width,
                parsed_date=pemerhati.current_time('string'),
                odeh=odeh, yallop=yallop,
                mabims2021=mabims2021, malaysia2013=malaysia2013,
                muhammadiyah=muhammadiyah, istanbul78=istanbul78, topo=topo,
                original_date=original_date, pilihan=pilihan,
                gambar_Hilal=img_base64, latitud=latitude,
                longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pemerhati.pressure,
                timezone=timezone, moon_rise=moon_rise, sun_rise=sun_rise)

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
                parsed_date=pemerhati.current_time('string'),
                original_date=original_date, pilihan=pilihan,
                elong_min=elong_min, elong_max=elong_max, conj=conj, opp=opp,
                sudut_min=sudut_min, sudut_max=sudut_max,
                latitud=latitude, longitud=longitude, elevation=elevation,
                temperature=temperature, pressure=pemerhati.pressure,
                timezone=timezone)

        elif request.form["pilihan"] == "dataHilalPenuh":
            moon_age = pemerhati.moon_age(t=maghrib_1, topo=topo)
            ijtimak = pemerhati.moon_conjunction(time_format='string', topo=topo)
            if "days" not in moon_age:
                efemeris_hilal = pemerhati.efemeris_hilal(topo, 'web')
                moon_set = pemerhati.moon_set(time_format='string')
                return render_template(
                    "sahabatfalakplus.html", efemeris_hilal=efemeris_hilal,
                    topo=topo, parsed_date=pemerhati.current_time('string'),
                    pilihan=pilihan,
                    latitud=latitude, longitud=longitude, elevation=elevation,
                    temperature=temperature, pressure=pemerhati.pressure,
                    original_date=original_date, maghrib=maghrib,
                    moon_set=moon_set, timezone=timezone)
            return render_template(
                    "sahabatfalakplus.html", efemeris_hilal=None, ijtimak=ijtimak,
                    topo=topo, parsed_date=pemerhati.current_time('string'),
                    pilihan=pilihan,
                    latitud=latitude, longitud=longitude, elevation=elevation,
                    temperature=temperature, pressure=pemerhati.pressure,
                    original_date=original_date, maghrib=maghrib, timezone=timezone)
    return render_template(
            "sahabatfalakplus.html")


@app.route("/panduan")
def panduan():
    return render_template("panduan.html")


@app.route("/tentangkami")
def tentangkami():
    return render_template("tentangkami.html")


if __name__ == "__main__":
    app.run(host='0.0.0.0', port=8000, debug=True)
