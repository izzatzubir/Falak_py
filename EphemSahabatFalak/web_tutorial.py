from flask import (Flask, render_template, request)
from kiraan_waktu_solat import Takwim
from negeri import (Pulau_Pinang)
from datetime import datetime
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
import io
import os
import base64
import pandas as pd
from dotenv import load_dotenv
from pathlib import Path
THIS_FOLDER = Path(__file__).parent.resolve()

load_dotenv()


app = Flask(__name__)


@app.route("/", methods=["POST", "GET"])
def home():
    if request.method == "POST":
        mapbox_access_token = os.getenv('MAPBOX_ACCESS_TOKEN')
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
            azimut_b = pemerhati.moon_azimuth(angle_format='skylib')
            azimut_degree = azimut_b.degrees
            azimut_string = azimut_b.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')

            altitud_b = pemerhati.moon_altitude(angle_format='skylib')
            altitud_degree = altitud_b.degrees
            altitud_string = altitud_b.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')

            azimut_b = f"{azimut_string} ({round(azimut_degree,4)})"
            altitud_b = f"{altitud_string} ({round(altitud_degree,4)})"

            fasa_bulan = pemerhati.moon_phase(angle_format='degree')
            fasa_bulan_asing = fasa_bulan
            fasa_bulan_asing = "{:.2f}".format(fasa_bulan_asing)
            if fasa_bulan <= 75:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Sabit Muda)"
            elif fasa_bulan > 75 and fasa_bulan <= 105:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Separa Muda)"
            elif fasa_bulan > 105 and fasa_bulan <= 165:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Hampir Purnama Muda)"
            elif fasa_bulan > 165 and fasa_bulan <= 195:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Purnama)"
            elif fasa_bulan > 195 and fasa_bulan <= 255:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Hampir Purnama Tua)"
            elif fasa_bulan > 255 and fasa_bulan <= 285:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Separa Tua)"
            else:
                fasa_bulan = f"{fasa_bulan_asing} (Bulan Sabit Tua)"

            jarak_lengkung = pemerhati.elongation_moon_sun(angle_format='skylib')
            jarak_lengkung_degree = jarak_lengkung.degrees
            jarak_lengkung_string = jarak_lengkung.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')
            jarak_lengkung = f"{jarak_lengkung_string} ({round(jarak_lengkung_degree,4)})"

            azimut_m = pemerhati.sun_azimuth(angle_format='skylib')
            azimut_degree = azimut_m.degrees
            azimut_string = azimut_m.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')
            altitud_m = pemerhati.sun_altitude(angle_format='skylib')
            altitud_degree = altitud_m.degrees
            altitud_string = altitud_m.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')

            azimut_m = f"{azimut_string} ({round(azimut_degree,4)})"
            altitud_m = f"{altitud_string} ({round(altitud_degree,4)})"

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
                original_date=original_date, timezone=timezone, mapbox_access_token=mapbox_access_token
            )

        elif request.form["pilihan"] == "efemerisKiblat":
            pilihan = "efemerisKiblat"
            pemerhati.second = 0
            pemerhati.time = pemerhati.current_time()
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
    running_remote = os.getenv("RUNNING_REMOTE")
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
                    timezone=timezone, running_remote=running_remote
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
            timezone=timezone, running_remote=running_remote
            )
    lokasi_dic = {
                    "A": (5.282574, 100.177892, 40),
                    "B": (5.284480, 100.17720, 40),
                    "C": (5.142472, 100.502861, 40),
                    "D": (5.551423, 100.528866, 40),
                    "E": (5.472839, 100.180822, 40)}
    takwim_multipoint = pd.read_csv(
                    THIS_FOLDER / "takwim_tahunan_penang_2024.csv", index_col='date')
    lokasi_pilihan = pd.read_csv(THIS_FOLDER / 'lokasi_pilihan_penang_2024.csv')
    perbandingan_julat = pd.read_csv(THIS_FOLDER / 'perbandingan_julat_penang_2024.csv')
    analisis_subuh = pd.read_csv(THIS_FOLDER / 'analisis_subuh_penang_2024.csv')
    analisis_syuruk = pd.read_csv(THIS_FOLDER / 'analisis_syuruk_penang_2024.csv')
    analisis_zohor = pd.read_csv(THIS_FOLDER / 'analisis_zohor_penang_2024.csv')
    analisis_asar = pd.read_csv(THIS_FOLDER / 'analisis_asar_penang_2024.csv')
    analisis_maghrib = pd.read_csv(THIS_FOLDER / 'analisis_maghrib_penang_2024.csv')
    analisis_isyak = pd.read_csv(THIS_FOLDER / 'analisis_isyak_penang_2024.csv')
    analisis_perbandingan_julat_subuh = pd.read_csv(
        THIS_FOLDER / 'analisis_perbandingan_julat_subuh_penang_2024.csv')
    analisis_perbandingan_julat_syuruk = pd.read_csv(
        THIS_FOLDER / 'analisis_perbandingan_julat_syuruk_penang_2024.csv')
    analisis_perbandingan_julat_zohor = pd.read_csv(
        THIS_FOLDER / 'analisis_perbandingan_julat_zohor_penang_2024.csv')
    analisis_perbandingan_julat_asar = pd.read_csv(
        THIS_FOLDER / 'analisis_perbandingan_julat_asar_penang_2024.csv')
    analisis_perbandingan_julat_maghrib = pd.read_csv(
        THIS_FOLDER / 'analisis_perbandingan_julat_maghrib_penang_2024.csv')
    analisis_perbandingan_julat_isyak = pd.read_csv(
        THIS_FOLDER / 'analisis_perbandingan_julat_isyak_penang_2024.csv')

    return render_template(
        "sahabatfalakpro.html", running_remote=running_remote, lokasi_dic=lokasi_dic,
        takwim_multipoint=takwim_multipoint,
        lokasi_pilihan=lokasi_pilihan,
        perbandingan_julat=perbandingan_julat,
        analisis_subuh=analisis_subuh,
        analisis_syuruk=analisis_syuruk,
        analisis_zohor=analisis_zohor,
        analisis_asar=analisis_asar,
        analisis_maghrib=analisis_maghrib,
        analisis_isyak=analisis_isyak,
        analisis_perbandingan_julat_subuh=analisis_perbandingan_julat_subuh,
        analisis_perbandingan_julat_syuruk=analisis_perbandingan_julat_syuruk,
        analisis_perbandingan_julat_zohor=analisis_perbandingan_julat_zohor,
        analisis_perbandingan_julat_asar=analisis_perbandingan_julat_asar,
        analisis_perbandingan_julat_maghrib=analisis_perbandingan_julat_maghrib,
        analisis_perbandingan_julat_isyak=analisis_perbandingan_julat_isyak)


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
                t=maghrib_1, angle_format='skylib', topo=topo)
            elongation_string = elongation.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')
            elongation = f"{elongation_string}° ({round(elongation.degrees, 4)}°)"

            azimuth_difference = pemerhati.daz(
                t=maghrib_1, angle_format='skylib')
            azimuth_difference_string = azimuth_difference.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')
            azimuth_difference = f"{azimuth_difference_string}° ({round(azimuth_difference.degrees,4)})°"

            altitude = pemerhati.moon_altitude(
                t=maghrib_1, angle_format='skylib', topo=topo)
            altitude_string = altitude.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')
            altitude = f"{altitude_string}° ({round(altitude.degrees,4)})°"

            illumination = pemerhati.moon_illumination(
                t=maghrib_1, topo=topo)
            moon_distance_ratio = pemerhati.moon_distance(t=maghrib_1, topo=topo, compare=True)
            moon_distance_ratio = f"{round(moon_distance_ratio,4)}"
            illumination = f"{round(illumination,4)}%"
            lag_time = pemerhati.lag_time()
            moon_rise = pemerhati.moon_rise(time_format='string')
            sun_rise = pemerhati.waktu_syuruk(time_format='string')
            crescent_width = pemerhati.lunar_crescent_width(
                t=maghrib_1, topo=topo, angle_format='skylib'
            )
            crescent_width_string = crescent_width.dstr(format=u'{0}{1}° {2:02}′ {3:02}.{4:0{5}}″')
            crescent_width = f"{crescent_width_string} ({round(crescent_width.degrees*60, 3)}') (arka minit)"
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
                efemeris_hilal_present = True
                return render_template(
                    "sahabatfalakplus.html", efemeris_hilal=efemeris_hilal,
                    efemeris_hilal_present=efemeris_hilal_present,
                    topo=topo, parsed_date=pemerhati.current_time('string'),
                    pilihan=pilihan,
                    latitud=latitude, longitud=longitude, elevation=elevation,
                    temperature=temperature, pressure=pemerhati.pressure,
                    original_date=original_date, maghrib=maghrib,
                    moon_set=moon_set, timezone=timezone)
            efemeris_hilal_present = False
            return render_template(
                    "sahabatfalakplus.html", efemeris_hilal=None, ijtimak=ijtimak,
                    efemeris_hilal_present=efemeris_hilal_present,
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
