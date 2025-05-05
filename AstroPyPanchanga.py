from datetime import datetime
import pytz
from timezonefinder import TimezoneFinder
from astropy.time import Time
from astropy.coordinates import Angle, get_body, EarthLocation, SkyCoord
from astropy.table import Table
import numpy as np
import astropy.units as u

from plot_panchanga import *

# The difference in the position of the vernal equinox between the sāyana and nirāyana rāśis
# Lahiri ayanāṃśa at 2000 as coordinates being used by Astropy is J2000
# Refernce https://www.astro.com/swisseph/slae/2000/slae_2000.pdf -
# use ayanāṃśa for January 2000 and subtract delta t = 63.83 to convert to terrestrial time (J2000)
#
ayanāṃśa = Angle("23d50m21.17s").degree

tithi_extent = 360 / 30

nakṣatra_extent = 360 / 27
pāda_extent = nakṣatra_extent / 4

rāśi_extent = 30


def calc_tithi(sun_lambda, moon_lambda, tithi_names_file="tithi_names.tex"):
    """
    Paramters
    ---------
    sun_lambda: float
        The value of geocentrictrueecliptic longitude of the sun in degree

    moon_lambda: float
        The value of geocentrictrueecliptic longitude of the moon in degree

    tithi_names_file -- str, default = "tithi_names.tex"
        File with names of tithis

    Returns
    ---------
    final_tithi: float
        Real value of tithis elapsed since śuklapakṣa prathama represnted as 0, amāvāsya = 29

    final_tithi_name: str
        Name of tithi
    """

    if sun_lambda > moon_lambda:
        tithi_theta = 360 - (sun_lambda - moon_lambda)
    else:
        tithi_theta = moon_lambda - sun_lambda

    final_tithi = tithi_theta / tithi_extent

    tithi_names_tab = Table.read(tithi_names_file, format="latex")
    tithi_names = tithi_names_tab["tithi_names"].data

    final_tithi_name = tithi_names[int(final_tithi)]

    return final_tithi, final_tithi_name


def cal_vāra(date, vāra_names_file="vaara_names.tex"):
    """
    Paramters
    ---------
    date: datetime.datetime
        The date of interest

    vāra_names_file -- str, default = "tithi_names.tex"
        File with names of vāras

    Returns
    ---------
    vāra: str
        The vāra corresponding to the input date
    """

    day_of_the_week_at_t = date.strftime("%A")

    vāra_names_tab = Table.read(vāra_names_file, format="latex")
    weekday = vāra_names_tab["weekday"].data
    vāra_names = vāra_names_tab["vaara"].data

    vāra = vāra_names[np.where(weekday == day_of_the_week_at_t)][0]

    return vāra


def cal_nakṣatra(moon_lambda, nakṣatra_names_file="nakshatra_names.tex"):
    """
    Paramters
    ---------
    moon_lambda: float
        The value of geocentrictrueecliptic longitude of the moon in degree

    nakṣatra_names_file -- str, default = "nakshatra_names.tex"
        File with names of vāras

    Returns
    ---------
    nakṣatra: str
        The vāra corresponding to the input date
    """
    nakṣatra_names_tab = Table.read(nakṣatra_names_file, format="latex")
    nakṣatra_names = nakṣatra_names_tab["names"].data

    nakṣatra_id = np.floor((moon_lambda - ayanāṃśa) / nakṣatra_extent).astype(int)

    nakṣatra = nakṣatra_names[nakṣatra_id]

    # pada

    pāda_remainder = (moon_lambda - ayanāṃśa) % nakṣatra_extent

    final_pāda = int(pāda_remainder / pāda_extent) + 1

    pāda_values = ["१", "२", "३", "४"]

    final_pāda_text = pāda_values[final_pāda - 1]

    return nakṣatra, final_pāda_text


def calc_yoga(sun_lambda, moon_lambda, yoga_names_file="yoga_names.tex"):
    """
    Paramters
    ---------
    sun_lambda: float
        The value of geocentrictrueecliptic longitude of the sun in degree

    moon_lambda: float
        The value of geocentrictrueecliptic longitude of the moon in degree

    yoga_names_file -- str, default = "yoga_names.tex"
        File with names of yogas

    Returns
    ---------
    yoga : str
        The vaue of yoga at given time
    """
    ### -------- Yoga --------###
    yoga_names_tab = Table.read(yoga_names_file, format="latex")
    yoga_names = yoga_names_tab["names"].data

    yoga_extent = nakṣatra_extent

    yoga_id = np.floor(
        ((sun_lambda + moon_lambda - 2 * ayanāṃśa) % 360) / yoga_extent
    ).astype(int)

    yoga = yoga_names[yoga_id]

    return yoga


def calc_karaṇa(sun_lambda, moon_lambda, karaṇa_names_file="karana_names.tex"):
    """
    Paramters
    ---------
    sun_lambda: float
        The value of geocentrictrueecliptic longitude of the sun in degree

    moon_lambda: float
        The value of geocentrictrueecliptic longitude of the moon in degree

    karaṇa_names_file -- str, default = "karana_names.tex"
        File with names of tithis

    Returns
    ---------
    """

    if sun_lambda > moon_lambda:
        karaṇa_theta = 360 - (sun_lambda - moon_lambda)
    else:
        karaṇa_theta = moon_lambda - sun_lambda

    karaṇa_id = int(karaṇa_theta / tithi_extent * 2)

    karaṇa_names_tab = Table.read(karaṇa_names_file, format="latex")
    karaṇa_names = karaṇa_names_tab["names"].data

    karaṇa = karaṇa_names[karaṇa_id]

    return karaṇa


def calc_pañcāṅga(
    location,
    time,
    time_format="%Y-%m-%d %H:%M:%S",
    plotfilename="panchanga_at_test_time.pdf",
):
    """
    Calculates nakshatra and tithi at input time and makes plot of grahas

    Paramters
    ---------
    location = enter address of your location. Eg: "Mumbai, India"
    time = time at which to calculate panchanga
    time_format = format of input time
    filename = name of file to write the plot of position of grahas in nakshatra and rashi

    Returns
    ---------
    Plot saved at filename
    """

    observing_location = EarthLocation.of_address(location)
    tz_obj = TimezoneFinder()
    timezone_location = tz_obj.timezone_at(
        lng=observing_location.lon.value, lat=observing_location.lat.value
    )

    fmt = time_format
    date_str = time
    tz = pytz.timezone(timezone_location)
    test_date = datetime.strptime(date_str, fmt)

    local_time = tz.localize(test_date, is_dst=None)
    test_date_utc = local_time.astimezone(pytz.utc)

    test_date_utc_time = Time(test_date_utc.strftime(fmt), format="iso", scale="utc")

    print("Time in UTC : ", test_date_utc_time)

    moon_coord = get_body("moon", test_date_utc_time, location=observing_location)
    sun_coord = get_body("sun", test_date_utc_time, location=observing_location)

    moon_lambda = moon_coord.geocentrictrueecliptic.lon.value
    sun_lambda = sun_coord.geocentrictrueecliptic.lon.value

    ### -------- tithi -------- ###

    tithi, tithi_name = calc_tithi(sun_lambda, moon_lambda)

    ### -------- vāra -------- ###

    vāra = cal_vāra(test_date)

    ### --------- nakṣatra ----------- ###

    nakṣatra, pāda = cal_nakṣatra(moon_lambda)

    ### --------- yoga ------------ ###

    yoga = calc_yoga(sun_lambda, moon_lambda)

    ### --------- karaṇa ------------- ###

    karaṇa = calc_karaṇa(sun_lambda, moon_lambda)

    ### --------- lagna ------------- ###

    rising_rāśi = SkyCoord(
        alt=0 * u.deg,
        az=90 * u.deg,
        frame="altaz",
        obstime=test_date_utc_time,
        location=observing_location,
    )

    ### --------- Plotting ------------ ###

    make_circle_plot(
        test_date_utc_time,
        location,
        date_str,
        nakṣatra_extent,
        rāśi_extent,
        ayanāṃśa,
        tithi,
        tithi_name,
        vāra,
        nakṣatra,
        pāda,
        yoga,
        karaṇa,
        moon_lambda,
        sun_lambda,
        rising_rāśi.geocentrictrueecliptic.lon.value,
        language="Devanagari",
        plotfile=plotfilename,
    )

    make_jatakam_plot(
        test_date_utc_time,
        location,
        date_str,
        rāśi_extent,
        ayanāṃśa,
        tithi_name,
        vāra,
        nakṣatra,
        pāda,
        yoga,
        karaṇa,
        moon_lambda,
        sun_lambda,
        rising_rāśi.geocentrictrueecliptic.lon.value,
        language="Devanagari",
    )

    return tithi_name, vāra, nakṣatra + pāda, karaṇa


location = "Kanchipuram, India"
date_str = "2026-04-30 06:00:00"

print(calc_pañcāṅga(location, date_str))
