from datetime import datetime
import pytz
from timezonefinder import TimezoneFinder
from astropy.time import Time
from astropy.coordinates import Angle, get_body, EarthLocation
from astropy.table import Table
import numpy as np

# The difference in the position of the vernal equinox between the sāyana and nirāyana rāśis
# Lahiri ayanāṃśa at 2000 as coordinates being used by Astropy is J2000
ayanāṃśa = (
    Angle("23d50m30s").degree
)  # I have guessed this value works for current panchanga. May have to verify this


tithi_extent = 360 / 30
nakṣatra_extent = 360 / 27


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

    return nakṣatra


def calc_pañcāṅga(
    location,
    time,
    time_format="%Y-%m-%d %H:%M:%S",
    filename="nakshatra_at_test_time.pdf",
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

    ### -------- Thithi -------- ###

    tithi, tithi_name = calc_tithi(sun_lambda, moon_lambda)

    ### -------- Vaara -------- ###

    vāra = cal_vāra(test_date)

    ### --------- Nakṣatra ----------- ###

    nakṣatra = cal_nakṣatra(moon_lambda)

    return tithi_name, vāra, nakṣatra


location = "Bengaluru, India"
date_str = "2025-04-15 00:13:00"

print(calc_pañcāṅga(location, date_str, filename="nakshatra_at_test_time.pdf"))
