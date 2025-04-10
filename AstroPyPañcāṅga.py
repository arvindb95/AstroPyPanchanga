from datetime import datetime
import pytz
from timezonefinder import TimezoneFinder
from astropy.time import Time
from astropy.coordinates import Angle, get_body, EarthLocation
from astropy.table import Table

# The difference in the position of the vernal equinox between the sāyana and nirāyana rāśis
# Ref : Indian Astronomy : An Introduction - S Balachandra Rao - Pg 35

ayanāṃśa = Angle("23d47m14.1s").degree


def calc_pañcāṅga(
    location,
    time,
    time_format="%Y-%m-%d %H:%M:%S",
    filename="nakshatra_at_test_time.pdf",
):
    """
    Calculates nakshatra and tithi at input time and makes plot of grahas

    Inputs :
    location = enter address of your location. Eg: "Mumbai, India"
    time = time at which to calculate panchanga
    time_format = format of input time
    filename = name of file to write the plot of position of grahas in nakshatra and rashi

    Output :
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

    moon_coord = get_body("moon", test_date_utc_time)
    sun_coord = get_body("sun", test_date_utc_time)

    moon_lambda = moon_coord.geocentrictrueecliptic.lon.value
    sun_lambda = sun_coord.geocentrictrueecliptic.lon.value

    ### -------- Thithi -------- ###

    each_tithi = 360 / 30

    if sun_lambda > moon_lambda:
        tithi_theta = 360 - (sun_lambda - moon_lambda)
    else:
        tithi_theta = moon_lambda - sun_lambda

    final_tithi = tithi_theta / each_tithi
    
    tithi_names_tab = Table.read("tithi_names.tex", format="latex")
    tithi_names = tithi_names_tab["tithi_names"].data
    
    final_tithi_name = tithi_names[int(final_tithi)]

    return final_tithi_name 

location = "Bengaluru, India"
date_str = "2025-03-30 19:55:00"

print(calc_pañcāṅga(location, date_str, filename="nakshatra_at_test_time.pdf"))


