import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.transforms as transforms
import matplotlib as mpl
from astropy.coordinates import Angle, get_body, EarthLocation
from astropy.table import Table
from astropy.time import Time
import astropy.units as u
from mpl_toolkits.axes_grid1 import make_axes_locatable

#### Colors for grahas ####
# Sun - Surya - Red - circle
# Moon - Chandra - White - semicircle
# Mars - Mangala - Red - triangle
# Mercury - Budha - Green - bow (dhanusha)
# Jupiter - Guru - Yellow - lotus
# Venus - Shukra - White - square
# Saturn - Shani - Black - snake
# Rahu - Blue - crocodile
# Ketu - Grey - sword


def calc_rahu_ketu_pos(
    known_eclipse_time, test_date_utc_time, moon_lambda_known_eclipse
):
    ketu_speed = 360 / 18.612958  # deg per year

    if test_date_utc_time.jd > known_eclipse_time.jd:
        time_diff = ((test_date_utc_time - known_eclipse_time).value * u.d).to(u.year)
        ketu_degrees_moved = (ketu_speed * time_diff.value) % 360
        ketu_lambda = moon_lambda_known_eclipse - ketu_degrees_moved
        if ketu_lambda < 0:
            ketu_lambda = 360 - ketu_lambda
    else:
        time_diff = ((known_eclipse_time - test_date_utc_time).value * u.d).to(u.year)
        ketu_degrees_moved = (ketu_speed * time_diff.value) % 360
        ketu_lambda = moon_lambda_known_eclipse + ketu_degrees_moved
        if ketu_degrees_moved > 360:
            ketu_lambda = 360 - ketu_lambda

    if ketu_lambda < 180:
        rahu_lambda = 180 + ketu_lambda
    else:
        rahu_lambda = (180 + ketu_lambda) - 360

    return rahu_lambda, ketu_lambda


def plot_moon_phase(day, drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )
    if 0 < day <= 7.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="w",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - day * 12 * np.pi / 180),
            2 * radius,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 7.5 < day <= 15:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="w",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - day * 12 * np.pi / 180),
            2 * radius,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 15 < day <= 22.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="darkslategrey",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (day - 15) * 12 * np.pi / 180),
            2 * radius,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 22.5 < day <= 30:
        circle = mpatches.Circle(
            center,
            radius,
            ec="darkslategrey",
            fc="darkslategrey",
            lw=0.5,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (day - 15) * 12 * np.pi / 180),
            2 * radius,
            fc="darkslategrey",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc="w",
            ec=None,
            transform=trans,
            clip_on=False,
            zorder=10,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    return 0


def plot_sun(drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    circle = mpatches.Circle(
        center, radius, ec="orange", fc="yellow", transform=trans, clip_on=False
    )
    ax.add_patch(circle)

    list_of_patches = []
    for theta in np.arange(0, 2 * np.pi, np.pi / 6):
        triangle = mpatches.RegularPolygon(
            [1.2 * radius * np.cos(theta), 1.2 * radius * np.sin(theta)],
            3,
            radius=0.2 * radius,
            ec="yellow",
            fc="orange",
            orientation=theta + np.pi / 6,
            transform=trans,
            clip_on=False,
            linewidth=0.5,
        )
        ax.add_patch(triangle)

    return 0


def plot_inner_graha_phase(graha, angle_to_sun, drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    if graha == "budha":
        bright_side_color = "green"
        dark_side_color = "darkgreen"
    elif graha == "shukra":
        bright_side_color = "white"
        dark_side_color = "lightgrey"
    else:
        print(
            "That is not an inner graha! Use the outer graha function to plot ", graha
        )
        return

    if 0 < (angle_to_sun / 12) <= 7.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (angle_to_sun / 12) * 12 * np.pi / 180),
            2 * radius,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 7.5 < (angle_to_sun / 12) <= 15:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
        )
        ellipse = mpatches.Ellipse(
            center,
            2 * radius * np.sin(np.pi / 2 - (angle_to_sun / 12) * 12 * np.pi / 180),
            2 * radius,
            fc=bright_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            90,
            270,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 15 < (angle_to_sun / 12) <= 22.5:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
        )
        ellipse = mpatches.Ellipse(
            center,
            2
            * radius
            * np.sin(np.pi / 2 - ((angle_to_sun / 12) - 15) * 12 * np.pi / 180),
            2 * radius,
            fc=bright_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            270,
            90,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    if 22.5 < (angle_to_sun / 12) <= 30:
        circle = mpatches.Circle(
            center,
            radius,
            ec=dark_side_color,
            fc=bright_side_color,
            lw=0.5,
            transform=trans,
            clip_on=False,
        )
        ellipse = mpatches.Ellipse(
            center,
            2
            * radius
            * np.sin(np.pi / 2 - ((angle_to_sun / 12) - 15) * 12 * np.pi / 180),
            2 * radius,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        wedge = mpatches.Wedge(
            center,
            radius,
            270,
            90,
            fc=dark_side_color,
            ec=None,
            transform=trans,
            clip_on=False,
        )
        ax.add_patch(circle)
        ax.add_patch(wedge)
        ax.add_patch(ellipse)

    return 0


def plot_outer_graha(graha, drawing_origin, radius, fig, ax):
    center = np.array([0, 0])
    trans = fig.dpi_scale_trans + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    if graha == "mangala":
        mcolor = "red"
    elif graha == "guru":
        mcolor = "yellow"
    elif graha == "shani":
        mcolor = "black"
    else:
        print("That is not an outer graha! Use the inner graha function to plot", graha)
        return

    circle = mpatches.Circle(
        center, radius, ec="w", fc=mcolor, lw=0.5, transform=trans, clip_on=False
    )
    ax.add_patch(circle)

    return 0


def plot_rahu(drawing_origin, scale, fig, ax):
    center = np.array([0, 0])
    trans = transforms.Affine2D().scale(scale) + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    from matplotlib.path import Path

    verts = [
        (-1, 1),  # Vert1
        (1.0, 1.0),  # Vert2
        (1.0, -1.0),  # Pivot1
        (0, -1.0),  # Vert3
        (-1.0, -1.0),  # Pivot2
        (-1.0, 1.0),
        (-1.0, 1.0),
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.CURVE3,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    patch = mpatches.PathPatch(
        path,
        facecolor="blue",
        edgecolor="lightblue",
        lw=0.5 * scale / 5,
        transform=trans,
        clip_on=False,
    )
    ax.add_patch(patch)

    return 0


def plot_ketu(drawing_origin, scale, fig, ax):
    center = np.array([0, 0])
    trans = transforms.Affine2D().scale(scale) + transforms.ScaledTranslation(
        drawing_origin[0], drawing_origin[1], ax.transData
    )

    from matplotlib.path import Path

    verts = [
        (-1, 1),  # Vert1
        (1.0, 1.0),  # Vert2
        (1.0, -1.0),  # Vert3
        (-1.0, -1.0),  # Vert4
        (0, 0),  # Vert5
        (-1.0, 1.0),  # Vert6
        (-1.0, 1.0),
    ]

    codes = [
        Path.MOVETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.LINETO,
        Path.CLOSEPOLY,
    ]

    path = Path(verts, codes)

    patch = mpatches.PathPatch(
        path,
        facecolor="dimgrey",
        edgecolor="lightgrey",
        lw=0.5 * scale / 5,
        transform=trans,
        clip_on=False,
    )
    ax.add_patch(patch)

    return 0


def make_circle_plot(
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
    lagna_lambda,
    nakṣatra_names_file="nakshatra_names.tex",
    rāśi_names_file="rashi_names.tex",
    plotfile="panchanga_at_test_time.pdf",
):
    ## Sanskrit typesetting using XeLaTeX ##
    mpl.use("pgf")

    ## TeX preamble ##
    preamble = r"""\usepackage{fontspec}
               \usepackage{polyglossia}
               \usepackage[T1]{fontenc}
               \usepackage{tikz}
               \setmainlanguage{english}
               \setotherlanguages{sanskrit}
               \newfontfamily\devanagarifont[Script=Devanagari]{Sanskrit 2003}
               \XeTeXgenerateactualtext 1
               \newcommand{\sam}[1]{\begin{sanskrit}#1\end{sanskrit}}"""

    params = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "xelatex",
        "pgf.preamble": preamble,
    }
    mpl.rcParams.update(params)

    plt.style.use("dark_background")

    observing_location = EarthLocation.of_address(location)

    fig, ax = plt.subplots(subplot_kw={"projection": "polar"}, figsize=(7, 5))

    ## Setup circle plot with nakṣatra names and rāśi names

    nakṣatra_edges = np.linspace(ayanāṃśa, ayanāṃśa + 26 * nakṣatra_extent, 27) % 360
    nakṣatra_centers = nakṣatra_edges + nakṣatra_extent / 2
    nakṣatra_edge_labels = []

    nakṣatra_names_tab = Table.read(nakṣatra_names_file, format="latex")
    nakṣatra_names = nakṣatra_names_tab["names"].data

    for i in range(len(nakṣatra_edges)):
        coord = nakṣatra_edges[i]
        if coord > 360:
            coord = nakṣatra_edges[i] - 360
        deg_coord = int(coord)
        min_coord = 60 * (coord % 1)
        nakṣatra_edge_labels.append(
            r"%d$^{{\circ}}$\,%d$^{{\prime}}$" % (deg_coord, min_coord)
        )
        ax.plot(np.deg2rad([coord, coord]), [1, 2.5], "grey", alpha=0.3)

        rotation_val = nakṣatra_centers[i]
        color_val = "grey"
        alpha_val = 0.3

        if 90 < nakṣatra_centers[i] < 265:
            rotation_val = nakṣatra_centers[i] + 180

        if nakṣatra_names[i] == nakṣatra:
            color_val = "white"
            alpha_val = 1.0

        ax.text(
            np.deg2rad(nakṣatra_centers[i]),
            1.5,
            nakṣatra_names[i],
            rotation=rotation_val,
            alpha=alpha_val,
            color=color_val,
            ha="center",
            va="center",
        )

    rāśi_edges = np.linspace(ayanāṃśa, ayanāṃśa + (rāśi_extent * 11), 12)
    rāśi_centers = rāśi_edges + rāśi_extent / 2
    rāśi_names_tab = Table.read(rāśi_names_file, format="latex")
    rāśi_names = rāśi_names_tab["names"].data

    for i in range(len(rāśi_edges)):
        ax.plot(np.deg2rad([rāśi_edges[i], rāśi_edges[i]]), [0, 1], "grey", alpha=0.3)

        rotation_val = rāśi_centers[i]
        color_val = "grey"
        alpha_val = 0.3

        if 90 < rāśi_centers[i] < 265:
            rotation_val = rāśi_centers[i] + 180

        if (
            rāśi_edges[np.floor((sun_lambda - ayanāṃśa) / 30).astype(int)]
            == rāśi_edges[i]
        ):
            color_val = "white"
            alpha_val = 1
        ax.text(
            np.deg2rad(rāśi_centers[i]),
            0.65,
            rāśi_names[i],
            rotation=rotation_val,
            alpha=alpha_val,
            color=color_val,
            ha="center",
            va="center",
        )

    ### Plotting grahas
    ### Plot Sun
    plot_sun(np.array([np.deg2rad(sun_lambda), 0.9]), 0.08, fig, ax)

    ## Plot moon with phase ##
    plot_moon_phase(tithi, np.array([np.deg2rad(moon_lambda), 1.9]), 0.08, fig, ax)

    ### Inner grahas ###

    mercury_lambda = get_body(
        "mercury", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    venus_lambda = get_body(
        "venus", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    if sun_lambda > mercury_lambda:
        mercury_angle_to_sun = 360 - (sun_lambda - mercury_lambda)
    else:
        mercury_angle_to_sun = mercury_lambda - sun_lambda

    if sun_lambda > venus_lambda:
        venus_angle_to_sun = 360 - (sun_lambda - venus_lambda)
    else:
        venus_angle_to_sun = venus_lambda - sun_lambda

    plot_inner_graha_phase(
        "budha", mercury_angle_to_sun, [np.deg2rad(mercury_lambda), 1.9], 0.05, fig, ax
    )
    plot_inner_graha_phase(
        "shukra", venus_angle_to_sun, [np.deg2rad(venus_lambda), 1.9], 0.05, fig, ax
    )

    ### Outer grahas ###

    mars_lambda = get_body(
        "mars", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    jupiter_lambda = get_body(
        "jupiter", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    saturn_lambda = get_body(
        "saturn", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    plot_outer_graha("mangala", [np.deg2rad(mars_lambda), 1.9], 0.05, fig, ax)
    plot_outer_graha("guru", [np.deg2rad(jupiter_lambda), 1.9], 0.05, fig, ax)
    plot_outer_graha("shani", [np.deg2rad(saturn_lambda), 1.9], 0.05, fig, ax)

    ##  rahu and ketu  ##

    # This eclipse happens when moon was in Rahu's postion
    known_eclipse_time = Time("2023-04-20 04:16:49", format="iso", scale="utc")

    moon_lambda_known_eclipse = get_body(
        "moon", known_eclipse_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    rahu_lambda, ketu_lambda = calc_rahu_ketu_pos(
        known_eclipse_time, test_date_utc_time, moon_lambda_known_eclipse
    )

    rahu_lambda, ketu_lambda = ketu_lambda, rahu_lambda

    plot_rahu([np.deg2rad(rahu_lambda), 1.9], 5, fig, ax)
    plot_ketu([np.deg2rad(ketu_lambda), 1.9], 5, fig, ax)

    ## lagna ##

    ax.text(
        np.deg2rad(lagna_lambda),
        0.9,
        r"\sam{लग्न}",
        color="yellow",
        rotation=lagna_lambda + 90,
        ha="center",
        va="center",
    )

    ################ legend for grahas #####################

    inv = ax.transData.inverted()

    xpos = 1100
    ypos = 450

    sun_legend_pos = inv.transform((xpos, ypos))
    sun_label_pos = inv.transform((xpos + 45, ypos - 5))
    moon_legend_pos = inv.transform((xpos, ypos - 50))
    moon_label_pos = inv.transform((xpos + 45, ypos - 55))
    mercury_legend_pos = inv.transform((xpos, ypos - 100))
    mercury_label_pos = inv.transform((xpos + 45, ypos - 105))
    venus_legend_pos = inv.transform((xpos, ypos - 150))
    venus_label_pos = inv.transform((xpos + 45, ypos - 155))
    mars_legend_pos = inv.transform((xpos, ypos - 200))
    mars_label_pos = inv.transform((xpos + 45, ypos - 205))
    jupiter_legend_pos = inv.transform((xpos, ypos - 250))
    jupiter_label_pos = inv.transform((xpos + 45, ypos - 255))
    saturn_legend_pos = inv.transform((xpos, ypos - 300))
    saturn_label_pos = inv.transform((xpos + 45, ypos - 305))
    rahu_legend_pos = inv.transform((xpos, ypos - 350))
    rahu_label_pos = inv.transform((xpos + 45, ypos - 355))
    ketu_legend_pos = inv.transform((xpos, ypos - 400))
    ketu_label_pos = inv.transform((xpos + 45, ypos - 405))

    plot_sun(sun_legend_pos, 0.08, fig, ax)
    plt.text(sun_label_pos[0], sun_label_pos[1], "\sam{सूर्यः } ")
    plot_moon_phase(tithi, moon_legend_pos, 0.08, fig, ax)
    plt.text(moon_label_pos[0], moon_label_pos[1], "\sam{चन्द्रः } ")
    plot_inner_graha_phase(
        "budha", mercury_angle_to_sun, mercury_legend_pos, 0.05, fig, ax
    )
    plt.text(mercury_label_pos[0], mercury_label_pos[1], "\sam{बुधः } ")
    plot_inner_graha_phase(
        "shukra", venus_angle_to_sun, venus_legend_pos, 0.05, fig, ax
    )
    plt.text(venus_label_pos[0], venus_label_pos[1], "\sam{शुक्रः } ")
    plot_outer_graha("mangala", mars_legend_pos, 0.05, fig, ax)
    plt.text(mars_label_pos[0], mars_label_pos[1], "\sam{मङ्गलः } ")
    plot_outer_graha("guru", jupiter_legend_pos, 0.05, fig, ax)
    plt.text(jupiter_label_pos[0], jupiter_label_pos[1], "\sam{बृहस्पतिः } ")
    plot_outer_graha("shani", saturn_legend_pos, 0.05, fig, ax)
    plt.text(saturn_label_pos[0], saturn_label_pos[1], "\sam{शनैश्चरः } ")

    plot_rahu(rahu_legend_pos, 5, fig, ax)
    plt.text(rahu_label_pos[0], rahu_label_pos[1], "\sam{राहुः} ")
    plot_ketu(ketu_legend_pos, 5, fig, ax)
    plt.text(ketu_label_pos[0], ketu_label_pos[1], "\sam{केतुः} ")

    ################ legend for panchanga #####################
    xpos_panch = -550
    ypos_panch = 500
    plt.text(
        inv.transform((xpos_panch, ypos_panch))[0],
        inv.transform((xpos_panch, ypos_panch - 100))[1],
        r"\sam{" + location + "}",
    )

    plt.text(
        inv.transform((xpos_panch, ypos_panch - 30))[0],
        inv.transform((xpos_panch, ypos_panch - 130))[1],
        r"\sam{" + date_str.split(" ")[0] + "}",
    )
    plt.text(
        inv.transform((xpos_panch, ypos_panch - 60))[0],
        inv.transform((xpos_panch, ypos_panch - 160))[1],
        r"\sam{" + date_str.split(" ")[1] + "}",
    )

    plt.text(
        inv.transform((xpos_panch + 50, ypos_panch - 160))[0],
        inv.transform((xpos_panch + 50, ypos_panch - 160))[1],
        r"\sam{पञ्चाङ्ग }",
        bbox=dict(facecolor="none", edgecolor="white"),
    )

    plt.text(
        inv.transform((xpos_panch, ypos_panch - 220))[0],
        inv.transform((xpos_panch, ypos_panch - 220))[1],
        tithi_name,
    )

    plt.text(
        inv.transform((xpos_panch, ypos_panch - 260))[0],
        inv.transform((xpos_panch, ypos_panch - 260))[1],
        vāra,
    )

    plt.text(
        inv.transform((xpos_panch, ypos_panch - 300))[0],
        inv.transform((xpos_panch, ypos_panch - 300))[1],
        nakṣatra + r"\hspace{5pt}" + pāda,
    )

    plt.text(
        inv.transform((xpos_panch, ypos_panch - 340))[0],
        inv.transform((xpos_panch, ypos_panch - 340))[1],
        yoga,
    )
    plt.text(
        inv.transform((xpos_panch, ypos_panch - 380))[0],
        inv.transform((xpos_panch, ypos_panch - 380))[1],
        karaṇa,
    )
    # plt.text(inv.transform((-200, 220))[0], inv.transform((-200, 220))[1], final_yoga)

    plt.title(r"\sam{ॐ }", fontsize=20)

    # ---------------------------------------

    ax.set_xticks(np.deg2rad(nakṣatra_edges))
    ax.set_xticklabels(nakṣatra_edge_labels, fontsize=5)
    ax.set_rmax(2)
    ax.set_yticks([1])
    ax.set_yticklabels([""])
    ax.grid(axis="x")

    plt.savefig(plotfile)

    return 0


def make_jatakam_plot(
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
    lagna_lambda,
    rāśi_names_file="rashi_names.tex",
    plotfile="jatakam_at_test_time.pdf",
):
    plt.style.use("default")
    ## Sanskrit typesetting using XeLaTeX ##
    mpl.use("pgf")

    ## TeX preamble ##
    preamble = r"""\usepackage{fontspec}
               \usepackage{polyglossia}
               \usepackage[T1]{fontenc}
               \usepackage{tikz}
               \setmainlanguage{english}
               \setotherlanguages{sanskrit}
               \newfontfamily\devanagarifont[Script=Devanagari]{Sanskrit 2003}
               \XeTeXgenerateactualtext 1
               \newcommand{\sam}[1]{\begin{sanskrit}#1\end{sanskrit}}"""

    params = {
        "font.family": "serif",
        "text.usetex": True,
        "pgf.rcfonts": False,
        "pgf.texsystem": "xelatex",
        "pgf.preamble": preamble,
    }
    mpl.rcParams.update(params)

    fig, ax = plt.subplots()
    divider = make_axes_locatable(ax)
    pañcāṅga_ax = divider.append_axes("left", 0.5, sharey=ax)

    fig.suptitle(r"\sam{{ॐ }}", color="b")
    ax.vlines([0, 2, 6, 8], ymin=0, ymax=8, color="b")
    ax.vlines([4, 4], ymin=[0, 6], ymax=[2, 8], color="b")
    ax.hlines([0, 2, 6, 8], xmin=0, xmax=8, color="b")
    ax.hlines([4, 4], xmin=[0, 6], xmax=[2, 8], color="b")

    coord_for_rāśi = [
        (3, 7),
        (5, 7),
        (7, 7),
        (7, 5),
        (7, 3),
        (7, 1),
        (5, 1),
        (3, 1),
        (1, 1),
        (1, 3),
        (1, 5),
        (1, 7),
    ]
    rāśi_tab = Table.read(rāśi_names_file, format="latex")

    rāśi_names = rāśi_tab["names"].data
    ### Inner grahas ###
    observing_location = EarthLocation.of_address(location)
    mercury_lambda = get_body(
        "mercury", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    venus_lambda = get_body(
        "venus", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    ### Outer grahas ###

    mars_lambda = get_body(
        "mars", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    jupiter_lambda = get_body(
        "jupiter", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    saturn_lambda = get_body(
        "saturn", test_date_utc_time, location=observing_location
    ).geocentrictrueecliptic.lon.value

    ##  rahu and ketu  ##

    # This eclipse happens when moon was in Rahu's postion
    known_eclipse_time = Time("2023-04-20 04:16:49", format="iso", scale="utc")

    moon_lambda_known_eclipse = get_body(
        "moon", known_eclipse_time, location=observing_location
    ).geocentrictrueecliptic.lon.value
    rahu_lambda, ketu_lambda = calc_rahu_ketu_pos(
        known_eclipse_time, test_date_utc_time, moon_lambda_known_eclipse
    )

    rahu_lambda, ketu_lambda = ketu_lambda, rahu_lambda

    lst_of_grahas_pos = [
        sun_lambda,
        moon_lambda,
        mercury_lambda,
        venus_lambda,
        mars_lambda,
        jupiter_lambda,
        saturn_lambda,
        rahu_lambda,
        ketu_lambda,
    ]
    graha_tab = Table.read("graha_names.tex", format="latex")

    graha_names = graha_tab["names"].data

    for graha_id in range(len(graha_names)):
        sel_rāśi_id = np.floor(
            (lst_of_grahas_pos[graha_id] - ayanāṃśa) / rāśi_extent
        ).astype(int)
        print(graha_names[graha_id], rāśi_names[sel_rāśi_id])

    ################ legend for panchanga #####################
    pañcāṅga_ax.text(
        0.2,
        7.5,
        r"\sam{" + location + "}",
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        7.0,
        r"\sam{" + date_str.split(" ")[0] + "}",
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        6.5,
        r"\sam{" + date_str.split(" ")[1] + "}",
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        1,
        4.0,
        r"\sam{पञ्चाङ्ग }",
        bbox=dict(facecolor="none", edgecolor="b"),
        color="b",
        va="center",
        ha="center",
    )

    pañcāṅga_ax.text(
        0.2,
        2.5,
        tithi_name,
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        2.0,
        vāra,
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        1.5,
        nakṣatra + r"\hspace{5pt}" + pāda,
        color="b",
        va="center",
    )

    pañcāṅga_ax.text(
        0.2,
        1.0,
        yoga,
        color="b",
        va="center",
    )
    pañcāṅga_ax.text(
        0.2,
        0.5,
        karaṇa,
        color="b",
        va="center",
    )

    ax.set_aspect("equal")
    # pañcāṅga_ax.set_aspect("equal")
    pañcāṅga_ax.axis("off")
    ax.axis("off")
    plt.savefig(plotfile, bbox_inches="tight")

    return 0
