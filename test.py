import matplotlib.pyplot as plt

fig, ax = plt.subplots()

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
rāśi_names = [
    "mesha",
    "vrishabha",
    "mithuna",
    "karka",
    "simha",
    "kanya",
    "tula",
    "vrshchika",
    "dhanusha",
    "makara",
    "kumbh",
    "mina",
]

for i in range(len(rāśi_names)):
    ax.text(
        coord_for_rāśi[i][0],
        coord_for_rāśi[i][1],
        rāśi_names[i],
        ha="center",
        va="center",
        color="b",
    )

ax.set_aspect("equal")
ax.axis("off")
plt.show()
