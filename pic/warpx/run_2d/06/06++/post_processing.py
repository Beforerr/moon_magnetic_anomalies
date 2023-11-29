# %%
from parameters import *
from pywarpx import picmi
import unyt
from unyt import unyt_quantity
import yt
from yt import derived_field

yt.enable_parallelism()
yt.set_log_level(40)

_check = 0
_rho_slice_plot = 1
_slice_plot = 1
_particle_plot = 0

simulation_length_yt = unyt_quantity.from_astropy(simulation_length)
width_x = 7 / 8 * simulation_length_yt
center_x = simulation_length_yt / 2
center = unyt.unyt_array([width_x / 2, 0 * width_x, 0 * width_x])
L_y_yt = unyt_quantity.from_astropy(L_y)
width = (width_x, L_y_yt)


from pathlib import Path

subdirectories = ["figures"]
for path in subdirectories:
    Path(path).mkdir(parents=True, exist_ok=True)


# %%
if _check:
    pass

plotfile = "diags/diag1/plt??????"
ts = yt.load(plotfile)

# %%
def slice_plot(
    fields,
    name,
    _Multipanel=False,
    nrows_ncols=None,
):
    for ds in ts.piter():
        slc = yt.SlicePlot(
            ds, "z", fields, origin="native"
        )  # Create a sliceplot object
        slc.set_axes_unit("km")
        if _Multipanel:
            fig = slc.export_to_mpl_figure(nrows_ncols)
            fig.tight_layout()
            fig.savefig("figures/{}_{}.png".format(str(ds), name))
        else:
            slc.save("figures/")


def slice_sub_plot(
    fields,
    name,
    _Multipanel=False,
    nrows_ncols=None,
):
    for ds in ts.piter():
        slc = yt.SlicePlot(
            ds,
            "z",
            fields,
            origin="native",
            center=center,
            width=width,
        )  # Create a sliceplot object
        slc.set_axes_unit("km")
        if _Multipanel:
            fig = slc.export_to_mpl_figure(nrows_ncols)
            fig.tight_layout()
            fig.savefig("figures/{}_{}_sub.png".format(str(ds), name))
        else:
            slc.save("figures/")


if _slice_plot:
    fields = [
        ("boxlib", "Bx"),
        ("boxlib", "By"),
        ("boxlib", "Bz"),
        ("mesh", "magnetic_field_strength"),
    ]
    # slice_plot(fields, "B", _Multipanel=True, nrows_ncols=(2, 2))
    # slice_sub_plot(fields, "B", _Multipanel=True, nrows_ncols=(2, 2))

    fields = [
        ("mesh", "magnetic_field_strength"),
        ("boxlib", "Ex"),
        ("boxlib", "Ey"),
        ("boxlib", "Ez"),
    ]
    _Multipanel = True
    _show = False
    nrows_ncols = (2, 2)
    for ds in ts.piter():
        slc = yt.SlicePlot(
            ds,
            "z",
            fields,
            origin="native",
            center=center,
            width=width,
        )  # Create a sliceplot object
        slc.set_axes_unit("km")
        slc.set_cmap(field=("mesh", "magnetic_field_strength"), cmap="RdBu")

        for field in [("boxlib", "Ex"), ("boxlib", "Ey"), ("boxlib", "Ez")]:
            slc.set_unit(field, "mV/m")
            slc.set_log(field, False)
            slc.set_zlim(field, zmin=(-1000, "mV/m"), zmax=(1000, "mV/m"))
            slc.set_cmap(field=field, cmap="RdBu")

        if _Multipanel:
            fig = slc.export_to_mpl_figure(nrows_ncols)
            fig.tight_layout()
            fig.savefig("figures/{}_BE_sub.png".format(str(ds)))
        if _show:
            fig.show()
        else:
            slc.save("figures/")


# %%

# %%
def rho_slice_plot(
    fields,
    directions,
    name,
    _Multipanel=False,
    nrows_ncols=None,
    add_field=None,
):
    for ds in ts.piter():
        if add_field:
            add_field(ds)
        for direction in directions:
            slc = yt.SlicePlot(ds, direction, fields, origin="native")
            slc.set_axes_unit("km")
            slc.set_cmap(field=("boxlib", "density"), cmap="RdBu")

            if _Multipanel:
                fig = slc.export_to_mpl_figure(nrows_ncols)
                fig.tight_layout()
                fig.savefig(
                    "figures/{}_Slice_{}_{}.png".format(str(ds), direction, name)
                )
            else:
                slc.save("figures/")


def rho_slice_sub_plot(
    fields,
    directions,
    name,
    _Multipanel=False,
    nrows_ncols=None,
    add_field=None,
):
    for ds in ts.piter():
        if add_field:
            add_field(ds)
        for direction in directions:
            slc = yt.SlicePlot(
                ds, direction, fields, origin="native", center=center, width=width
            )
            slc.set_axes_unit("km")
            slc.set_cmap(field=("boxlib", "density"), cmap="RdBu")

            if _Multipanel:
                fig = slc.export_to_mpl_figure(nrows_ncols)
                fig.tight_layout()
                fig.savefig(
                    "figures/{}_Slice_sub_{}_{}.png".format(str(ds), direction, name)
                )
            else:
                slc.save("figures/")


def _density_e(field, data):
    return data["boxlib", "rho_electron"] / -picmi.constants.q_e * unyt.m**-3


def _density_p(field, data):
    return data["boxlib", "rho_proton"] / picmi.constants.q_e * unyt.m**-3


def _density(field, data):
    return (
        data["boxlib", "rho_proton"] / picmi.constants.q_e
        + data["boxlib", "rho_electron"] / picmi.constants.q_e
    ) * unyt.m**-3


if _rho_slice_plot:
    # directions = ["z"]
    # fields = [("boxlib", "rho_electron"), ("boxlib", "rho_proton"), ("boxlib", "rho")]
    # rho_slice_plot(fields, directions, "rho", _Multipanel=True, nrows_ncols=(3, 1))
    directions = ["z"]
    fields = [
        ("boxlib", "density_electron"),
        ("boxlib", "density_proton"),
        ("boxlib", "density"),
    ]

    def add_density_field(ds):
        ds.add_field(
            name=("boxlib", "density_electron"),
            function=_density_e,
            sampling_type="local",
            units="cm**-3",
            display_name=r"\rho_e",
        )
        ds.add_field(
            name=("boxlib", "density_proton"),
            function=_density_p,
            sampling_type="local",
            units="cm**-3",
            display_name=r"\rho_p",
        )
        ds.add_field(
            name=("boxlib", "density"),
            function=_density,
            sampling_type="local",
            units="cm**-3",
            display_name=r"\rho",
        )

    _Multipanel = True
    nrows_ncols = (3, 1)
    _save = False

    for ds in ts.piter():
        add_density_field(ds)
        for direction in directions:
            slc = yt.SlicePlot(
                ds, direction, fields, origin="native", center=center, width=width
            )
            slc.set_axes_unit("km")
            for field in fields:
                slc.set_cmap(field=field, cmap="RdBu")

            if _Multipanel:
                fig = slc.export_to_mpl_figure(nrows_ncols)
                fig.tight_layout()

            if _save:
                if _Multipanel:
                    fig.savefig(
                        "figures/{}_Slice_sub_{}_density.png".format(str(ds), direction)
                    )
                else:
                    slc.save("figures/")

# %%
plotfile = "diags/diag2/plt??????"
ts = yt.load(plotfile)

# %%
def get_species(ds):
    psset = set()
    for ps in ds.particle_types:
        if ps == "all":
            continue
        psset.add(ps)
    return list(psset)


@yt.particle_filter(requires=["particle_position_x"], filtered_type="electron")
def filtered_electron(pfilter, data):
    filter = data[(pfilter.filtered_type, "particle_position_x")] < simulation_length_yt
    return filter


@yt.particle_filter(requires=["particle_position_x"], filtered_type="proton")
def filtered_proton(pfilter, data):
    filter = data[(pfilter.filtered_type, "particle_position_x")] < simulation_length_yt
    return filter


# %%


def particle_plot(filter=False):
    for ds in ts.piter():
        # pslist = get_species(ds)
        if filter:
            ds.add_particle_filter("filtered_electron")
            ds.add_particle_filter("filtered_proton")

        for ps in ds.particle_types:

            if filter:
                if "filtered" not in ps:
                    continue

            ds.index
            ds.field_info[ps, "particle_momentum_x"].take_log = False
            ds.field_info[ps, "particle_momentum_y"].take_log = False
            ds.field_info[ps, "particle_momentum_z"].take_log = False

            name = "figures/{}/".format(ps)

            p = yt.ParticlePlot(
                ds,
                (ps, "particle_position_x"),
                (ps, "particle_position_y"),
                (ps, "particle_weight"),
                origin="native",
                center=center,
                width=width,
            )
            p.save(name)
            p = yt.ParticlePlot(
                ds,
                (ps, "particle_position_x"),
                (ps, "particle_momentum_x"),
                (ps, "particle_weight"),
            )
            p.save(name)
            p = yt.ParticlePlot(
                ds,
                (ps, "particle_position_x"),
                (ps, "particle_momentum_y"),
                (ps, "particle_weight"),
            )
            p.save(name)
            p = yt.ParticlePlot(
                ds,
                (ps, "particle_position_x"),
                (ps, "particle_momentum_z"),
                (ps, "particle_weight"),
            )
            p.save(name)
            # p = yt.ParticlePlot(
            #     ds,
            #     (ps, "particle_momentum_x"),
            #     (ps, "particle_momentum_y"),
            #     (ps, "particle_weight"),
            # )
            # p.save(name)
            # p = yt.ParticlePlot(
            #     ds,
            #     (ps, "particle_momentum_x"),
            #     (ps, "particle_momentum_z"),
            #     (ps, "particle_weight"),
            # )
            # p.save(name)
            # p = yt.ParticlePlot(
            #     ds,
            #     (ps, "particle_momentum_z"),
            #     (ps, "particle_momentum_y"),
            #     (ps, "particle_weight"),
            # )
            # p.save(name)


if _particle_plot:
    # particle_plot(filter=False)
    particle_plot(filter=True)
