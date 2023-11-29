# %%
from pywarpx import picmi
import unyt
import yt
from yt import derived_field

yt.enable_parallelism()

_check = 0
_rho_slice_plot = 0
_slice_plot = 1
_particle_plot = 0

# %%
@derived_field(
    name=("boxlib", "density_electron"),
    sampling_type="local",
    units="cm**-3",
    display_name=r"\rho_e",
)
def _density_e(field, data):
    return data["boxlib", "rho_electron"] / -picmi.constants.q_e * unyt.m**-3


@derived_field(
    name=("boxlib", "density_proton"),
    sampling_type="local",
    units="cm**-3",
    display_name=r"\rho_p",
)
def _density_p(field, data):
    return data["boxlib", "rho_proton"] / picmi.constants.q_e * unyt.m**-3


@derived_field(
    name=("boxlib", "density"),
    sampling_type="local",
    units="cm**-3",
    display_name=r"\rho",
)
def _density(field, data):
    return (
        data["boxlib", "rho_proton"] / picmi.constants.q_e
        + data["boxlib", "rho_electron"] / picmi.constants.q_e
    ) * unyt.m**-3


# %%
plotfile = "diags/diag1/plt??????"
ts = yt.load(plotfile)

from pathlib import Path

subdirectories = ["figures"]
for path in subdirectories:
    Path(path).mkdir(parents=True, exist_ok=True)

# %%
if _check:
    pass


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
        if _Multipanel:
            fig = slc.export_to_mpl_figure(nrows_ncols)
            fig.tight_layout()
            fig.savefig("figures/{}_{}.png".format(str(ds), name))
        else:
            slc.save("figures/")


if _slice_plot:
    fields = [
        ("boxlib", "Bx"),
        ("boxlib", "By"),
        ("boxlib", "Bz"),
        ("mesh", "magnetic_field_strength"),
    ]
    slice_plot(fields, "B", _Multipanel=True, nrows_ncols=(2, 2))

    fields = [
        ("boxlib", "Ex"),
        ("boxlib", "Ey"),
        ("boxlib", "Ez"),
    ]
    slice_plot(fields, "E", _Multipanel=True, nrows_ncols=(3, 1))


# %%


def rho_slice_plot(
    fields,
    directions,
    name,
    _Multipanel=False,
    nrows_ncols=None,
):
    for ds in ts.piter():
        for direction in directions:
            slc = yt.SlicePlot(ds, direction, fields, origin="native")

            slc.set_cmap(field=("boxlib", "density"), cmap="RdBu")

            if _Multipanel:
                fig = slc.export_to_mpl_figure(nrows_ncols)
                fig.tight_layout()
                fig.savefig(
                    "figures/{}_Slice_{}_{}.png".format(str(ds), direction, name)
                )
            else:
                slc.save("figures/")


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
    rho_slice_plot(fields, directions, "density", _Multipanel=True, nrows_ncols=(3, 1))


# %%
def get_species(ds):
    psset = set()
    for ps in ds.particle_types:
        if ps == "all":
            continue
        psset.add(ps)
    return list(psset)


# %%
def particle_plot():
    for ds in ts.piter():
        # pslist = get_species(ds)

        for ps in ds.particle_types:
            ds.index
            for direction in ["x", "y", "z"]:
                ds.field_info[
                    ps, "particle_momentum_{}".format(direction)
                ].take_log = False
            name = "figures/{}/".format(ps)

            p = yt.ParticlePlot(
                ds,
                (ps, "particle_position_x"),
                (ps, "particle_position_y"),
                (ps, "particle_weight"),
                origin="native",
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
            p = yt.ParticlePlot(
                ds,
                (ps, "particle_momentum_x"),
                (ps, "particle_momentum_y"),
                (ps, "particle_weight"),
            )
            p.save(name)
            p = yt.ParticlePlot(
                ds,
                (ps, "particle_momentum_x"),
                (ps, "particle_momentum_z"),
                (ps, "particle_weight"),
            )
            p.save(name)
            p = yt.ParticlePlot(
                ds,
                (ps, "particle_momentum_z"),
                (ps, "particle_momentum_y"),
                (ps, "particle_weight"),
            )
            p.save(name)


if _particle_plot:
    particle_plot()
