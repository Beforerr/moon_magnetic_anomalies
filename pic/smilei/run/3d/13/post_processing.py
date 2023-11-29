# %% [markdown]
# # Solar wind interaction with Lunar magnetic anomalies

# %% [markdown]
# ## Set up environment and configure setting

# %%
# Smilei job
# ! qdel 2538
# ! spack env activate smilei-oneapi && $HOME/Smilei/build/smilei-oneapi/smilei test.py
# ! qsub $HOME/pic/zijin_oneapi.pbs -V
# ! showq && qstat

# %%
import happi
import numpy as np
import yt
from yt import derived_field
from unyt.dimensions import length
import xarray as xr
import holoviews as hv
import hvplot.xarray  # noqa
from holoviews import opts
import matplotlib.pyplot as plt
from utilities import *
from _utils import *

yt.set_log_level(40)

SPP = SmileiPostProcessing()
SPP.print_simulation()
S = SPP.SmileiSimulation

# %% [markdown]
# ### Set metadata

# %%
# Set fields metadata: Ex Ey Ez Bx By Bz Bx_m By_m Bz_m Jx Jy Jz Rho Jx_solar_ion Jy_solar_ion Jz_solar_ion Rho_solar_proton Jx_solar_electron Jy_solar_electron Jz_solar_electron Rho_solar_electron

magnetic_field = {}

magnetic_field["yt"] = ["magnetic_field_x", "magnetic_field_y", "magnetic_field_z"]
magnetic_field["smilei"] = ["Bx", "By", "Bz"]
magnetic_field["units"] = ["dimensionless"] * dimensions
B_labels = [r"$B_{x}$", r"$B_{y}$", r"$B_{z}$"]

electric_field = {}
electric_field["yt"] = ["electric_field_x", "electric_field_y", "electric_field_z"]
electric_field["smilei"] = ["Ex", "Ey", "Ez"]
electric_field["units"] = ["dimensionless"] * dimensions
electric_field["labels"] = [r"$E_{x}$", r"$E_{y}$", r"$E_{z}$"]

rho_fields = {}
rho_fields['yt'] = ["Rho_solar_electron", "Rho_solar_proton"]
rho_fields['smilei'] = ["-Rho_solar_electron", "Rho_solar_proton"]
rho_fields['units'] = ["dimensionless"] * 2  # rho_units = ['kg/m**3', 'kg/m**3']
rho_fields['labels'] = [r"$\rho_{solar-electron}$", r"$\rho_{solar-proton}$"]

fields_smilei = [
    *magnetic_field["smilei"],
    *electric_field["smilei"],
    *rho_fields['smilei'],
]
fields_yt = [*magnetic_field["yt"], *electric_field["yt"], *rho_fields['yt']]
fields_units = [*magnetic_field["units"], *electric_field["units"], *rho_fields['units']]


# %%
def _reference_quantities_yt():
    from unyt import unyt_quantity

    global L_r_yt, d_i_yt, V_r_yt, T_r_yt
    d_i_yt = unyt_quantity.from_astropy(d_i)

    L_r_yt = unyt_quantity.from_astropy(L_r)
    V_r_yt = unyt_quantity.from_astropy(V_r)
    T_r_yt = unyt_quantity.from_astropy(T_r)
    pass


# %%
bbox = [
    [0, grid_length[0]],
    [-grid_length[1] / 2, grid_length[1] / 2],
    [-grid_length[2] / 2, grid_length[2] / 2],
] * L_r

# timesteps = S.Field(0, fields_smilei[0]).getTimesteps().astype(int)

# %%
# def getFrequency(diag=0):
#     """Get the frequency of a probe"""
#     fields = S.probeInfo(diag)["fields"]
#     field = fields[0]
#     times = S.Probe(diag, field).getTimes() * T_r
#     return 1 / (times[1] - times[0])


# # %%
# def plot_local_probe(location, fields, probeNumber=0, path="figures/"):
#     from scipy import signal

#     Diag = {}
#     for field in fields:
#         Diag[field] = S.Probe(0, field, subset={"axis1": location}, label=field)
#     happi.multiPlot(*Diag.values())
#     data = {}
#     ds = xr.Dataset()
#     fs = getFrequency().value

#     def _plt():
#         for field in fields:
#             Diag = S.Probe(probeNumber, field, subset={"axis1": location}, label=field)
#             data = np.array(Diag.getData())
#             f, t, Sxx = signal.spectrogram(data, fs)
#             ic(f.shape, t.shape, Sxx.shape)
#             plt.figure(tight_layout=True)
#             plt.pcolormesh(t, f, Sxx, shading="gouraud")
#         return

#     def _pplt():
#         fig, axs = pplt.subplots(ncols=len(fields))
#         for index, field in enumerate(fields):
#             Diag = S.Probe(probeNumber, field, subset={"axis1": location}, label=field)
#             data = np.array(Diag.getData())
#             f, t, Sxx = signal.spectrogram(data, fs)
#             axs[index].pcolormesh(t, f, Sxx, N=200, colorbar="t")

#     def _xr():
#         for field in fields:
#             Diag = S.Probe(probeNumber, field, subset={"axis1": location}, label=field)
#             data = np.array(Diag.getData())
#             f, t, Sxx = signal.spectrogram(data, fs)
#             time = xr.DataArray(t, dims=["Time"], attrs={"units": "s"})
#             frequency = xr.DataArray(f, dims=["Frequency"], attrs={"units": "Hz"})
#             ds[field] = xr.DataArray(
#                 Sxx,
#                 name=field,
#                 coords=[frequency, time],
#                 attrs={"long_name": field, "units": ""},
#             )
#         return ds

#     def _xr_pplt():
#         ds = _xr()
#         fig, axs = pplt.subplots(ncols=len(fields), suptitle="Spectrogram")
#         for index, field in enumerate(fields):
#             axs[index].contourf(ds[field], N=200, colorbar="t")
#         fig.save(path + "B/spectrogram")

#     def _hvplot():
#         ds = _xr()
#         f = {}
#         opts.defaults(opts.Image(colorbar=False))
#         for index, field in enumerate(fields):
#             f[field] = ds[field].hvplot()
#             f[field].opts(title=field)
#             if index != 0:
#                 f[field].opts(yaxis="bare")
#         layout = hv.Layout(list(f.values()))
#         hvplot.save(layout, path + "B/spectrogram.html", resources=INLINE)
#         return layout

#     return _xr_pplt()
#     return _hvplot()
#     return _plt()
#     return _pplt()


fields = magnetic_field["smilei"]

locations = [
    # cell_length[0],
    grid_length[0] * 1 / 2,
    # grid_length[0] - cell_length[0],
]
l = [plot_local_probe(location, fields) for location in locations]
# l[0]

# %%
def plot_temp(timestep=0):
    for species in solar_species:
        diag_w = species.name + "_weight"
        weight = S.ParticleBinning(diag_w, sum={"y": "all", "z": "all"})
        data = weight.getData(timestep=timestep)
        p = {}
        for weight in ["weight_vx_px", "weight_vy_py", "weight_vz_pz"]:

            diag_p = species.name + "_" + weight
            p[weight] = S.ParticleBinning(diag_p, sum={"y": "all", "z": "all"})
            temp =S.ParticleBinning(diag_p, sum={"y": "all", "z": "all"})
            plt.figure()
            p[weight].plot(timestep=timestep)
            temp.plot()
    pass


plot_temp()
# %%
def TrackParticles(S):
    """
    >>> TrackParticles(S)
    """
    track = {}
    for species in S.namelist.Species:
        track[species.name] = S.TrackParticles(
            species=species.name, axes=["x", "y", "z", "px", "py", "pz", "Id"]
        )
        # track[species].plot()
    pass
    # track[species.name] = S.TrackParticles(species =species.name,axes = ["x","px"])


TrackParticles(S)

#%%
def getParticleData(timestep=None):
    """ """
    bbox = [
        [0, grid_length[0]],
        [0, grid_length[1]],
        [0, grid_length[2]],
    ]
    raw_data = {}
    not_nan = {}
    particle_field_list = [
        "particle_mass",
        "particle_position_x",
        "particle_position_y",
        "particle_position_z",
        "particle_velocity_x",
        "particle_velocity_y",
        "particle_velocity_z",
    ]
    attribute_list = [
        "w",
        "x",
        "y",
        "z",
        "px",
        "py",
        "pz",
        "Ex",
        "Ey",
        "Ez",
        "Bx",
        "By",
        "Bz",
        "Id",
    ]
    for species in solar_species:
        Diag = S.TrackParticles(
            species=species.name,
            axes=attribute_list,
        )
        timesteps = Diag.getTimesteps()
        if timestep == None:
            timestep = timesteps[-1]
        raw_data[species.name] = Diag.getData(timestep=timestep)
        temp = raw_data[species.name][attribute_list[0]][0]
        not_nan[species.name] = ~np.isnan(temp)
        for p in ["px", "py", "pz"]:
            raw_data[species.name][p][0] = (
                raw_data[species.name][p][0] / species.mass_norm
            )

    yt_data = {
        (species.name, field): raw_data[species.name][attribute][0][
            not_nan[species.name]
        ]
        for (field, attribute) in zip(particle_field_list, attribute_list)
        for species in solar_species
    }

    ds = yt.load_particles(
        yt_data,
        length_unit=L_r_yt,
        velocity_unit=V_r_yt,
        time_unit=T_r_yt,
        bbox=bbox,
        sim_time=timestep,
        periodicity=(False, False, False),
        unit_system="mks",
    )
    ds.timestep = timestep
    ds.unit_registry.add("d_i", d_i.value, length, tex_repr="d_i")
    return ds


# %%
def particlePlot_yt(ds):
    p = yt.ParticlePlot(
        ds,
        ("all", "particle_position_x"),
        ("all", "particle_position_y"),
        ("all", "particle_mass"),
    )
    p.set_unit(("all", "particle_position_x"), "d_i")
    p.set_unit(("all", "particle_position_y"), "d_i")
    p.save()

    for species in solar_species:
        for velocity in [
            "particle_velocity_x",
            "particle_velocity_y",
            "particle_velocity_z",
        ]:
            p = yt.ParticlePlot(
                ds,
                (species.name, "particle_position_x"),
                (species.name, velocity),
                (species.name, "particle_mass"),
            )
            p.set_unit((species.name, "particle_position_x"), "d_i")
            p.set_unit((species.name, velocity), "c")
            p.set_xlim(0, 0.625)
            p.save(species.name)


# %% [markdown]
# ## Set up `yt`

# %%
from yt import derived_field
from unyt.dimensions import length

# %% [markdown]
# ### Set derived_field

# %%
@derived_field(
    name=("stream", "magnetic_field_strength"),
    units="dimensionless",
    sampling_type="cell",
)
def _magnetic_field_strength(field, data):
    return np.sqrt(
        data["stream", "magnetic_field_x"] ** 2
        + data["stream", "magnetic_field_y"] ** 2
        + data["stream", "magnetic_field_z"] ** 2
    )


@derived_field(
    name=("stream", "electric_field_strength"),
    units="dimensionless",
    sampling_type="cell",
)
def _magnetic_field_strength(field, data):
    return np.sqrt(
        data["stream", "electric_field_x"] ** 2
        + data["stream", "electric_field_y"] ** 2
        + data["stream", "electric_field_z"] ** 2
    )


# %% [markdown]
# ### Load data

# %%
def getData(timestep):
    data = {}
    for field_yt, field_smilei in zip(fields_yt, fields_smilei):
        data[field_yt] = S.Field(0, field_smilei).getData(timestep=timestep)[0]
    domain_dimensions = data[fields_yt[0]].shape
    data = {k: (v, u) for (k, v), u in zip(data.items(), fields_units)}
    ds = yt.load_uniform_grid(
        data,
        domain_dimensions,
        length_unit="m",
        bbox=bbox,
        nprocs=24,
        sim_time=timestep,
        periodicity=(False, False, False),
        unit_system="mks",
    )
    ds.timestep = timestep
    ds.unit_registry.add("d_i", d_i.value, length, tex_repr="d_i")
    return ds


# %% [markdown]
# ## Plot Fields

# %% [markdown]
# ### Plot 2-D species density profiles
#
# Two-dimensional (top) electron and (bottom) ion charge density profiles, scaled to the initial density, $n_{sw}$ ,and along the dipole axis (Y direction) at z = 0 after the simulation has reached quasi-steady state. The solar wind is flowing perpendicular (in the −X direction) to the lunar surface. Superimposed in black are magnetic field lines.

# %%
def plot_density_profile_2D_yt(ds, axis="z"):
    slc = yt.SlicePlot(
        ds, axis, [("stream", field) for field in rho_fields['yt']], origin="native"
    )
    slc.set_axes_unit("d_i")
    slc.annotate_timestamp()
    for field, label in zip(rho_fields['yt'], rho_fields['labels']):
        slc.set_colorbar_label(("stream", field), label)
        slc.set_log(("stream", field), False)
        slc.set_cmap(("stream", field), "doom")
    if axis == "y":
        fig = slc.export_to_mpl_figure((2, 1))
    else:
        fig = slc.export_to_mpl_figure((1, 2))

    fig.tight_layout()
    fig.savefig("figures/" + "rho/yt_rho_" + axis + "_" + str(ds.timestep).zfill(10))
    return fig


# %%
def plot_density_profile_2D_happi():
    def rho_transform(rho):
        return rho / S.namelist.n_solar

    # subset
    Rho_solar_electron_z0 = S.Field(
        0,
        "-Rho_solar_electron",
        subset={"z": grid_length[2] / 2},
        data_transform=rho_transform,
        title=r"$\rho_{solar-electron}$",
    )
    Rho_solar_ion_z0 = S.Field(
        0,
        "Rho_solar_proton",
        subset={"z": grid_length[2] / 2},
        data_transform=rho_transform,
        title=r"$\rho_{solar-proton}$",
    )
    happi.multiPlot(
        Rho_solar_electron_z0,
        Rho_solar_ion_z0,
        shape=[1, 2],
        saveAs="figures/rho/rho_sub_",
    )
    # average
    Rho_solar_electron_z0 = S.Field(
        0,
        "-Rho_solar_electron",
        average={"z": "all"},
        data_transform=rho_transform,
        title=r"$\rho_{solar-electron}$",
    )
    Rho_solar_ion_z0 = S.Field(
        0,
        "Rho_solar_proton",
        average={"z": "all"},
        data_transform=rho_transform,
        title=r"$\rho_{solar-proton}$",
    )
    happi.multiPlot(
        Rho_solar_electron_z0,
        Rho_solar_ion_z0,
        shape=[1, 2],
        saveAs="figures/rho/rho_avg_",
    )


# %% [markdown]
# ### Plot 1-D species density profiles
# Profiles along the direction parallel to the solar wind flow and through the center of the dipole. The upper panel presents the density profiles, normalized to the initial density $n_{sw}$. The remaining panels hold the magnetic and kinetic pressure profiles for the electron (middle) and ion (bottom) populations, in code units.

# %%
def plot_density_profile_1D_happi(timestep):
    def rho_transform(rho):
        return rho / S.namelist.n_solar

    data = {}
    for field, label in zip(rho_fields['smilei'], rho_fields['labels']):
        data[field] = S.Field(
            0,
            field,
            timesteps=timestep,
            label=label,
            subset={"z": grid_length[2] / 2, "y": grid_length[2] / 2},
            data_transform=rho_transform,
            xlabel=r"Distance above the surface ($d_i$)",
            xmin=grid_length[0],
            xmax=0,
            ylabel=r"$\rho$",
        )
    return happi.multiPlot(
        *[data[field] for field in rho_fields['smilei']], saveAs="figures/rho/rho_sub_1D_"
    )


# %% [markdown]
# ### Plot magnetic field

# %%
def plot_mag_field_yt(ds, axis="z"):
    slc = yt.SlicePlot(
        ds,
        axis,
        [("stream", "magnetic_field_strength")]
        + [("stream", field) for field in magnetic_field["yt"]],
        origin="native",
    )
    slc.set_axes_unit("d_i")
    slc.set_colorbar_label(("stream", "magnetic_field_strength"), r"$B_m$")
    for field, label in zip(magnetic_field["yt"], B_labels):
        slc.set_colorbar_label(("stream", field), label)
        slc.set_log(("stream", field), False)
    fig = slc.export_to_mpl_figure((2, 2))
    fig.tight_layout()
    fig.savefig("figures/B/yt_B_" + str(ds.timestep).zfill(10))
    return fig


# %%
def plot_mag_field_2D_happi():
    plt.figure(figsize=(12, 12), tight_layout=True)
    B_0 = S.Field(
        0, "(Bx**2+By**2+Bz**2)**0.5", subset={"z": grid_length[2] / 2}, data_log=True
    )
    B_x0 = S.Field(0, "Bx", subset={"z": grid_length[2] / 2}, vsym=True)
    B_y0 = S.Field(0, "By", subset={"z": grid_length[2] / 2}, vsym=True)
    B_z0 = S.Field(0, "Bz", subset={"z": grid_length[2] / 2}, vsym=True)
    happi.multiPlot(
        B_0,
        B_x0,
        B_y0,
        B_z0,
        shape=[2, 2],
        saveAs="figures/B/smilei_B_2D_",
    )


def plot_mag_field_1D_happi():
    plt.figure(figsize=(12, 12), tight_layout=True)
    B_0 = S.Field(
        0,
        "(Bx**2+By**2+Bz**2)**0.5",
        subset={"y": grid_length[2] / 2, "z": grid_length[2] / 2},
        data_log=True,
    )
    B_x0 = S.Field(
        0, "Bx", subset={"y": grid_length[2] / 2, "z": grid_length[2] / 2}, vsym=True
    )
    B_y0 = S.Field(
        0, "By", subset={"y": grid_length[2] / 2, "z": grid_length[2] / 2}, vsym=True
    )
    B_z0 = S.Field(
        0, "Bz", subset={"y": grid_length[2] / 2, "z": grid_length[2] / 2}, vsym=True
    )
    happi.multiPlot(
        B_0,
        B_x0,
        B_y0,
        B_z0,
        shape=[2, 2],
        saveAs="figures/B/smilei_B_1D_",
    )


# %% [markdown]
# #### Plot magnetic field streamlines

# %%
def plot_streamlines(ds):
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from yt.visualization.api import Streamlines

    # Define c: the center of the box, N: the number of streamlines,
    # scale: the spatial scale of the streamlines relative to the boxsize,
    # and then pos: the random positions of the streamlines.
    c = ds.domain_center
    N = 100
    scale = ds.domain_width[0]
    pos_dx = np.random.random((N, 3)) * scale - scale / 2.0
    pos = c + pos_dx

    streamlines = Streamlines(
        ds,
        pos,
        ("stream", "Bx"),
        ("stream", "By"),
        ("stream", "Bz"),
    )
    streamlines.integrate_through_volume()

    # Create a 3D plot, trace the streamlines through the 3D volume of the plot
    fig = plt.figure()
    ax = Axes3D(fig, auto_add_to_figure=False)
    fig.add_axes(ax)
    for stream in streamlines.streamlines:
        stream = stream[np.all(stream != 0.0, axis=1)]
        ax.plot3D(stream[:, 0], stream[:, 1], stream[:, 2])

    # Save the plot to disk.
    plt.savefig("streamlines.png")


# %% [markdown]
# ### Plot electric field

# %%
def plot_electric_field_yt(ds, axis="z"):
    slc = yt.SlicePlot(
        ds,
        axis,
        [("stream", "electric_field_strength")]
        + [("stream", field) for field in electric_field["yt"]],
        origin="native",
    )
    slc.set_axes_unit("d_i")
    for field, label in zip(electric_field["yt"], electric_field["labels"]):
        slc.set_colorbar_label(("stream", field), label)
        slc.set_log(("stream", field), False)
    fig = slc.export_to_mpl_figure((2, 2))
    fig.tight_layout()
    fig.savefig("figures/E/yt_E_" + str(ds.timestep).zfill(10))
    return fig


# %% [markdown]
# ## Plot gyroradius

# %%
from astropy import units as u
from astropy.constants import m_e
from plasmapy.formulary import gyroradius


def plot_gyroradii_1D(timestep):
    pxx = S.ParticleBinning("pxx", sum={"x": "all", "y": "all"})
    pyy = S.ParticleBinning("pyy", sum={"x": "all", "y": "all"})
    v_perp = 0.0
    m = 0.0
    q = 0.0
    B = 0.0

    r_g_solar_electron = gyroradius(B, electron, v_perp)
    r_g_solar_ion = gyroradius(B, ion, v_perp)
    pass


# %% [markdown]
# ## Plot plasma beta

# %%
import astropy.units as u
from plasmapy.formulary import beta


def plot_beta_1D(timestep):
    """
    Plot the
    """
    T = u.K
    n = u.m**-3
    B = u.T
    beta(T, n, B)
    pass


# %% [markdown]
# ## Plot velocity distribution

# %%
def plot_electric_drift_1D(timestep):
    pass


def plot_magnetic_drift_1D(timestep):
    pass


# %%
def plot_particle_distribution_vx_vy_vz_happi(timestep=0, center=0):
    """use `happi` to visualize particle distribution
    >>> plot_particle_distribution_vx_vy_vz_happi()
    """
    import proplot as pplt
    import matplotlib.pyplot as plt

    Nt = S.ParticleBinning(0).getTimesteps()[-1]  # the last timestep
    diag_numbers = S.getDiags("ParticleBinning")[0][-1] + 1
    diag_numbers = len(S._diag_numbers["ParticleBinning"])
    for diag in range(diag_numbers):
        f_initial = S.ParticleBinning(
            diag,
            # data_log=True,
            timesteps=0,
            label="initial",
            sum={"x": "all", "y": "all", "z": "all"},
        )
        f_final = S.ParticleBinning(
            diag,
            # data_log=True,
            timesteps=Nt,
            label="final",
            sum={"x": "all", "y": "all", "z": "all"},
        )
    pass


plot_particle_distribution_vx_vy_vz_happi()

# %%
def plot_particle_distribution_x_vx_vy_vz_happi(timestep=None):
    """use `happi` to visualize particle distribution
    >>> plot_particle_distribution_x_vx_vy_vz_happi()
    """
    Nt = S.ParticleBinning(0).getTimesteps()[-1]  # the last timestep
    diag_numbers = len(S._diag_numbers["ParticleBinning"])
    for diag in range(diag_numbers):
        f_initial = S.ParticleBinning(
            diag,
            # data_log=True,
            timesteps=0,
            label="initial",
            sum={"y": "all", "z": "all"},
        )
        f_initial.plot()
        f_final = S.ParticleBinning(
            diag,
            # data_log=True,
            timesteps=Nt,
            label="final",
            sum={"y": "all", "z": "all"},
        )
        f_final.plot()
    pass


# %%
def plot_particle_distribution_x_vx_vy_vz_proplot(timestep=None):
    """use `xarray` and `proplot` to visualize particle distribution
    >>> plot_particle_distribution_x_vx_vy_vz_proplot()
    """
    pass


# %%


ds = xr.Dataset()
ds.coords["x"] = S.ParticleBinning(0, sum={"y": "all", "z": "all"}).getAxis("x")
ds.coords["timestep"] = S.ParticleBinning(
    0, sum={"y": "all", "z": "all"}
).getTimesteps()

species = ["solar_electron", "solar_proton"]
velocities = ["vx", "vy", "vz"]

for species in species:
    for velocity in velocities:
        name = species + "_" + velocity
        ds.coords[name] = S.ParticleBinning(name, sum={"y": "all", "z": "all"}).getAxis(
            velocity
        )
        ds[name + "_distribution"] = (
            ("timestep", "x", name),
            S.ParticleBinning(name, sum={"y": "all", "z": "all"}).getData(),
        )
# %%
def plot_particle_distribution_vx_vy_vz_layout(path="figures/"):
    """
    >>>plot_particle_distribution_vx_vy_vz_layout()
    """
    from bokeh.resources import INLINE

    opts.defaults(opts.Image(invert_xaxis=True, xaxis="bare"))
    f1 = ds["solar_electron_vx_distribution"].hvplot(x="x", y="solar_electron_vx")
    f2 = ds["solar_proton_vx_distribution"].hvplot(x="x", y="solar_proton_vx")
    f3 = ds["solar_electron_vy_distribution"].hvplot(x="x", y="solar_electron_vy")
    f4 = ds["solar_proton_vy_distribution"].hvplot(x="x", y="solar_proton_vy")
    f5 = ds["solar_electron_vz_distribution"].hvplot(x="x", y="solar_electron_vz")
    f6 = ds["solar_proton_vz_distribution"].hvplot(x="x", y="solar_proton_vz")

    f1.opts(colorbar_position="top")
    f2.opts(yaxis="right", colorbar_position="top")
    f3.opts(colorbar=False)
    f4.opts(yaxis="right", colorbar=False)
    f5.opts(colorbar=False, xaxis="bottom")
    f6.opts(yaxis="right", colorbar=False, xaxis="bottom")
    layout = f1 + f2 + f3 + f4 + f5 + f6
    layout.cols(2)
    hvplot.save(layout, path + "particle_distribution.html", resources=INLINE)
    return layout


plot_particle_distribution_vx_vy_vz_layout()

#%%
def plot_particle_distribution_vx_vy_vz_gridspace():
    """
    >>> plot_particle_distribution_vx_vy_vz_gridspace()
    """
    from bokeh.resources import INLINE

    species = ["solar_electron", "solar_proton"]
    velocities = ["vx", "vy", "vz"]

    def particle_distribution(species, velocity):
        name = species + "_" + velocity
        return ds[name + "_distribution"].hvplot(
            x="x", y=name, flip_xaxis=True, colorbar=False
        )

    dict_2D = {(p, f): particle_distribution(p, f) for p in species for f in velocities}

    gridspace = hv.GridSpace(dict_2D, kdims=["species", "velocities"])

    hvplot.save(gridspace, "gridspace.html", resources=INLINE)
    return gridspace


# %% [markdown]
# ## Visualize the tracked macro-particles

# %% [markdown]
# ## Job

# %%
if __name__ == "__main__":

    # temp()
    plot_mag_field_1D_happi()
    plot_mag_field_2D_happi()
# ts = {timestep: getData(timestep) for timestep in timesteps}
# plot_mag_field_yt(ts[0])
# plot_density_profile_2D_yt(ts[0])
# plot_density_profile_1D_happi(timesteps[0])

# SPP.plot_scalar()
# plot_density_profile_2D_happi()

# [plot_mag_field_yt(ds) for ds in ts.values()]
# [plot_electric_field_yt(ds) for ds in ts.values()]

# %time [plot_density_profile_2D_yt(ds, "z") for ds in ts.values()]
# %time [plot_density_profile_2D_yt(ds, "y") for ds in ts.values()]
# %time [plot_density_profile_2D_yt(ds, "x") for ds in ts.values()]
# [plot_density_profile_1D_happi(timestep) for timestep in timesteps]

# %time Parallel(n_jobs=2)(delayed(plot_mag_field_yt)(ds) for ds in ts.values())
# %time Parallel(n_jobs=2)(delayed(plot_density_profile_1D_happi)(timestep) for timestep in timesteps)

# ds=getParticleData()
# particlePlot_yt(ds)


# %% [markdown]
#

# %% [markdown]
# #TODO Parallel visualization in python

# %% [markdown]
# ## Code test

# %%
import doctest

doctest.testmod()

# %%
