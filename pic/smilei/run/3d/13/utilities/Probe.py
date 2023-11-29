from .parameters import *
from happi._SmileiSimulation import SmileiSimulation
import proplot as pplt
import matplotlib.pyplot as plt


def getFrequency(simulation:SmileiSimulation,diag=0):
    """Get the frequency of a probe"""
    fields = simulation.probeInfo(diag)["fields"]
    field = fields[0]
    times = simulation.Probe(diag, field).getTimes() * T_r
    return 1 / (times[1] - times[0])

def plot_local_probe(simulation:SmileiSimulation,location, fields, probeNumber=0, path="figures/"):
    from scipy import signal
    import happi
    import numpy as np
    import xarray as xr
    from holoviews import opts

    Diag = {}
    for field in fields:
        Diag[field] = simulation.Probe(0, field, subset={"axis1": location}, label=field)
    happi.multiPlot(*Diag.values())
    data = {}
    ds = xr.Dataset()
    fs = getFrequency().value

    def _plt():
        for field in fields:
            Diag = simulation.Probe(probeNumber, field, subset={"axis1": location}, label=field)
            data = np.array(Diag.getData())
            f, t, Sxx = signal.spectrogram(data, fs)
            plt.figure(tight_layout=True)
            plt.pcolormesh(t, f, Sxx, shading="gouraud")
        return

    def _pplt():
        fig, axs = pplt.subplots(ncols=len(fields))
        for index, field in enumerate(fields):
            Diag = simulation.Probe(probeNumber, field, subset={"axis1": location}, label=field)
            data = np.array(Diag.getData())
            f, t, Sxx = signal.spectrogram(data, fs)
            axs[index].pcolormesh(t, f, Sxx, N=200, colorbar="t")

    def _xr():
        for field in fields:
            Diag = simulation.Probe(probeNumber, field, subset={"axis1": location}, label=field)
            data = np.array(Diag.getData())
            f, t, Sxx = signal.spectrogram(data, fs)
            time = xr.DataArray(t, dims=["Time"], attrs={"units": "s"})
            frequency = xr.DataArray(f, dims=["Frequency"], attrs={"units": "Hz"})
            ds[field] = xr.DataArray(
                Sxx,
                name=field,
                coords=[frequency, time],
                attrs={"long_name": field, "units": ""},
            )
        return ds

    def _xr_pplt():
        ds = _xr()
        fig, axs = pplt.subplots(ncols=len(fields), suptitle="Spectrogram")
        for index, field in enumerate(fields):
            axs[index].contourf(ds[field], N=200, colorbar="t")
        fig.save(path + "B/spectrogram")

    def _hvplot():
        from bokeh.resources import INLINE
        import holoviews as hv
        import hvplot.xarray  # noqa
        ds = _xr()
        f = {}
        opts.defaults(opts.Image(colorbar=False))
        for index, field in enumerate(fields):
            f[field] = ds[field].hvplot()
            f[field].opts(title=field)
            if index != 0:
                f[field].opts(yaxis="bare")
        layout = hv.Layout(list(f.values()))
        hvplot.save(layout, path + "B/spectrogram.html", resources=INLINE)
        return layout

    return _xr_pplt()
    return _hvplot()
    return _plt()
    return _pplt()

