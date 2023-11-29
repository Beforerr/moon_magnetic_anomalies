# Open SMILEI simulation
from utilities.parameters import *


class SmileiPostProcessing(object):
    """Object for post-processing a Smilei simulation"""

    def __init__(self, results_path="."):
        import happi

        self.SmileiSimulation = happi.Open(
            results_path=results_path, verbose=False, show=True
        )
        self.results_path = results_path
        self.grid_length = self.SmileiSimulation.namelist.Main.grid_length
        from pathlib import Path

        subdirectories = ["figures", "figures/rho", "figures/B", "figures/E", "data"]
        for path in subdirectories:
            Path(results_path + "/" + path).mkdir(parents=True, exist_ok=True)

    def print_simulation(self):
        """Print simulation parameters"""

        S = self.SmileiSimulation
        timestep = S.namelist.Main.timestep
        print("timestep:", timestep, "(code unit)")
        if S.namelist.Main.number_of_timesteps:
            number_of_timesteps = S.namelist.Main.number_of_timesteps
        if S.namelist.Main.simulation_time:
            simulation_time = S.namelist.Main.simulation_time
            number_of_timesteps = int(simulation_time / timestep)
            print("simulation_time:", simulation_time, "(code unit)")
        print("number_of_timesteps:", number_of_timesteps)

        print("geometry:", S.namelist.Main.geometry)  # print the simulation dimensions
        print("number_of_cells:", S.namelist.Main.number_of_cells)
        print("cell_length:", S.namelist.Main.cell_length, "(code unit)")
        print("grid_length:", S.namelist.Main.grid_length, "(code unit)")

        for species in S.namelist.Species:  # Species Information
            print(
                "species " + species.name + " has mass " + str(species.mass),
                " has temperature ",
                species.temperature,
            )
        for external_field in S.namelist.ExternalField:  # Field Information
            print("An external field " + external_field.field + " was applied")

        print("Fields: {0}".format(S.getDiags("Fields")))
        print("Probes: {0}".format(S.getDiags("Probes")))
        print("ParticleBinning: {0}".format(S.getDiags("ParticleBinning")))
        print("Screen: {0}".format(S.getDiags("Screen")))
        print("RadiationSpectrum: {0}".format(S.getDiags("RadiationSpectrum")))

        if "T_r" in globals():
            print("Reference time:", T_r)
            print("timestep:", timestep * T_r)
            print("simulation_time:", simulation_time * T_r)
        if "L_r" in globals():
            print("Reference length:", L_r)
            print("cell_length:", S.namelist.Main.cell_length * L_r)
            print("grid_length:", S.namelist.Main.grid_length * L_r)

    def plot_scalar(self):
        """Plot Scalar field"""
        import happi

        S = self.SmileiSimulation
        Utot = S.Scalar("Utot", label="Total energy")
        Ukin = S.Scalar("Ukin", label="Total kinetic energy")
        Uelm = S.Scalar("Uelm", label="Total electromagnetic energy")
        Ubal = S.Scalar("Ubal", label="Balance energy")
        Ukin_ = {}
        for species in S.namelist.Species:
            Ukin_[species.name] = S.Scalar(
                "Ukin_" + species.name, label="Total " + species.name + " energy"
            )
        happi.multiPlot(
            Utot,
            Ukin,
            Uelm,
            Ubal,
            *Ukin_.values(),
            saveAs=self.results_path + "/" + "figures/scalars",
        )

    # def export_to_vtk(self):
    #     """Export the results in VTK format"""
    #     for field in self.fields_smilei:
    #         self.SmileiSimulation.Field(0, field).toVTK()

    # def export_to_yt(self):
    #     """Save the field file from `Smilei` format to `yt` compatible format."""
    #     import yt

    #     timesteps = (
    #         self.SmileiSimulation.Field(0, self.fields_smilei[0])
    #         .getTimesteps()
    #         .astype(int)
    #     )
    #     for timestep in timesteps:
    #         data = {}
    #         field_types = {}
    #         for field_yt, field_smilei in zip(self.fields_yt, self.fields_smilei):
    #             data[field_yt] = yt.YTArray(
    #                 self.SmileiSimulation.Field(0, field_smilei).getData(
    #                     timestep=timestep
    #                 )[0],
    #                 "dimensionless",
    #             )
    #             field_types[field_yt] = "stream"
    #         ds = {}
    #         yt.save_as_dataset(
    #             ds,
    #             self.results_path + "/" + "data/yt_data_" + str(timestep) + ".h5",
    #             data,
    #             field_types=field_types,
    #         )
