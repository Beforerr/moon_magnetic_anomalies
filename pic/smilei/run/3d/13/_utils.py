# Open SMILEI simulation
class SmileiPostProcessing(object):
    """Object for post-processing a Smilei simulation"""

    def __init__(
        self, results_path=".", fields_smilei=None, fields_yt=None, fields_units=None
    ):
        import happi
        
        self.SmileiSimulation = happi.Open(results_path=results_path, verbose=False)
        self.results_path = results_path
        self.fields_smilei = fields_smilei
        self.fields_yt = fields_yt
        self.fields_units = fields_units
        # self.timesteps = self.SmileiSimulation.Field(0, fields_smilei[0]).getTimesteps().astype(int)
        self.grid_length = self.SmileiSimulation.namelist.Main.grid_length
        from pathlib import Path
        subdirectories = ["figures", "figures/rho", "figures/B", "figures/E", "data"]
        for path in subdirectories:
            Path(results_path + "/" + path).mkdir(parents=True, exist_ok=True)

    def print_simulation(self):
        """Print simulation parameters"""
        S = self.SmileiSimulation
        timestep = S.namelist.Main.timestep
        print("timestep:", timestep)  # print the timestep
        if S.namelist.Main.number_of_timesteps:
            number_of_timesteps = S.namelist.Main.number_of_timesteps
        if S.namelist.Main.simulation_time:
            simulation_time = S.namelist.Main.simulation_time
            number_of_timesteps = int(simulation_time / timestep)
            print("simulation_time:", simulation_time)
        print("number_of_timesteps:", number_of_timesteps)
        print("geometry:", S.namelist.Main.geometry)  # print the simulation dimensions
        print("number_of_cells:", S.namelist.Main.number_of_cells)
        print("cell_length:", S.namelist.Main.cell_length)
        print("grid_length:", S.namelist.Main.grid_length)
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

    def plot_scalar(self):
        """Plot Scalar field"""
        import happi

        S = self.SmileiSimulation
        Utot = S.Scalar("Utot", label="Total energy")
        Ukin = S.Scalar("Ukin", label="Total kinetic energy")
        Uelm = S.Scalar("Uelm", label="Total electromagnetic energy")
        Ubal = S.Scalar("Ubal", label="Balance energy")
        Ukin_ion = S.Scalar("Ukin_solar_ion", label="Total solar ion energy")
        Ukin_electron = S.Scalar(
            "Ukin_solar_electron", label="Total solar electron energy"
        )
        happi.multiPlot(
            Utot,
            Ukin,
            Uelm,
            Ubal,
            Ukin_ion,
            Ukin_electron,
            saveAs=self.results_path + "/" + "figures/scalars.png",
        )

    def export_to_vtk(self):
        """Export the results in VTK format"""
        for field in self.fields_smilei:
            self.SmileiSimulation.Field(0, field).toVTK()

    def export_to_yt(self):
        """Save the field file from `Smilei` format to `yt` compatible format."""
        import yt

        timesteps = (
            self.SmileiSimulation.Field(0, self.fields_smilei[0]).getTimesteps().astype(int)
        )
        for timestep in timesteps:
            data = {}
            field_types = {}
            for field_yt, field_smilei in zip(self.fields_yt, self.fields_smilei):
                data[field_yt] = yt.YTArray(
                    self.SmileiSimulation.Field(0, field_smilei).getData(timestep=timestep)[0],
                    "dimensionless",
                )
                field_types[field_yt] = "stream"
            ds ={}
            yt.save_as_dataset(
                ds,
                self.results_path + "/" + "data/yt_data_" + str(timestep) + ".h5",
                data,
                field_types=field_types,
            )

