'''Generate the HDF5 file for LMA field

>>> from icecream import ic
>>> Bx_py, By_py, Bz_py = LMA_field_plasmapy()
>>> Bx_np, By_np, Bz_np = LMA_field_np()
>>> ic(abs(Bx_np.value).min())
>>> ic(((Bx_np.value-Bx_py.value)/Bx_np.value).min())'
>>> np.divide(Bx_np, Bx_py, out=np.ones_like(Bx_np),  where=Bx_py!=0).max()
>>> Bx_np.units==Bx_py.units'
'''
from parameters import *

def domain(grid_length, cell_length, number_of_cells, length_unit=1):
    return [np.linspace(0.-cell_length[i]*2, grid_length[i]+cell_length[i]*3, number_of_cells[i]+6)*length_unit for i in range(3)]

def LMA_field_np():
    '''Generate the HDF5 file for LMA field
    '''
    from numpy import pi

    r_0_x = s_d; r_0_y = L_y / 2;r_0_z = L_z / 2

    m_d_x, m_d_y, m_d_z = m_d

    x, y, z = domain(grid_length_norm, cell_length_norm, number_of_cells)
    r_x, r_y, r_z = np.meshgrid(x, y, z, indexing='ij')*L_r
    
    r = np.sqrt((r_x-r_0_x)**2+(r_y- r_0_y)**2+(r_z-r_0_z)**2)
    mr_dot = (r_x - r_0_x)*m_d_x + (r_y - r_0_y)*m_d_y +  (r_z - r_0_z)*m_d_z
    Bx = mu0 / (4 * pi)*(3*(r_x-r_0_x)*(mr_dot)/r**5 - m_d_x/r**3)/B_r
    By = mu0 / (4 * pi)*(3*(r_y-r_0_y)*(mr_dot)/r**5 - m_d_y/r**3)/B_r
    Bz = mu0 / (4 * pi)*(3*(r_z-r_0_z)*(mr_dot)/r**5 - m_d_z/r**3)/B_r
    return Bx, By, Bz

def LMA_field_plasmapy():
    '''Generate the HDF5 file for LMA field
    '''
    from unyt import unyt_array
    from plasmapy.formulary import magnetostatics
    from plasmapy.plasma.sources import Plasma3D
    r_0 = unyt_array([s_d, L_y / 2, L_z / 2])

    dipole = magnetostatics.MagneticDipole(
        m_d.to_astropy() , r_0.to_astropy()
    )

    domain_x, domain_y, domain_z = domain(grid_length_norm, cell_length_norm, number_of_cells, L_r)

    plasma = Plasma3D(
        domain_x=domain_x.to_astropy(), domain_y=domain_y.to_astropy(), domain_z=domain_z.to_astropy()
    )
    plasma.add_magnetostatic(dipole)


    Bx = unyt_array.from_astropy(plasma.magnetic_field[0, :, :, :])/B_r
    By = unyt_array.from_astropy(plasma.magnetic_field[1, :, :, :])/B_r
    Bz = unyt_array.from_astropy(plasma.magnetic_field[2, :, :, :])/B_r

    return Bx, By, Bz


def plot_LMA_mag_field_yt(results_path = ".", axis="z"):
    ''' Plot the initial LMA magnetic field
    >>> results_path = "run/08/"
    >>> plot_LMA_mag_field_yt(results_path)
    '''
    import h5py
    data = {}
    bbox = [
    [0, grid_length_norm[0]],
    [-grid_length_norm[1] / 2, grid_length_norm[1] / 2],
    [-grid_length_norm[2] / 2, grid_length_norm[2] / 2],
    ] * L_r
    # bbox = [
    # [0-cell_length[0]*2, grid_length[0]+cell_length[0]*2],
    # [-grid_length[1] / 2-cell_length[1]*2, grid_length[1] / 2+cell_length[1]*2],
    # [-grid_length[2] / 2-cell_length[2]*2, grid_length[2] / 2+cell_length[2]*2] 
    # ] * L_r
    with h5py.File(results_path+'Fields_input.h5', 'r') as hf:
        for field_yt, field_smilei in zip(B_fields_yt, B_fields_smilei):
            data[field_yt] = hf[field_smilei][:]
    domain_dimensions = data[B_fields_yt[0]].shape
    ic(domain_dimensions)
    data = {k: (v, u) for (k, v), u in zip(data.items(), B_units)}
    ds = yt.load_uniform_grid(
        data,
        domain_dimensions,
        length_unit="m",
        bbox=bbox,
        nprocs=24,
        sim_time=0,
        periodicity=(False, False, False),
        unit_system="mks",
    )
    ds.unit_registry.add("d_i", float(d_i), length, tex_repr="d_i")
    slc = yt.SlicePlot(
        ds,
        axis,
        [("stream", "magnetic_field_strength")]
        + [("stream", field) for field in B_fields_yt],
        origin="native",
    )
    slc.set_axes_unit("d_i")
    slc.set_colorbar_label(("stream", "magnetic_field_strength"), r"$B_m$")
    for field, label in zip(B_fields_yt, B_labels):
        slc.set_colorbar_label(("stream", field), label)
        slc.set_log(("stream", field), False)
    fig = slc.export_to_mpl_figure((2, 2))
    fig.tight_layout()
    fig.savefig(results_path + "figures/B/yt_LMA_B_")
    return fig

def main():
    import h5py
    Bx, By, Bz = LMA_field_np()
    with h5py.File(LMA_filename, 'w') as hf:
        hf.create_dataset("Bx",  data=Bx)
        hf.create_dataset("By",  data=By)
        hf.create_dataset("Bz",  data=Bz)
        print(LMA_filename + " is created.")

if __name__ == '__main__':
    main()
    # with h5py.File('Fields_input.h5', 'r') as hf:
    #     data = hf['B_x'][:]