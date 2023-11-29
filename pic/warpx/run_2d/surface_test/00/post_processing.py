# %%
from pywarpx import picmi
import yt
from yt import derived_field
yt.enable_parallelism()

# %%
@derived_field(name="density_e-", sampling_type="local",display_name=r"\rho_e")
def _density_e(field, data):
    return data['boxlib','rho_e-']/ -picmi.constants.q_e

@derived_field(name="density_p+", sampling_type="local",display_name=r"\rho_p")
def _density_p(field, data):
    return data['boxlib','rho_p+']/picmi.constants.q_e

@derived_field(name="density", sampling_type="local",display_name=r"\rho")
def _density(field, data):
    return data['boxlib','rho_p+']/picmi.constants.q_e + data['boxlib','rho_e-']/picmi.constants.q_e

# %%
diag_name = 'plt'
plotfile = './diags/plotfiles/{}??????'.format(diag_name)
ts = yt.load(plotfile)

from pathlib import Path

subdirectories = ["figures"]
for path in subdirectories:
    Path(path).mkdir(parents=True, exist_ok=True)

# %%
_check = 1
if _check:
    from icecream import ic
    ds = ts[0]
    ic(ds.time_unit)
    ic(ds.current_time) 
    # ic(ds.field_list)
    # ic(ds.derived_field_list)
    

# %%
def _SlicePlot(fields,name, _Multipanel = False, nrows_ncols = None, ):
    for ds in ts.piter():
        slc = yt.SlicePlot(ds, "z", fields,origin="native") # Create a sliceplot object
        if _Multipanel:
            fig = slc.export_to_mpl_figure(nrows_ncols)
            fig.tight_layout()
            fig.savefig('figures/{}_{}.png'.format(str(ds), name))
        else:
            slc.save('figures/')

fields = [('boxlib','Bx'),('boxlib','By'),('boxlib','Bz'),('mesh', 'magnetic_field_strength'),]
_SlicePlot(fields,'B', _Multipanel = True, nrows_ncols = (2,2))

fields = [('boxlib','Ex'),('boxlib','Ey'),('boxlib','Ez'),]
_SlicePlot(fields,'E', _Multipanel = True, nrows_ncols = (3,1))


# %%

def _rhoSlicePlot(fields, directions, name , _Multipanel = False, nrows_ncols = None, ):
    for ds in ts.piter():
        for direction in directions:
            slc = yt.SlicePlot(ds, direction, fields ,origin="native")

            if _Multipanel:
                fig = slc.export_to_mpl_figure(nrows_ncols)
                fig.tight_layout()
                fig.savefig('figures/{}_Slice_{}_{}.png'.format(str(ds),direction,name))
            else:
                slc.save('figures/')

# directions = ["x", "y", "z"]
directions = ["z"]
fields = [('boxlib','rho_e-'),('boxlib','rho_p+'), ('boxlib','rho')]
_rhoSlicePlot(fields,directions,'rho',_Multipanel = True, nrows_ncols = (3,1))

fields = [('gas','density_e-'),('gas','density_p+'), ('gas','density')]
_rhoSlicePlot(fields,directions,'density',_Multipanel = True, nrows_ncols = (3,1))

# %%
def _ParticlePhasePlot():
    for ds in ts.piter():
        p=  yt.ParticlePhasePlot( ds.all_data(), ('e-','particle_position_x'), ('e-','particle_position_y'), ('e-','particle_weight'))
        p.save('figures/')

        p=  yt.ParticlePhasePlot( ds.all_data(), ('e-','particle_position_x'), ('e-','particle_position_z'), ('e-','particle_weight'))
        p.save('figures/')

# _ParticlePhasePlot()