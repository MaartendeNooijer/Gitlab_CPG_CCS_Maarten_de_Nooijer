import numpy as np
import os
from Darts_Packages_Adjusted.cpg_reservoir import CPG_Reservoir, save_array, read_arrays, check_arrays, make_burden_layers, make_full_cube
#from darts.reservoirs.cpg_reservoir import CPG_Reservoir, save_array, read_arrays, check_arrays, make_burden_layers, make_full_cube
from darts.discretizer import load_single_float_keyword
from darts.engines import value_vector

from darts.tools.gen_cpg_grid import gen_cpg_grid

from darts.models.cicd_model import CICDModel

from darts.engines import value_vector, index_vector

def fmt(x):
    return '{:.3}'.format(x)

#####################################################

class Model_CPG(CICDModel):
    def __init__(self):
        super().__init__()

    def init_reservoir(self):
        if self.idata.generate_grid:
            if self.idata.grid_out_dir is None:
                self.idata.gridname = None
                self.idata.propname = None
            else:  # save generated grid to grdecl files
                os.makedirs(self.idata.grid_out_dir, exist_ok=True)
                self.idata.gridname = os.path.join(self.idata.grid_out_dir, 'grid.grdecl')
                self.idata.propname = os.path.join(self.idata.grid_out_dir, 'reservoir.in')
                arrays = gen_cpg_grid(nx=self.idata.geom.nx, ny=self.idata.geom.ny, nz=self.idata.geom.nz,
                                  dx=self.idata.geom.dx, dy=self.idata.geom.dy, dz=self.idata.geom.dz,
                                  start_z=self.idata.geom.start_z,
                                  permx=self.idata.rock.permx, permy=self.idata.rock.permy, permz=self.idata.rock.permz,
                                  poro=self.idata.rock.poro,
                                  gridname=self.idata.gridname, propname=self.idata.propname)
        else:
            # read grid and rock properties
            arrays = read_arrays(self.idata.gridfile, self.idata.propfile)
            check_arrays(arrays)
            # if self.physics_type == 'deadoil':  # set inactive cells with small porosity (isothermal case)
            #     arrays['ACTNUM'][arrays['PORO'] < self.idata.geom.min_poro] = 0
            # elif self.physics_type == 'geothermal':  # process cells with small poro (thermal case)
            #     for arr in ['PORO', 'PERMX', 'PERMY', 'PERMZ']:
            #         arrays[arr][arrays['PORO'] < self.idata.geom.min_poro] = self.idata.geom.min_poro

        if self.idata.geom.burden_layers > 0:
            # add over- and underburden layers
            make_burden_layers(number_of_burden_layers=self.idata.geom.burden_layers,
                               initial_thickness=self.idata.geom.burden_init_thickness,
                               property_dictionary=arrays,
                               burden_layer_prop_value=self.idata.rock.burden_prop)

        self.reservoir = CPG_Reservoir(self.timer, arrays, faultfile=self.idata.faultfile, minpv=self.idata.geom.minpv) #NEW #was self.reservoir = CPG_Reservoir(self.timer, arrays, minpv=self.idata.geom.minpv)
        self.reservoir.discretize()


        # store modified arrrays (with burden layers) for output to grdecl
        self.reservoir.input_arrays = arrays

        volume = np.array(self.reservoir.mesh.volume, copy=False)
        poro = np.array(self.reservoir.mesh.poro, copy=False)
        print("Pore volume = " + str(sum(volume[:self.reservoir.mesh.n_blocks] * poro)))

        # imitate open-boundaries with a large volume
        bv = self.idata.geom.bound_volume   # volume, will be assigned to each boundary cell [m3]
        self.reservoir.set_boundary_volume(xz_minus=bv, xz_plus=bv, yz_minus=bv, yz_plus=bv)
        self.reservoir.apply_volume_depth()


        # Get the SATNUM values stored in the reservoir #NEW
        satnum_full = np.array(self.reservoir.satnum)  # Full SATNUM, including inactive cells #NEW

        # Get the mapping of active cells #NEW
        global_to_local = np.array(self.reservoir.discr_mesh.global_to_local, copy=False)  # NEW

        # Select only active cell SATNUM values #NEW
        # satnum_active = satnum_full[global_to_local >= 0]  # NEW

        mask_shale = (satnum_full == 3) & (global_to_local >= 0)
        mask_sand = ((satnum_full == 1) | (satnum_full == 2)) & (global_to_local >= 0)

        self.reservoir.conduction[mask_shale] = self.idata.rock.conduction_shale
        self.reservoir.conduction[mask_sand] = self.idata.rock.conduction_sand

        self.reservoir.hcap[mask_shale] = self.idata.rock.hcap_shale
        self.reservoir.hcap[mask_sand] = self.idata.rock.hcap_sand


        # add hcap and rcond to be saved into mesh.vtu
        l2g = np.array(self.reservoir.discr_mesh.local_to_global, copy=False)
        g2l = np.array(self.reservoir.discr_mesh.global_to_local, copy=False)

        self.reservoir.global_data.update({'heat_capacity': make_full_cube(self.reservoir.hcap, l2g, g2l),
                                           'rock_conduction': make_full_cube(self.reservoir.conduction, l2g, g2l)})
                                           # 'rock_conduction': make_full_cube(self.reservoir.conduction, l2g, g2l) })

        self.set_physics()

        # time stepping and convergence parameters
        sim = self.idata.sim  # short name
        self.set_sim_params(first_ts=sim.first_ts, mult_ts=sim.mult_ts, max_ts=sim.max_ts, runtime=sim.runtime,
                            tol_newton=sim.tol_newton, tol_linear=sim.tol_linear)
        if hasattr(sim, 'linear_type'):
            self.params.linear_type = sim.linear_type

        self.timer.node["initialization"].stop()
        self.check_conductivity()


    def set_wells(self):
        # read perforation data from a file
        if hasattr(self.idata, 'schfile'):
            # apply to the reservoir; add wells and perforations, 1-based indices
            for wname, wdata in self.idata.well_data.wells.items():
                self.reservoir.add_well(wname)
                for perf_tuple in wdata.perforations:
                    perf = perf_tuple[1]
                    # adjust to account for added overburden layers
                    perf_ijk_new = (perf.loc_ijk[0], perf.loc_ijk[1], perf.loc_ijk[2] + self.idata.geom.burden_layers)
                    self.reservoir.add_perforation(wname,
                                                   cell_index=perf_ijk_new,
                                                   well_index=perf.well_index, well_indexD=perf.well_indexD,
                                                   multi_segment=perf.multi_segment, verbose=True)
        else:
            # add wells and perforations, 1-based indices
            for wname, wdata in self.idata.well_data.wells.items():
                self.reservoir.add_well(wname)
                for k in range(1 + self.idata.geom.burden_layers,  self.reservoir.nz+1-self.idata.geom.burden_layers):
                    self.reservoir.add_perforation(wname,
                                                   cell_index=(wdata.location.I, wdata.location.J, k),
                                                   well_index=None, multi_segment=False, verbose=True)

    def set_initial_pressure_from_file(self, fname : str):
        # set initial pressure
        p_cpp = value_vector()
        load_single_float_keyword(p_cpp, fname, 'PRESSURE', -1)
        p_file = np.array(p_cpp, copy=False)
        p_mesh = np.array(self.reservoir.mesh.pressure, copy=False)
        try:
            actnum = np.array(self.reservoir.actnum, copy=False) # CPG Reservoir
        except:
            actnum = self.reservoir.global_data['actnum']  #Struct reservoir
        p_mesh[:self.reservoir.mesh.n_res_blocks * 2] = p_file[actnum > 0]

    def well_is_inj(self, wname : str):  # determine well control by its name
        return "INJ" in wname

    def do_after_step(self):
        # save to grdecl file after each time step
        #self.reservoir.save_grdecl(self.get_arrays(), os.path.join(out_dir, 'res_' + str(ti+1)))
        self.physics.engine.report()
        #self.print_well_rate()

    def check_conductivity(self):
        """
        Verifies the conductivity values assigned based on SATNUM.
        Ensures that shale and sandstone have the correct properties.
        """
        # Get SATNUM values again for active cells
        satnum_full = np.array(self.reservoir.satnum)  # Full-grid SATNUM
        global_to_local = np.array(self.reservoir.discr_mesh.global_to_local, copy=False)

        # Use only active cells
        satnum_active = satnum_full[global_to_local >= 0]

        # Get conductivity values
        conduction_values = np.array(self.reservoir.conduction, copy=False)  # Active-cell conductivity

        # Expected values
        expected_shale = self.idata.rock.conduction_shale
        expected_sandstone = self.idata.rock.conduction_sand

        # Count occurrences
        unique_vals, counts = np.unique(conduction_values, return_counts=True)
        shale_count = np.sum(conduction_values == expected_shale)
        sandstone_count = np.sum(conduction_values == expected_sandstone)

        # Print results
        print(f"\n **Conductivity Check Report** ")
        print(f"Unique conductivity values in the reservoir: {dict(zip(unique_vals, counts))}")
        print(f"Shale assignments (SATNUM=3): {shale_count} cells → Expected: {expected_shale} kJ/m/day/K")
        print(
            f"Sandstone assignments (SATNUM=1 or 2): {sandstone_count} cells → Expected: {expected_sandstone} kJ/m/day/K")

        # Debugging: Check if any unexpected values exist
        unexpected_conductivity = [val for val in unique_vals if val not in [expected_shale, expected_sandstone]]
        if unexpected_conductivity:
            print(
                f"Warning: Found unexpected conductivity values: {unexpected_conductivity}. Check SATNUM processing!")

        print("Conductivity verification complete!\n")




