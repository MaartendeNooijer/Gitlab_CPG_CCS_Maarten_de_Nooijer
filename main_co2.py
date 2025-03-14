import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time
import os, sys
from darts.engines import value_vector, index_vector
from darts.engines import redirect_darts_output, value_vector
from darts.tools.plot_darts import *
from darts.tools.logging import redirect_all_output, abort_redirection
from sympy.physics.units import pressure

from darts.physics.super.property_container import PropertyContainer

#from model_geothermal import ModelGeothermal
# from model_deadoil import ModelDeadOil
from model_co2 import ModelCCS

from darts.engines import redirect_darts_output

def run(physics_type : str, case: str, out_dir: str, export_vtk=True, redirect_log=True, platform='cpu'):
    '''
    :param physics_type: "geothermal" or "dead_oil"
    :param case: input grid name
    :param out_dir: directory name for outpult files
    :param export_vtk:
    :return:
    '''
    print('Test started', 'physics_type:', physics_type, 'case:', case, 'platform=', platform)
    # base_results_dir = "results"
    # out_dir = os.path.join(base_results_dir, out_dir)
    os.makedirs(out_dir, exist_ok=True)
    log_filename = os.path.join(out_dir, 'run.log')
    if redirect_log:
        log_stream = redirect_all_output(log_filename)



    if physics_type == 'geothermal':
        m = ModelGeothermal(iapws_physics=True)
    elif physics_type == 'deadoil':
        m = ModelDeadOil()
    elif physics_type == 'ccs':
        m = ModelCCS()
    else:
        print('Error: wrong physics specified:', physics_type)
        exit(1)

    # physics_type = 'ccs'
    #
    # m = ModelCCS()
    #
    m.physics_type = physics_type

    m.set_input_data(case=case)

    m.init_reservoir()

    # init_reservoir.reservoir.mesh.op_num

    m.init(output_folder=out_dir, platform=platform)
    #m.reservoir.mesh.init_grav_coef(0)
    m.save_data_to_h5(kind = 'solution')
    m.set_well_controls()

    # Retrieve SATNUM (facies ID) from the reservoir
    satnum_array = np.array(m.reservoir.satnum, copy=False) - 1 #
    m.reservoir.mesh.op_num = index_vector([int(x) for x in satnum_array] + [0, 0])


    # m.reservoir.save_grdecl(m.get_arrays(ith_step=0), os.path.join(out_dir, 'res_init')) #NEW not sure if it works for CCS #doesn't work for test case

    ret = m.run_simulation()
    if ret != 0:
        exit(1)

    m.reservoir.centers_to_vtk(out_dir)

    #m.reservoir.save_grdecl(m.get_arrays(ith_step=1), os.path.join(out_dir, 'res_last'))
    
    m.print_timers()
    #m.print_stat()

    if export_vtk:
        # read h5 file and write vtk
        m.reservoir.create_vtk_wells(output_directory=out_dir)
        for ith_step in range(len(m.idata.sim.time_steps)):
            m.output_to_vtk(ith_step=ith_step)

    def add_columns_time_data(time_data):
        molar_mass_co2 = 44.01 #kg/kmol
        time_data['Time (years)'] = time_data['time'] / 365.25
        for k in time_data.keys():
            if physics_type == 'ccs' and 'V rate (m3/day)' in k:
                time_data[k.replace('V rate (m3/day)', 'V rate (kmol/day)')] = time_data[k]
                time_data[k.replace('V rate (m3/day)', 'V rate (ton/day)')] = time_data[k] * molar_mass_co2 / 1000 #tons
                time_data.drop(columns=k, inplace=True)

            if physics_type == 'ccs' and 'V  volume (m3)' in k:
                time_data[k.replace('V  volume (m3)', 'V volume (kmol)')] = time_data[k]
                time_data[k.replace('V  volume (m3)', 'V volume (Mt/year)')] = time_data[k] * molar_mass_co2 / 1000 / 1000000 #Mt per year
                time_data.drop(columns=k, inplace=True)

    # from add_colums_output import add_columns_time_data, add_columns_time_data_report

    time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)
    add_columns_time_data(time_data)
    time_data.to_pickle(os.path.join(out_dir, 'time_data.pkl'))

    time_data_report = pd.DataFrame.from_dict(m.physics.engine.time_data_report)
    add_columns_time_data(time_data_report)
    time_data_report.to_pickle(os.path.join(out_dir, 'time_data_report.pkl'))

    # time_data = pd.DataFrame.from_dict(m.physics.engine.time_data)
    # add_columns_time_data(time_data)
    # time_data.to_pickle(os.path.join(out_dir, 'time_data.pkl'))
    #
    # time_data_report = pd.DataFrame.from_dict(m.physics.engine.time_data_report)
    # add_columns_time_data(time_data_report)
    # time_data_report.to_pickle(os.path.join(out_dir, 'time_data_report.pkl'))

    writer = pd.ExcelWriter(os.path.join(out_dir, 'time_data.xlsx'))
    time_data.to_excel(writer, sheet_name='time_data')
    writer.close()

    # filter time_data_report and write to xlsx
    # list the column names that should be removed
    press_gridcells = time_data_report.filter(like='reservoir').columns.tolist()
    chem_cols = time_data_report.filter(like='Kmol').columns.tolist()
    #remove columns from data
    time_data_report.drop(columns=press_gridcells + chem_cols, inplace=True)
    # add time in years
    time_data_report['Time (years)'] = time_data_report['time'] / 365.25
    writer = pd.ExcelWriter(os.path.join(out_dir, 'time_data_report.xlsx'))
    time_data_report.to_excel(writer, sheet_name='time_data_report')
    writer.close()

    #failed, sim_time = check_performance_local(m=m, case=case, physics_type=physics_type) #NEW, was activated

    # if redirect_log:
    #     abort_redirection(log_stream)
    # print('Failed' if failed else 'Ok')

    return  time_data, time_data_report, m.idata.well_data.wells.keys(), m.well_is_inj #NEW, originalaly also failed, sim_time included

##########################################################################################################
def plot_results(wells, well_is_inj, time_data_list, time_data_report_list, label_list, physics_type, out_dir):
    plt.rc('font', size=12)

    for well_name in wells:
        # if well_is_inj(well_name):
        #     continue

        # ax = None
        # for time_data_report, label in zip(time_data_report_list, label_list):
        #     print(time_data_report)
        #     ax = plot_gas_rate_darts_volume(well_name, time_data_report, ax=ax)  # , label=label) #NEW was ax = plot_gas_rate_darts(well_name, time_data_report, ax=ax)
        # ax.set(xlabel="Days", ylabel="Gas Rate [m3/day]")
        # plt.tight_layout()
        # plt.savefig(os.path.join(out_dir, 'well_gas_rate_' + well_name + '_' + case + '.png'))
        # plt.close()

        # use time_data here as we are going to compute a cumulative plot
        ax = None
        for time_data, label in zip(time_data_list, label_list):
            ax = plot_total_inj_gas_rate_darts_volume(time_data, ax=ax)  # , label=label) #NEW was ax = plot_total_inj_gas_rate_darts(time_data, ax=ax)
        ax.set(xlabel="Days", ylabel="Total Inj Gas Rate [Ton/Day]")
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, 'total_inj_gas_rate_' + well_name + '_' + case + '.png'))
        plt.close()

    for well_name in wells:
        ax = None
        for time_data_report, label in zip(time_data_list, label_list): #NEW was time_data_report_list to check per year
            ax = plot_bhp_darts(well_name, time_data_report, ax=ax)#, label=label)
        ax.set(xlabel="Days", ylabel="BHP [bar]")
        plt.savefig(os.path.join(out_dir, 'well_' + well_name + '_bhp_' + case + '.png'))
        plt.tight_layout()
        plt.close()

##########################################################################################################
# for CI/CD
def check_performance_local(m, case, physics_type):
    import platform

    os.makedirs('ref', exist_ok=True)

    pkl_suffix = ''
    if os.getenv('TEST_GPU') != None and os.getenv('TEST_GPU') == '1':
        pkl_suffix = '_gpu'
    elif os.getenv('ODLS') != None and os.getenv('ODLS') == '-a':
        pkl_suffix = '_iter'
    else:
        pkl_suffix = '_odls'
    print('pkl_suffix=', pkl_suffix)

    file_name = os.path.join('ref', 'perf_' + platform.system().lower()[:3] + pkl_suffix +
                             '_' + case + '_' + physics_type + '.pkl')
    overwrite = 0
    if os.getenv('UPLOAD_PKL') == '1':
        overwrite = 1

    is_plk_exist = os.path.isfile(file_name)

    failed = m.check_performance(perf_file=file_name, overwrite=overwrite, pkl_suffix=pkl_suffix)

    if not is_plk_exist or overwrite == '1':
        m.save_performance_data(file_name=file_name, pkl_suffix=pkl_suffix)
        return False, 0.0

    if is_plk_exist:
        return (failed > 0), -1.0 #data[-1]['simulation time']
    else:
        return False, -1.0

def run_test(args: list = [], platform='cpu'):
    if len(args) > 1:
        case = args[0]
        physics_type = args[1]

        out_dir = 'results_' + physics_type + '_' + case
        ret = run(case=case, physics_type=physics_type, out_dir=out_dir, platform=platform)
        return ret[0], ret[1] #failed_flag, sim_time
    else:
        print('Not enough arguments provided')
        return True, 0.0
##########################################################################################################

if __name__ == '__main__':
    platform = 'cpu'
    if os.getenv('TEST_GPU') != None and os.getenv('TEST_GPU') == '1':
            platform = 'gpu'

    physics_list = []
    #physics_list += ['geothermal']
    #physics_list += ['deadoil']
    physics_list += ['ccs']

    cases_list = []
    #cases_list += ['generate_5x3x4']
    #cases_list += ['generate_51x51x1']
    #cases_list += ['generate_100x100x100']
    #cases_list += ['case_40x40x10']
    cases_list += ['20x20x10']
    #cases_list += ['brugge']
    #cases_list += ['case_test_own_FM3_CO1_G1_TS1_PIX']
    #cases_list += ['case_test_own_FM3_CO1_G1_TS1_PIX_2']
    #cases_list += ['case_test_own_FM3_CO1_G1_TS1_PIX_3']
    #cases_list += ["G1_TS1_FM1_20x20x2_5"]
    #cases_list += ["G1_TS1_FM1_25x25x2_5"]
    # cases_list += ["G1_TS1_FM1_50x50x5"]
    #cases_list += ["G1_TS2_FM2_20x20x2_5"]
    #cases_list += ["G1_TS2_FM2_25x25x2_5"]
    # cases_list += ["G1_TS2_FM2_50x50x5"]
    #cases_list += ["G1_TS3_FM3_20x20x2_5"]
    #cases_list += ["G1_TS3_FM3_25x25x2_5"]
    # cases_list += ["G1_TS3_FM3_50x50x5"]

    well_controls = []
    well_controls += ['wrate']
    #well_controls += ['wbhp']
    #well_controls += ['wperiodic']

    for physics_type in physics_list:
        for case_geom in cases_list:
            for wctrl in well_controls:
                if physics_type == 'deadoil' and wctrl == 'wrate':
                    continue
                case = case_geom + '_' + wctrl
                out_dir = 'results_' + physics_type + '_' + case
                base_results_dir = "results"
                out_dir = os.path.join(base_results_dir, out_dir)
                time_data, time_data_report, wells, well_is_inj = run(physics_type=physics_type, #NEW, originally failed, sim_time, were included
                                                                                        case=case, out_dir=out_dir,
                                                                                        redirect_log=False,
                                                                                        export_vtk= True,
                                                                                        platform=platform)

                # one can read well results from pkl file to add/change well plots without re-running the model
                pkl1_dir = '.'
                pkl_fname = 'time_data.pkl'
                pkl_report_fname = 'time_data_report.pkl'
                time_data_list = [time_data]
                time_data_report_list = [time_data_report]
                label_list = [None]

                # compare the current results with another run
                #pkl1_dir = r'../../../open-darts_dev/models/cpg_sloping_fault/results_' + physics_type + '_' + case_geom
                #time_data_1 = pd.read_pickle(os.path.join(pkl1_dir, pkl_fname))
                #time_data_report_1 = pd.read_pickle(os.path.join(pkl1_dir, pkl_report_fname))
                #time_data_list = [time_data_1, time_data]
                #time_data_report_list = [time_data_report_1, time_data_report]
                #label_list = ['1', 'current']

                plot_results(wells=wells, well_is_inj=well_is_inj,
                             time_data_list=time_data_list, time_data_report_list=time_data_report_list, label_list=label_list,
                             physics_type=physics_type, out_dir=out_dir)
