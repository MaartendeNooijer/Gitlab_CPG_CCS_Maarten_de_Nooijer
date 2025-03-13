import numpy as np

from darts.input.input_data import InputData
from darts.models.darts_model import DartsModel
from set_case import set_input_data
from darts.engines import value_vector

from darts.physics.super.physics import Compositional
from darts.physics.super.property_container import PropertyContainer

from darts.physics.properties.basic import PhaseRelPerm, ConstFunc
from darts.physics.properties.density import Garcia2001
from darts.physics.properties.viscosity import Fenghour1998, Islam2012
from darts.physics.properties.eos_properties import EoSDensity, EoSEnthalpy

from dartsflash.libflash import NegativeFlash
from dartsflash.libflash import CubicEoS, AQEoS, FlashParams, InitialGuess
from dartsflash.components import CompData
from model_cpg import Model_CPG, fmt

from dataclasses import dataclass, field

from scipy.special import erf

# region Dataclasses
@dataclass
class Corey:
    nw: float
    ng: float
    swc: float
    sgc: float
    krwe: float
    krge: float
    labda: float
    p_entry: float
    pcmax: float
    c2: float
    def modify(self, std, mult):
        i = 0
        for attr, value in self.__dict__.items():
            if attr != 'type':
                setattr(self, attr, value * (1 + mult[i] * float(getattr(std, attr))))
            i += 1

    def random(self, std):
        for attr, value in self.__dict__.items():
            if attr != 'type':
                std_in = value * float(getattr(std, attr))
                param = np.random.normal(value, std_in)
                if param < 0:
                    param = 0
                setattr(self, attr, param)



class ModelCCS(Model_CPG):
    def __init__(self):
        self.zero = 1e-10
        super().__init__()

    # def set_physics(self): #ORIGINAL! #NEW
    #     """Physical properties"""
    #     # Fluid components, ions and solid
    #     components = ["H2O", "CO2"]
    #     phases = ["Aq", "V"]
    #     #nc = len(components)
    #     comp_data = CompData(components, setprops=True)
    #
    #     pr = CubicEoS(comp_data, CubicEoS.PR)
    #     # aq = Jager2003(comp_data)
    #     aq = AQEoS(comp_data, AQEoS.Ziabakhsh2012)
    #
    #     flash_params = FlashParams(comp_data)
    #
    #     # EoS-related parameters
    #     flash_params.add_eos("PR", pr)
    #     flash_params.add_eos("AQ", aq)
    #     flash_params.eos_order = ["AQ", "PR"]
    #
    #
    #     # Flash-related parameters
    #     # flash_params.split_switch_tol = 1e-3
    #     temperature = None
    #     # if temperature is None:  # if None, then thermal=True
    #     #     thermal = True
    #     # else:
    #     #     thermal = False
    #
    #     """ properties correlations """
    #     property_container = PropertyContainer(phases_name=phases, components_name=components, Mw=comp_data.Mw,
    #                                            temperature=temperature, min_z=self.zero/10)
    #
    #     property_container.flash_ev = NegativeFlash(flash_params, ["AQ", "PR"], [InitialGuess.Henry_AV])
    #
    #     property_container.density_ev = dict([('V', EoSDensity(pr, comp_data.Mw)),
    #                                           ('Aq', Garcia2001(components))])
    #
    #     property_container.viscosity_ev = dict([('V', Fenghour1998()),
    #                                             ('Aq', Islam2012(components))])
    #
    #     property_container.enthalpy_ev = dict([('V', EoSEnthalpy(pr)),
    #                                            ('Aq', EoSEnthalpy(aq))])
    #     property_container.conductivity_ev = dict([('V', ConstFunc(10.)),
    #                                                ('Aq', ConstFunc(180.)), ])
    #
    #     property_container.rel_perm_ev = dict([('V', PhaseRelPerm("gas")), #Original Activated, #NEW
    #                                             ('Aq', PhaseRelPerm("oil"))])  # This is correct, shouldn't it be water?
    #
    #     property_container.output_props = {"satA": lambda: property_container.sat[0],
    #                                        "satV": lambda: property_container.sat[1], #Question, shouldn't this be flipped around initially?
    #                                        "xCO2": lambda: property_container.x[0, 1],
    #                                        "yH2O": lambda: property_container.x[1, 0]}
    #
    #     self.physics = Compositional(components, phases, self.timer, n_points=self.idata.obl.n_points,
    #                                  min_p=self.idata.obl.min_p, max_p=self.idata.obl.max_p,
    #                                  min_z=self.idata.obl.min_z, max_z=self.idata.obl.max_z,
    #                                  min_t=self.idata.obl.min_t, max_t=self.idata.obl.max_t,
    #                                  thermal=self.idata.obl.thermal, cache=self.idata.obl.cache)
    #
    #     self.physics.add_property_region(property_container)
    #
    #     return

    def set_physics(self): #NEW
        """Assign physical properties, including relative permeability, based on facies (SATNUM)."""

        # Define components and phases
        components = ["H2O", "CO2"]
        phases = ["Aq", "V"]
        comp_data = CompData(components, setprops=True)

        # Define Equations of State (EOS) for CO2 and H2O
        pr = CubicEoS(comp_data, CubicEoS.PR)
        aq = AQEoS(comp_data, AQEoS.Ziabakhsh2012)
        flash_params = FlashParams(comp_data)
        flash_params.add_eos("PR", pr)
        flash_params.add_eos("AQ", aq)
        flash_params.eos_order = ["AQ", "PR"]

        # Retrieve SATNUM (facies ID) from the reservoir
        satnum_array = np.array(self.reservoir.satnum, copy=False)#
        #m.reservoir.mesh.op_num[:16000] = satnum_array # Get SATNUM array

        # # Define relative permeability models for each facies (SATNUM values)
        # facies_rel_perm = {
        #     1: Corey(nw=2.0, ng=2.0, swc=0.15, sgc=0.05, krwe=1.0, krge=0.9, labda=2, p_entry=0.1, pcmax=200, c2=1.5),
        #     # Channel Sands
        #     2: Corey(nw=1.8, ng=2.2, swc=0.10, sgc=0.07, krwe=0.9, krge=0.8, labda=2, p_entry=0.2, pcmax=250, c2=1.5),
        #     # Overbank Sands
        #     3: Corey(nw=1.5, ng=2.5, swc=0.20, sgc=0.10, krwe=0.8, krge=0.7, labda=2, p_entry=0.3, pcmax=300, c2=1.5)
        #     # Background Shales
        # }


        facies_rel_perm = { #based on fluidflower paper
                # **Channel Sand (Facies 1)**
                0: Corey(nw=1.5, ng=1.5, swc=0.10, sgc=0.10, krwe=1.0, krge=1.0, labda=2., p_entry=0.025602, pcmax=300, c2=1.5),
                # **Overbank Sand (Facies 2)**
                1: Corey(nw=1.5, ng=1.5, swc=0.12, sgc=0.10, krwe=1.0, krge=1.0, labda=2., p_entry=0.038706, pcmax=300, c2=1.5),
                # **Shale (Facies 3)**
                2: Corey(nw=1.5, ng=1.5, swc=0.32, sgc=0.10, krwe=1.0, krge=1.0, labda=2., p_entry=1.935314, pcmax=300, c2=1.5)}

        # Initialize physics model
        self.physics = Compositional(
            components=components,
            phases=phases,
            timer=self.timer,
            n_points=self.idata.obl.n_points,
            min_p=self.idata.obl.min_p,
            max_p=self.idata.obl.max_p,
            min_z=self.idata.obl.min_z,
            max_z=self.idata.obl.max_z,
            min_t=self.idata.obl.min_t,
            max_t=self.idata.obl.max_t,
            thermal=self.idata.obl.thermal,
            cache=self.idata.obl.cache
        )

        temperature = None

        # Assign properties per facies using `add_property_region`
        #for facies, corey_params in facies_rel_perm.items(): #original
        for i, (facies, corey_params) in enumerate(facies_rel_perm.items()): #NEW

            property_container = PropertyContainer(phases_name=phases, components_name=components,Mw=comp_data.Mw[:2], #NEW, was property_container = PropertyContainer(phases_name=phases, components_name=components,Mw=comp_data.Mw,
                                                   temperature=temperature, min_z=self.zero/10)

            # Assign phase behavior properties
            property_container.flash_ev = NegativeFlash(flash_params, ["AQ", "PR"], [InitialGuess.Henry_AV]) #NEW, was property_container.flash_ev = NegativeFlash(flash_params, ["AQ", "PR"], [InitialGuess.Henry_AV])

            property_container.density_ev = dict([('V', EoSDensity(eos=pr, Mw=comp_data.Mw[:2])),
                                                  ('Aq', Garcia2001(components)), ]) # NEW, was property_container.density_ev = {'V': EoSDensity(pr, comp_data.Mw),'Aq': Garcia2001(components)}


            property_container.viscosity_ev = dict([('V', Fenghour1998()),
                                                    ('Aq', Islam2012(components)), ]) #NEW, was property_container.viscosity_ev = {'V': Fenghour1998(), 'Aq': Islam2012(components)}

            # if thermal:
            #     property_container.enthalpy_ev = dict([('V', EoSEnthalpy(eos=pr)),
            #                                            ('Aq', EoSEnthalpy(eos=aq)), ]) #NEW
            #
            #     property_container.conductivity_ev = dict([('V', ConstFunc(10)),
            #                                                ('Aq', ConstFunc(180.)), ]) #NEW

            #Original:
            property_container.enthalpy_ev = {'V': EoSEnthalpy(pr),'Aq': EoSEnthalpy(aq)}

            property_container.conductivity_ev = {'V': ConstFunc(10.0), 'Aq': ConstFunc(180.0)}


            # Assign relative permeability and capillary pressure

            property_container.rel_perm_ev = dict([('V', ModBrooksCorey(corey_params, 'V')),
                                                   ('Aq', ModBrooksCorey(corey_params, 'Aq'))]) #NEW
            # #Original:
            # property_container.rel_perm_ev = {
            #     'V': ModBrooksCorey(corey_params, 'V'),  # Vapor/Gas phase relative permeability
            #     'Aq': ModBrooksCorey(corey_params, 'Aq')  # Aqueous/Water phase relative permeability
            # }

            property_container.capillary_pressure_ev = ModCapillaryPressure(corey_params)

            self.physics.add_property_region(property_container, i)

            property_container.output_props = {"satA": lambda ii = i: self.physics.property_containers[ii].sat[0],
                                               "satV": lambda ii = i: self.physics.property_containers[ii].sat[1],
                                               "xCO2": lambda ii = i: self.physics.property_containers[ii].x[0, 1],
                                               "yH2O": lambda ii = i: self.physics.property_containers[ii].x[1, 0]} #NEW

            #Original:
            # property_container.output_props = {"satA": lambda: property_container.sat[0],
            #                                    "satV": lambda: property_container.sat[1],
            #                                    "xCO2": lambda: property_container.x[0, 1],
            #                                    "yH2O": lambda: property_container.x[1, 0]}

            # **Correctly add property region using `add_property_region()`**


        return

    def get_arrays(self, ith_step):
        '''
        :return: dictionary of current unknown arrays (p, T)
        '''
        # Find index of properties to output


        ev_props = self.physics.property_operators[next(iter(self.physics.property_operators))].props_name
        # If output_properties is None, all variables and properties from property_operators will be passed
        props_names = list(ev_props)
        props_names = props_names + ['pressure', 'temperature']

        # print(props_names)#NEW

        timesteps, property_array = self.output_properties(output_properties=props_names, timestep=ith_step)

        #print(property_array) #NEW

        return property_array

    def set_input_data(self, case=''):
        self.idata = InputData(type_hydr='thermal', type_mech='none', init_type='uniform')
        set_input_data(self.idata, case)
        self.idata.faultfile =  self.idata.faultfile #NEW

        self.idata.geom.burden_layers = 0

        # well controls
        wdata = self.idata.well_data
        wells = wdata.wells  # short name
        # set default injection composition
        #wdata.inj = value_vector([self.zero])  # injection composition - water
        y2d = 365.25

        rate_kg_year = 0.1e9 #0.1 Mt
        rate_kg_day = rate_kg_year/y2d
        kg_to_mol = 44.01 #kg/kmol
        rate_kmol_day = rate_kg_day/kg_to_mol

        if 'wbhp' in case:
            for w in wells:
                wdata.add_inj_bhp_control(name=w, bhp=250, comp_index=1, temperature=300)  # kmol/day | bars | K
                #wdata.add_prd_rate_control(time=10 * y2d, name=w, rate=0., comp_index=0, bhp_constraint=70)  # STOP WELL
        elif 'wrate' in case:
            for w in wells:
                wdata.add_inj_rate_control(name=w, rate=rate_kmol_day, comp_index=1, bhp_constraint=250, temperature=300) #rate=5e5/8 is approximately 1 Mt per year #was rate = 6e6 # kmol/day | bars | K
                #wdata.add_prd_rate_control(time=10 * y2d, name=w, rate=0., comp_index=0, bhp_constraint=70)  # STOP WELL

        self.idata.obl.n_points = 1000
        self.idata.obl.zero = 1e-11
        self.idata.obl.min_p = 1.
        self.idata.obl.max_p = 400.
        self.idata.obl.min_t = 273.15
        self.idata.obl.max_t = 373.15
        self.idata.obl.min_z = self.idata.obl.zero
        self.idata.obl.max_z = 1 - self.idata.obl.zero
        self.idata.obl.cache = False
        self.idata.obl.thermal = True

    def set_initial_conditions(self):
        self.temperature_initial_ = 273.15 + 76.85  # K
        self.initial_values = {"pressure": 100.,
                            "H2O": 0.99995,
                            "temperature": self.temperature_initial_
                            }
        super().set_initial_conditions()

        mesh = self.reservoir.mesh
        depth = np.array(mesh.depth, copy=True) #- 2000
        # set initial pressure
        pressure_grad = 100
        pressure = np.array(mesh.pressure, copy=False)
        pressure[:] = depth[:pressure.size] / 1000 * pressure_grad + 1
        # set initial temperature
        temperature_grad = 30
        temperature = np.array(mesh.temperature, copy=False)
        temperature[:] = depth[:pressure.size] / 1000 * temperature_grad + 273.15 + 20
        temperature[:] = 350

    def set_well_controls(self, time: float = 0., verbose=True):
        '''
        :param time: simulation time, [days]
        :return:
        '''
        inj_stream_base = [self.zero * 100]
        eps_time = 1e-15
        for w in self.reservoir.wells:
            # find next well control in controls list for different timesteps
            wctrl = None
            for wctrl_t in self.idata.well_data.wells[w.name].controls:
                if np.fabs(wctrl_t[0] - time) < eps_time:  # check time
                    wctrl = wctrl_t[1]
                    break
            if wctrl is None:
                continue
            if wctrl.type == 'inj':  # INJ well
                inj_stream = inj_stream_base
                if self.physics.thermal:
                    inj_stream += [wctrl.inj_bht]
                if wctrl.mode == 'rate': # rate control
                    w.control = self.physics.new_rate_inj(wctrl.rate, inj_stream, wctrl.comp_index)
                    w.constraint = self.physics.new_bhp_inj(wctrl.bhp_constraint, inj_stream)
                elif wctrl.mode == 'bhp': # BHP control
                    w.control = self.physics.new_bhp_inj(wctrl.bhp, inj_stream)
                else:
                    print('Unknown well ctrl.mode', wctrl.mode)
                    exit(1)
            elif wctrl.type == 'prod':  # PROD well
                if wctrl.mode == 'rate': # rate control
                    w.control = self.physics.new_rate_prod(wctrl.rate, wctrl.comp_index)
                    w.constraint = self.physics.new_bhp_prod(wctrl.bhp_constraint)
                elif wctrl.mode == 'bhp': # BHP control
                    w.control = self.physics.new_bhp_prod(wctrl.bhp)
                else:
                    print('Unknown well ctrl.mode', wctrl.mode)
                    exit(1)
            else:
                print('Unknown well ctrl.type', wctrl.type)
                exit(1)
            if verbose:
                print('set_well_controls: time=', time, 'well=', w.name, w.control, w.constraint)

        # check
        for w in self.reservoir.wells:
            assert w.control is not None, 'well control is not initialized for the well ' + w.name
            if verbose and w.constraint is not None and 'rate' in str(type(w.control)):
                print('A constraint for the well ' + w.name + ' is not initialized!')



class ModBrooksCorey: #NEW
    def __init__(self, corey, phase):

        self.phase = phase

        if self.phase == "Aq":
            self.k_rw_e = corey.krwe
            self.swc = corey.swc
            self.sgc = 0
            self.nw = corey.nw
        else:
            self.k_rg_e = corey.krge
            self.sgc = corey.sgc
            self.swc = 0
            self.ng = corey.ng

    def evaluate(self, sat):
        if self.phase == "Aq":
            Se = (sat - self.swc)/(1 - self.swc - self.sgc)
            if Se > 1:
                Se = 1
            elif Se < 0:
                Se = 0
            k_r = self.k_rw_e * Se ** self.nw
        else:
            Se = (sat - self.sgc) / (1 - self.swc - self.sgc)
            if Se > 1:
                Se = 1
            elif Se < 0:
                Se = 0
            k_r = self.k_rg_e * Se ** self.ng

        return k_r


class ModCapillaryPressure: #NEW
    def __init__(self, corey):
        self.swc = corey.swc
        self.p_entry = corey.p_entry
        self.labda = corey.labda
        # self.labda = 3
        self.eps = 1e-10
        self.pcmax = corey.pcmax
        self.c2 = corey.c2

    def evaluate(self, sat):
        sat_w = sat[1]
        # sat_w = sat
        Se = (sat_w - self.swc)/(1 - self.swc)
        if Se < self.eps:
            Se = self.eps
        # pc = self.p_entry * self.eps ** (1/self.labda) * Se ** (-1/self.labda)  # for p_entry to non-wetting phase
        pc_b = self.p_entry * Se ** (-1/self.c2) # basic capillary pressure
        pc = self.pcmax * erf((pc_b * np.sqrt(np.pi)) / (self.pcmax * 2)) # smoothened capillary pressure
        # if Se > 1 - self.eps:
        #     pc = 0

        # pc = self.p_entry
        Pc = np.array([0, pc], dtype=object)  # V, Aq
        return Pc


class BrooksCorey: #NEW
    def __init__(self, wetting: bool):
        self.sat_wr = 0.15
        # self.sat_nwr = 0.1

        self.lambda_w = 4.2
        self.lambda_nw = 3.7

        self.wetting = wetting

    def evaluate(self, sat_w):
        # From Brooks-Corey (1964)
        Se = (sat_w - self.sat_wr)/(1-self.sat_wr)
        if Se > 1:
            Se = 1
        elif Se < 0:
            Se = 0

        if self.wetting:
            k_r = Se**((2+3*self.lambda_w)/self.lambda_w)
        else:
            k_r = (1-Se)**2 * (1-Se**((2+self.lambda_nw)/self.lambda_nw))

        if k_r > 1:
            k_r = 1
        elif k_r < 0:
            k_r = 0

        return k_r

