# poro_shale_threshold = 0.1 #self.idata.rock.poro_shale_threshold  # short name #Old
# poro = np.array(self.reservoir.mesh.poro) #Old
#
# self.reservoir.hcap[poro > poro_shale_threshold] = self.idata.rock.hcap_sand  # Naar kijken #Old
# self.reservoir.hcap[poro <= poro_shale_threshold] = self.idata.rock.hcap_shale #Naar kijken #Old
#
# self.reservoir.conduction[poro > poro_shale_threshold] = self.idata.rock.conduction_sand  # Naar kijken #Old
# self.reservoir.conduction[poro <= poro_shale_threshold] = self.idata.rock.conduction_shale  # Naar kijken #Old

# nb = len(np.array(self.reservoir.mesh.poro, copy=False))
#
# op_num_fixed = self.reservoir.mesh.op_num[:nb]
# self.reservoir.conduction[op_num_fixed == 2] = self.idata.rock.conduction_shale
# self.reservoir.conduction[op_num_fixed == 1] = self.idata.rock.conduction_sand
# self.reservoir.conduction[op_num_fixed == 0] = self.idata.rock.conduction_sand
#
# self.reservoir.hcap[op_num_fixed == 2] = self.idata.rock.hcap_shale
# self.reservoir.hcap[op_num_fixed == 1] = self.idata.rock.hcap_sand
# self.reservoir.hcap[op_num_fixed == 0] = self.idata.rock.hcap_sand



# Get SATNUM values ONLY for active cells #NEW
# Get the SATNUM values stored in the reservoir #NEW
#satnum_full = np.array(self.reservoir.satnum)  # Full SATNUM, including inactive cells #NEW

# Get the mapping of active cells #NEW
#global_to_local = np.array(self.reservoir.discr_mesh.global_to_local, copy=False)  # NEW

# Select only active cell SATNUM values #NEW
#satnum_active = satnum_full[global_to_local >= 0]  # NEW

# mask_shale = (satnum_full == 3) & (global_to_local >= 0)
# mask_sand = ((satnum_full == 1) | (satnum_full == 2)) & (global_to_local >= 0)
#
# self.reservoir.conduction[mask_shale] = self.idata.rock.conduction_shale
# self.reservoir.conduction[mask_sand] = self.idata.rock.conduction_sand
#
# self.reservoir.hcap[mask_shale] = self.idata.rock.hcap_shale
# self.reservoir.hcap[mask_sand] = self.idata.rock.hcap_sand







# # 🏗 Assign rock properties based on SATNUM: #NEW
# self.reservoir.conduction[satnum_active == 3] = self.idata.rock.conduction_shale  # SH = 3 #NEW
# self.reservoir.conduction[satnum_active == 2] = self.idata.rock.conduction_sand
# self.reservoir.conduction[satnum_active == 1] = self.idata.rock.conduction_sand
#
# #self.reservoir.conduction[(satnum_active == 1) | (satnum_active == 2)] = self.idata.rock.conduction_sand  # SS = 1 or 2 #NEW
#
# self.reservoir.hcap[satnum_active == 3] = self.idata.rock.hcap_shale
# self.reservoir.hcap[satnum_active == 2] = self.idata.rock.hcap_sand # Shale heat capacity #NEW
# self.reservoir.hcap[satnum_active == 1] = self.idata.rock.hcap_sand
