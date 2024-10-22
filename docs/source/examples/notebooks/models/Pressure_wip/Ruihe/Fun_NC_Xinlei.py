# For GEM-2 NC paper
import csv, random, os
import pybamm as pb;import pandas as pd   ;import numpy as np;import os;
import matplotlib.pyplot as plt;import os;#import imageio;import timeit
from scipy.io import savemat,loadmat;from pybamm import constants,exp,sqrt;
import matplotlib as mpl; 
from multiprocessing import Queue, Process, set_start_method
from queue import Empty
import openpyxl
import traceback
import random;import time, signal

###################################################################
#############    From Patrick from Github        ##################
###################################################################
class TimeoutFunc(object):
    def __init__(self, func, timeout=None, timeout_val=None):
        assert callable(func), 'Positional argument 1 must be a callable method'
        self.func = func
        self.timeout = timeout
        self.timeout_val = timeout_val
    def _queued_f(self, *args, queue, **kwargs):
        Result_List = self.func(*args, **kwargs)
        try:
            queue.put(Result_List)
        except BrokenPipeError as exc:
            pass

    def __call__(self, *args, **kwargs):
        q = Queue(1)
        p = Process(target=self._queued_f, args=args,
            kwargs={**kwargs, 'queue': q})
        p.daemon = True
        p.start()
        try:
            Result_List = q.get(timeout=self.timeout)
        except Empty as exc: 
            # todo: need to incooperate other cases such as pb.expression_tree.exceptions.ModelError,
            #      pb.expression_tree.exceptions.SolverError
            p.terminate()
            Call_ref = RioCallback() 
            # careful! the length of Result_List should be same as what you get from main function!
            Result_List = [         
                self.timeout_val,
                self.timeout_val,
                Call_ref,
                self.timeout_val]
        return Result_List  
###################################################################
#############    New functions from P3 - 221113        ############
###################################################################
# DEFINE my callback:
class RioCallback(pb.callbacks.Callback):
    def __init__(self, logfile=None):
        self.logfile = logfile
        self.success  = True
        if logfile is None:
            # Use pybamm's logger, which prints to command line
            self.logger = pb.logger
        else:
            # Use a custom logger, this will have its own level so set it to the same
            # level as the pybamm logger (users can override this)
            self.logger = pb.get_new_logger(__name__, logfile)
            self.logger.setLevel(pb.logger.level)
    
    def on_experiment_error(self, logs):
        self.success  = False
    def on_experiment_infeasible(self, logs):
        self.success  = False

class Experiment_error_infeasible(ValueError):
    pass


# define to kill too long runs - Jianbo Huang kindly writes this
class TimeoutError(Exception):
	pass
def handle_signal(signal_num, frame):
	raise TimeoutError

def save_rows_to_csv(file_path, rows, header):
    with open(file_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)  # Write parameter names as the header row
        writer.writerows(rows)
def generate_combinations(Bounds, Num_tot):
    lower_bounds = []
    upper_bounds = []
    for bound in Bounds:
        lower_bounds.append(bound[0])
        upper_bounds.append(bound[1])
    combinations = []
    for _ in range(Num_tot):
        combination = []
        for lower, upper in zip(lower_bounds, upper_bounds):
            value = random.uniform(lower, upper)
            combination.append(value)
        combinations.append(combination)
    return combinations

def save_combinations_to_csv(combinations, parameter_names, filename):
    with open(filename, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(parameter_names)  # Write parameter names as the first row
        for combination in combinations:
            writer.writerow(combination)

def Get_Scan_files(
        BasicPath,Target_name,model_options,
        parameter_names,para_short_name,
        Pack,   num,
        rows_per_file,Bundle):
    
    import itertools
    para_dict_Same = {
        "Cycles within RPT":1,
        "RPT temperature":25,
        
        "Para_Set": "OKane2023",
        "Model option":model_options,
        "Current solvent concentration in the reservoir [mol.m-3]":4541.0,
        "Current electrolyte concentration in the reservoir [mol.m-3]":1000,
        "Ratio of Li-ion concentration change in " 
        "electrolyte consider solvent consumption":1.0,
        'EC initial concentration in electrolyte [mol.m-3]':4541.0,
        'Typical EC concentration in electrolyte [mol.m-3]':4541.0, 
        }
    unchange_key2 = list(para_dict_Same.keys())
    unchange_val2 = list(para_dict_Same.values())
    short_pack = [lst for lst in Pack if len(lst) > 1]
    selected_indices = [i for i, lst in enumerate(Pack) if len(lst) > 1]
    shortList_para_short_name = [para_short_name[i] for i in selected_indices]
    shortList_para_short_name.insert(0,"No")
    really_change_val =  [
        list(comb) for comb in itertools.product(*short_pack)]

    change_val =  [
        list(comb) for comb in itertools.product(*Pack)]
    combinations = [[i+1,*elem, *unchange_val2] for i,elem in enumerate(change_val)]
    comb_short   = [[i+1,*elem] for i,elem in enumerate(really_change_val)]
    parameter_names = [*parameter_names,*unchange_key2]
    print("Total cases number is",len(combinations))
    if Bundle:
        # Specify the total number of cases
        total_cases = len(combinations)
        # Specify the number of rows per CSV file, rows_per_file
        # Calculate the number of files needed
        num_files = (total_cases - 1) // rows_per_file + 1
        # Create the target folder
        folder_path = os.path.join(BasicPath, "Get_Random_sets", Target_name)
        os.makedirs(folder_path, exist_ok=True)
        # Write data to each CSV file
        for i in range(num_files):
            file_name = f"Bundle_{i+1}.csv"
            file_path = os.path.join(folder_path, file_name)
            start_row = i * rows_per_file
            end_row = min(start_row + rows_per_file, total_cases)
            rows = combinations[start_row:end_row]
            save_rows_to_csv(file_path, rows, parameter_names)
        filename = BasicPath+f"/Get_Random_sets/{Target_name}/"+f'{Target_name}.csv'
        filename_short = BasicPath+f"/Get_Random_sets/{Target_name}/"+f'{Target_name}_s.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        save_combinations_to_csv(comb_short, shortList_para_short_name, filename_short)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
        print(f"CSV files created in folder '{Target_name}'.")
    else:
        filename = BasicPath+"/Get_Random_sets/"+f'{Target_name}.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
    return len(combinations)

def get_list_from_tuple(d, num):
    if d[1] > d[0]:
        if d[1] > 100 * d[0]:
            result_list = (np.exp(np.linspace(np.log(d[0]), np.log(d[1]), num=num))).tolist()
        else:
            result_list = (np.linspace(d[0], d[1], num=num)).tolist()
    else:
        result_list = []
    return result_list
def Get_Scan_Orth_Latin(
        BasicPath,Target_name,model_options,
        parameter_names,para_short_name,
        Pack, num,
        rows_per_file,Bundle):
    
    import itertools; from pyDOE import lhs
    para_dict_Same = {
        "Cycles within RPT":1,
        "RPT temperature":25,
        #"Mesh list":[5,5,5,60,20],  
        "Para_Set": "OKane2023",
        "Model option":model_options,
        "Current solvent concentration in the reservoir [mol.m-3]":4541.0,
        "Current electrolyte concentration in the reservoir [mol.m-3]":1000,
        "Ratio of Li-ion concentration change in " 
        "electrolyte consider solvent consumption":1.0,
        'EC initial concentration in electrolyte [mol.m-3]':4541.0,
        'Typical EC concentration in electrolyte [mol.m-3]':4541.0, 
        }
    unchange_key2 = list(para_dict_Same.keys())
    unchange_val2 = list(para_dict_Same.values())
    
    Pack_tuple = []; Pack_tuple_index = []
    Pack_list = [];  Pack_list_index  = []
    for i,item in enumerate(Pack):
        if isinstance(item, tuple):
            Pack_tuple.append(item)
            Pack_tuple_index.append(i)
        elif isinstance(item, list):
            Pack_list.append(item)
            Pack_list_index.append(i)
    com_tuple = []; comb_tu_list =[]
    if len(Pack_tuple) > 1:
        for tuple_i in Pack_tuple:
            com_tuple.append( get_list_from_tuple(tuple_i, num) )
        # apply Latin Hypercube:
        #print(com_tuple)
        samples = lhs(len(com_tuple), samples=num)
        for sample in samples:
            combination = []
            for i, candidate_list in enumerate(com_tuple):
                index = int(sample[i] * num)
                combination.append(candidate_list[index])
            comb_tu_list.append(combination)
    else:
        print("error! Pack_tuple must has 2 elements")
    # apply product sampling:
    comb_li_list = [list(comb) for comb in itertools.product(*Pack_list)]
    #print(comb_tu_list)
    #print(comb_li_list)
    Big_Comb = []
    for comb_tu in comb_tu_list:
        for comb_li in comb_li_list:
            big_comb = [0] * (len(comb_tu)+len(comb_li))
            for comb_tu_i,index in zip(comb_tu,Pack_tuple_index):
                big_comb[index] = comb_tu_i
            for comb_li_i,index in zip(comb_li,Pack_list_index):
                big_comb[index] = comb_li_i
            Big_Comb.append(big_comb)
    #print(Big_Comb)
    Big_Comb
    combinations = [[i+1,*elem, *unchange_val2] for i,elem in enumerate(Big_Comb)]

    parameter_names = [*parameter_names,*unchange_key2]
    print("Total cases number is",len(combinations))
    if Bundle:
        # Specify the total number of cases
        total_cases = len(combinations)
        # Specify the number of rows per CSV file, rows_per_file
        # Calculate the number of files needed
        num_files = (total_cases - 1) // rows_per_file + 1
        # Create the target folder
        folder_path = os.path.join(BasicPath, "Get_Random_sets2", Target_name)
        os.makedirs(folder_path, exist_ok=True)
        # Write data to each CSV file
        for i in range(num_files):
            file_name = f"Bundle_{i+1}.csv"
            file_path = os.path.join(folder_path, file_name)
            start_row = i * rows_per_file
            end_row = min(start_row + rows_per_file, total_cases)
            rows = combinations[start_row:end_row]
            save_rows_to_csv(file_path, rows, parameter_names)
        filename = BasicPath+f"/Get_Random_sets2/{Target_name}/"+f'{Target_name}.csv'
        # filename_short = BasicPath+f"/Get_Random_sets2/{Target_name}/"+f'{Target_name}_s.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        # save_combinations_to_csv(comb_short, shortList_para_short_name, filename_short)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
        print(f"CSV files created in folder '{Target_name}'.")
    else:
        filename = BasicPath+"/Get_Random_sets2/"+f'{Target_name}.csv'
        save_combinations_to_csv(combinations, parameter_names, filename)
        print(f"Combinations saved to '{Target_name}.csv' file.") 
    return len(combinations)


# Define a function to calculate concentration change, whether electrolyte being squeezed out or added in
def Cal_new_con_Update(Sol,Para):   # subscript r means the reservoir
    # Note: c_EC_r is the initial EC  concentraiton in the reservoir; 
    #       c_e_r  is the initial Li+ concentraiton in the reservoir;
    #############################################################################################################################  
    ###############################           Step-1 Prepare parameters:        #################################################
    #############################################################################################################################  
    L_p   =Para["Positive electrode thickness [m]"]
    L_n   =Para["Negative electrode thickness [m]"]
    L_s   =Para["Separator thickness [m]"]
    L_y   =Para["Electrode width [m]"]   # Also update to change A_cc and therefore Q
    L_z   =Para["Electrode height [m]"]
    c_EC_r_old=Para["Current solvent concentration in the reservoir [mol.m-3]"]   # add two new parameters for paper 2
    c_e_r_old =Para["Current electrolyte concentration in the reservoir [mol.m-3]"]
    # Initial EC concentration in JR
    c_EC_JR_old =Para["Bulk solvent concentration [mol.m-3]"]  # Be careful when multiplying with pore volume to get amount in mole. Because of electrolyte dry out, that will give more amount than real.   
    # LLI due to electrode,Ratio of EC and lithium is 1:1 -> EC amount consumed is LLINegSEI[-1]
    LLINegSEI = (
        Sol["Loss of lithium to negative SEI [mol]"].entries[-1] ## revised 241022
        - Sol["Loss of lithium to negative SEI [mol]"].entries[0] )## revised 241022
    LLINegSEIcr = (
        Sol["Loss of lithium to negative SEI on cracks [mol]"].entries[-1]  ## revised 241022
        - 
        Sol["Loss of lithium to negative SEI on cracks [mol]"].entries[0]  ## revised 241022
        )

    ## Current version does not have that 241021,changed to nly LLINegLiP

    LLINegLiP = (
        Sol["Loss of lithium to negative lithium plating [mol]"].entries[-1] 
        - Sol["Loss of lithium to negative lithium plating [mol]"].entries[0])
    
    # LLINegDeadLiP = (
    #     Sol["Loss of lithium to dead negative lithium plating [mol]"].entries[-1] 
    #     - Sol["Loss of lithium to dead negative lithium plating [mol]"].entries[0])
    # LLINegLiP = (
    #     Sol["Loss of lithium to negative lithium plating [mol]"].entries[-1] 
    #     - Sol["Loss of lithium to negative lithium plating [mol]"].entries[0])
    cLi_Xavg  = Sol["X-averaged electrolyte concentration [mol.m-3]"].entries[-1] 
    # Pore volume change with time:
    PoreVolNeg_0 = Sol["X-averaged negative electrode porosity"].entries[0]*L_n*L_y*L_z;
    PoreVolSep_0 = Sol["X-averaged separator porosity"].entries[0]*L_s*L_y*L_z;
    PoreVolPos_0 = Sol["X-averaged positive electrode porosity"].entries[0]*L_p*L_y*L_z;
    PoreVolNeg_1 = Sol["X-averaged negative electrode porosity"].entries[-1]*L_n*L_y*L_z;
    PoreVolSep_1 = Sol["X-averaged separator porosity"].entries[-1]*L_s*L_y*L_z;
    PoreVolPos_1 = Sol["X-averaged positive electrode porosity"].entries[-1]*L_p*L_y*L_z;
    #############################################################################################################################  
    ##### Step-2 Determine How much electrolyte is added (from the reservoir to JR) or squeezed out (from JR to reservoir) ######
    #######################   and finish electrolyte mixing     #################################################################  
    Vol_Elely_Tot_old = Para["Current total electrolyte volume in whole cell [m3]"] 
    Vol_Elely_JR_old  = Para["Current total electrolyte volume in jelly roll [m3]"] 
    if Vol_Elely_Tot_old - Vol_Elely_JR_old < 0:
        print('Model error! Electrolyte in JR is larger than in the cell!')
    Vol_Pore_tot_old  = PoreVolNeg_0 + PoreVolSep_0 + PoreVolPos_0    # pore volume at start time of the run
    Vol_Pore_tot_new  = PoreVolNeg_1 + PoreVolSep_1 + PoreVolPos_1    # pore volume at end   time of the run, intrinsic variable 
    Vol_Pore_decrease = Vol_Elely_JR_old  - Vol_Pore_tot_new # WHY Vol_Elely_JR_old not Vol_Pore_tot_old here? Because even the last state of the last solution (n-1) and the first state of the next solution (n) can be slightly different! 
    # EC:lithium:SEI=2:2:1     for SEI=(CH2OCO2Li)2, but because of too many changes are required, change to 2:1:1 for now
    # update 230703: Simon insist we should do: EC:lithium:SEI  =  2:2:1 and change V_SEI accordingly 
    # Because inner and outer SEI partial molar volume is the same, just set one for whole SEI
    VmolSEI   = Para["Outer SEI partial molar volume [m3.mol-1]"] # 9.8e-5,
    VmolLiP   = Para["Lithium metal partial molar volume [m3.mol-1]"] # 1.3e-05
    VmolEC    = Para["EC partial molar volume [m3.mol-1]"]
    #################   KEY EQUATION FOR DRY-OUT MODEL                   #################
    # UPDATE 230525: assume the formation of dead lithium doesn’t consumed EC
    # Mark: either with 2 or not, will decide how fast electrolyte being consumed!
    # Increase that to increase dry-out
    Vol_EC_consumed  =  ( LLINegSEI + LLINegSEIcr   ) * 2 * VmolEC    
   

    
    Vol_Elely_need   = Vol_EC_consumed - Vol_Pore_decrease
    Vol_SEILiP_increase = 0.5*(
        (LLINegSEI+LLINegSEIcr) * VmolSEI 
        #+ LLINegLiP * VmolLiP
        )    #  volume increase due to SEI+total LiP 
    #################   KEY EQUATION FOR DRY-OUT MODEL                   #################

    Test_V = Vol_SEILiP_increase - Vol_Pore_decrease  #  This value should always be zero, but now not, which becomes the source of error!
    Test_V2= (Vol_Pore_tot_old - Vol_Elely_JR_old) / Vol_Elely_JR_old * 100; # Porosity errors due to first time step
    
    # Start from here, there are serveral variables to be determined:
    #   1) Vol_Elely_add, or Vol_Elely_squeezed; depends on conditions. easier for 'squeezed' condition
    #   2) Vol_Elely_Tot_new, should always equals to Vol_Elely_Tot_old -  Vol_EC_consumed;
    #   3) Vol_Elely_JR_new, for 'add' condition: see old code; for 'squeezed' condition, equals to pore volume in JR
    #   4) Ratio_Dryout, for 'add' condition: see old code; for 'squeezed' condition, equals to Vol_Elely_Tot_new/Vol_Pore_tot_new 
    #   5) Ratio_CeEC_JR and Ratio_CeLi_JR: for 'add' condition: see old code; for 'squeezed' condition, equals to 1 (unchanged)
    #   6) c_e_r_new and c_EC_r_new: for 'add' condition: equals to old ones (unchanged); for 'squeezed' condition, need to carefully calculate     
    #   7) Width_new: for 'add' condition: see old code; for 'squeezed' condition, equals to L_y (unchanged)
    Vol_Elely_Tot_new = Vol_Elely_Tot_old - Vol_EC_consumed;

    
    if Vol_Elely_need < 0:
        print('Electrolyte is being squeezed out, check plated lithium (active and dead)')
        Vol_Elely_squeezed = - Vol_Elely_need;   # Make Vol_Elely_squeezed>0 for simplicity 
        Vol_Elely_add = 0.0;
        Vol_Elely_JR_new = Vol_Pore_tot_new; 
        Ratio_Dryout = Vol_Elely_Tot_new / Vol_Elely_JR_new
        Ratio_CeEC_JR = 1.0;    # the concentrations in JR don't need to change
        Ratio_CeLi_JR = 1.0; 
        SqueezedLiMol   =  Vol_Elely_squeezed*cLi_Xavg;
        SqueezedECMol   =  Vol_Elely_squeezed*c_EC_JR_old;
        Vol_Elely_reservoir_old = Vol_Elely_Tot_old - Vol_Elely_JR_old; 
        LiMol_reservoir_new = Vol_Elely_reservoir_old*c_e_r_old  + SqueezedLiMol;
        ECMol_reservoir_new = Vol_Elely_reservoir_old*c_EC_r_old + SqueezedECMol;
        c_e_r_new= LiMol_reservoir_new / (Vol_Elely_reservoir_old+Vol_Elely_squeezed);
        c_EC_r_new= ECMol_reservoir_new / (Vol_Elely_reservoir_old+Vol_Elely_squeezed);
        Width_new = L_y; 
    else:   # this means Vol_Elely_need >= 0, therefore the folliwing script should works for Vol_Elely_need=0 as well!
        # Important: Calculate added electrolyte based on excessive electrolyte, can be: 1) added as required; 2) added some, but less than required; 3) added 0 
        Vol_Elely_squeezed = 0;  
        if Vol_Elely_Tot_old > Vol_Elely_JR_old:                             # This means Vol_Elely_JR_old = Vol_Pore_tot_old ()
            if Vol_Elely_Tot_old-Vol_Elely_JR_old >= Vol_Elely_need:         # 1) added as required
                Vol_Elely_add     = Vol_Elely_need;  
                Vol_Elely_JR_new  = Vol_Pore_tot_new;  # also equals to 'Vol_Pore_tot_new', or Vol_Pore_tot_old - Vol_Pore_decrease
                Ratio_Dryout = 1.0;
            else:                                                            # 2) added some, but less than required;                                                         
                Vol_Elely_add     = Vol_Elely_Tot_old - Vol_Elely_JR_old;   
                Vol_Elely_JR_new  = Vol_Elely_Tot_new;                       # This means Vol_Elely_JR_new <= Vol_Pore_tot_new
                Ratio_Dryout = Vol_Elely_JR_new/Vol_Pore_tot_new;
        else:                                                                # 3) added 0 
            Vol_Elely_add = 0;
            Vol_Elely_JR_new  = Vol_Elely_Tot_new; 
            Ratio_Dryout = Vol_Elely_JR_new/Vol_Pore_tot_new;

        # Next: start mix electrolyte based on previous defined equation
        # Lithium amount in liquid phase, at initial time point
        TotLi_Elely_JR_Old = Sol["Total lithium in electrolyte [mol]"].entries[-1]; # remember this is only in JR
        # Added lithium and EC amount in the added lithium electrolyte: - 
        AddLiMol   =  Vol_Elely_add*c_e_r_old
        AddECMol   =  Vol_Elely_add*c_EC_r_old
        # Total amount of Li and EC in current electrolyte:
        # Remember Li has two parts, initial and added; EC has three parts, initial, consumed and added
        TotLi_Elely_JR_New   = TotLi_Elely_JR_Old + AddLiMol
        TotECMol_JR   = Vol_Elely_JR_old*c_EC_JR_old - LLINegSEI + AddECMol  # EC:lithium:SEI=2:2:1     for SEI=(CH2OCO2Li)2
        Ratio_CeEC_JR  = TotECMol_JR    /   Vol_Elely_JR_new   / c_EC_JR_old
        Ratio_CeLi_JR  = TotLi_Elely_JR_New    /   TotLi_Elely_JR_Old   /  Ratio_Dryout # Mark, change on 21-11-19
        c_e_r_new  = c_e_r_old;
        c_EC_r_new = c_EC_r_old;
        Width_new   = Ratio_Dryout * L_y;
        
    # Collect the above parameter in Data_Pack to shorten the code  
    Data_Pack   = [
        Vol_EC_consumed, 
        Vol_Elely_need, 
        Test_V, 
        Test_V2, 
        Vol_Elely_add, 
        Vol_Elely_Tot_new, 
        Vol_Elely_JR_new, 
        Vol_Pore_tot_new, 
        Vol_Pore_decrease, 
        c_e_r_new, c_EC_r_new,
        Ratio_Dryout, Ratio_CeEC_JR, 
        Ratio_CeLi_JR,
        Width_new, 
        ]  # 16 in total
    #print('Loss of lithium to negative SEI', LLINegSEI, 'mol') 
    #print('Loss of lithium to dead lithium plating', LLINegDeadLiP, 'mol') 
    #############################################################################################################################  
    ###################       Step-4 Update parameters here        ##############################################################
    #############################################################################################################################
    Para.update(   
        {'Bulk solvent concentration [mol.m-3]':  
         c_EC_JR_old * Ratio_CeEC_JR  })
    Para.update(
        {'EC initial concentration in electrolyte [mol.m-3]':
         c_EC_JR_old * Ratio_CeEC_JR },) 
    Para.update(   
        {'Ratio of Li-ion concentration change in electrolyte consider solvent consumption':  
        Ratio_CeLi_JR }, check_already_exists=False) 
    Para.update(   
        {'Current total electrolyte volume in whole cell [m3]':  Vol_Elely_Tot_new  }, 
        check_already_exists=False)
    Para.update(   
        {'Current total electrolyte volume in jelly roll [m3]':  Vol_Elely_JR_new  }, 
        check_already_exists=False)
    Para.update(   
        {'Ratio of electrolyte dry out in jelly roll':Ratio_Dryout}, 
        check_already_exists=False)
    Para.update(   {'Electrode width [m]':Width_new})    
    Para.update(   
        {'Current solvent concentration in the reservoir [mol.m-3]':c_EC_r_new}, 
        check_already_exists=False)     
    Para.update(   
        {'Current electrolyte concentration in the reservoir [mol.m-3]':c_e_r_new}, 
        check_already_exists=False)             
    return Data_Pack,Para

def Get_Last_state(Model, Sol):
    dict_short = {}
    list_short = []
    # update 220808 to satisfy random model option:
    for var, equation in Model.initial_conditions.items():
        #print(var._name)
        list_short.append(var._name)
    # delete Porosity times concentration and Electrolyte potential then add back
    list_short.remove("Porosity times concentration [mol.m-3]");
    list_short.remove("Electrolyte potential [V]");
    list_short.extend(
        ("Negative electrode porosity times concentration [mol.m-3]",
        "Separator porosity times concentration [mol.m-3]",
        "Positive electrode porosity times concentration [mol.m-3]",   

        "Negative electrolyte potential [V]",
        "Separator electrolyte potential [V]",
        "Positive electrolyte potential [V]",))
    for list_short_i in list_short:
        dict_short.update( { list_short_i : Sol.last_state[list_short_i].data  } )
    return list_short,dict_short


# Define a function to calculate based on previous solution
def Run_Model_Base_On_Last_Solution( 
    Model  , Sol , Para_update, ModelExperiment, 
    Update_Cycles,Temper_i ,mesh_list,submesh_strech):
    # Use Sulzer's method: inplace = false
    # Important line: define new model based on previous solution
    Ratio_CeLi = Para_update[
        "Ratio of Li-ion concentration change in electrolyte consider solvent consumption"]
    if isinstance(Sol, pb.solvers.solution.Solution):
        list_short,dict_short = Get_Last_state(Model, Sol)
    elif isinstance(Sol, list):
        [dict_short, getSth] = Sol
        list_short = []
        for var, equation in Model.initial_conditions.items():
            list_short.append(var._name)
    else:
        print("!! Big problem, Sol here is neither solution or list")

    dict_short["Negative electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Negative electrode porosity times concentration [mol.m-3]"] * Ratio_CeLi )# important: update sol here!
    dict_short["Separator porosity times concentration [mol.m-3]"] = (
        dict_short["Separator porosity times concentration [mol.m-3]"] * Ratio_CeLi )# important: update sol here!
    dict_short["Positive electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Positive electrode porosity times concentration [mol.m-3]"] * Ratio_CeLi )# important: update sol here!
    Model_new = Model.set_initial_conditions_from(dict_short, inplace=False)
    Para_update.update(   {'Ambient temperature [K]':Temper_i });   # run model at 45 degree C
    
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(mesh_list[0]),  
        var.x_s: int(mesh_list[1]),  
        var.x_p: int(mesh_list[2]),  
        var.r_n: int(mesh_list[3]),  
        var.r_p: int(mesh_list[4]),  
        }
    submesh_types = Model.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": int(submesh_strech)})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 
    Call_Age = RioCallback()  # define callback
    
    # update 231208 - try 3 times until give up 
    i_run_try = 0
    str_err = "Initialize only"
    while i_run_try<3:
        try:
            if i_run_try < 2:
                Simnew = pb.Simulation(
                    Model_new,
                    experiment = ModelExperiment, 
                    parameter_values=Para_update, 
                    solver = pb.CasadiSolver(),
                    var_pts = var_pts,
                    submesh_types=submesh_types )
                Sol_new = Simnew.solve(
                    calc_esoh=False,
                    save_at_cycles = Update_Cycles,
                    callbacks=Call_Age)
                if Call_Age.success == False:
                    raise Experiment_error_infeasible("Self detect")
            else:   # i_run_try = 2, 
                Simnew = pb.Simulation(
                    Model_new,
                    experiment = ModelExperiment, 
                    parameter_values=Para_update, 
                    solver = pb.CasadiSolver(return_solution_if_failed_early=True),
                    var_pts = var_pts,
                    submesh_types=submesh_types )
                Sol_new = Simnew.solve(
                    calc_esoh=False,
                    # save_at_cycles = Update_Cycles,
                    callbacks=Call_Age)
                Succeed_AGE_cycs = len(Sol_new.cycles)
                if Succeed_AGE_cycs < Update_Cycles:
                    str_err = f"Partially succeed to run the ageing set for {Succeed_AGE_cycs} cycles the {i_run_try+1}th time"
                    print(str_err)
                else: 
                    str_err = f"Fully succeed to run the ageing set for {Succeed_AGE_cycs} cycles the {i_run_try+1}th time"
                    print(str_err)
        except (
            pb.expression_tree.exceptions.ModelError,
            pb.expression_tree.exceptions.SolverError
            ) as e:
            i_run_try += 1
            Sol_new = "Model error or solver error"
            str_err = f"Fail to run the ageing set due to {Sol_new} for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        except Experiment_error_infeasible as custom_error:
            i_run_try += 1
            Sol_new = "Experiment error or infeasible"
            DeBug_List = [
                Model, Model_new, Call_Age, Simnew,  Sol , Sol_new, Para_update, ModelExperiment, 
                Update_Cycles,Temper_i ,mesh_list,submesh_strech, var_pts,
                submesh_types,  list_short, dict_short, 
            ]
            str_err= f"{Sol_new}: {custom_error} for ageing set for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        else:
            i_run_try += 1
            # add 230221 - update 230317 try to access Throughput capacity more than once
            i_try = 0
            while i_try<3:
                try:
                    getSth2 = Sol_new['Throughput capacity [A.h]'].entries[-1]
                except:
                    i_try += 1
                    print(f"Fail to read Throughput capacity for the {i_try}th time")
                else:
                    break
            # update 240110
            if isinstance(Sol, pb.solvers.solution.Solution):
                i_try = 0
                while i_try<3:
                    try:
                        getSth = Sol['Throughput capacity [A.h]'].entries[-1]
                    except:
                        i_try += 1
                        print(f"Fail to read Throughput capacity for the {i_try}th time")
                    else:
                        break
            elif isinstance(Sol, list):
                print("Must be the first run after restart, already have getSth")
            else:
                print("!! Big problem, Sol here is neither solution or list")
            # update 23-11-16 change method to get throughput capacity:
            # # (1) add old ones, as before; (2): change values of the last one
            Sol_new['Throughput capacity [A.h]'].entries += getSth
            cyc_number = len(Sol_new.cycles)
            if not Update_Cycles == 1: # the solution is imcomplete in this case
                thr_1st = np.trapz(
                    abs(Sol_new.cycles[0]["Current [A]"].entries), 
                    Sol_new.cycles[0]["Time [h]"].entries) # in A.h
                thr_end = np.trapz(
                    abs(Sol_new.cycles[-1]["Current [A]"].entries), 
                    Sol_new.cycles[-1]["Time [h]"].entries) # in A.h
                thr_tot = (thr_1st+thr_end) / 2 * cyc_number
            else:
                thr_tot = abs(
                    Sol_new['Throughput capacity [A.h]'].entries[-1] - 
                    Sol_new['Throughput capacity [A.h]'].entries[0]  )
            Sol_new['Throughput capacity [A.h]'].entries[-1] = getSth + thr_tot # only the last one is true
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
            print(f"Succeed to run the ageing set for {cyc_number} cycles the {i_run_try}th time")
            break # terminate the loop of trying to solve if you can get here. 
    
    # print("Solved this model in {}".format(ModelTimer.time()))
    Result_list = [Model_new, Sol_new,Call_Age,DeBug_List]
    return Result_list

# Update 231117 new method to get throughput capacity to avoid problems of empty solution
def Get_ThrCap(sol_PRT): 
    thr_tot = 0
    for cycle in sol_PRT.cycles:
        for step in cycle.steps:
            # print(type(step))
            if not isinstance(step,pb.solvers.solution.EmptySolution):
                thr_i = np.trapz(
                    abs(step["Current [A]"].entries), 
                    step["Time [h]"].entries)
                thr_tot += thr_i
    return thr_tot

def Run_Model_Base_On_Last_Solution_RPT( 
    Model  , Sol,  Para_update, 
    ModelExperiment ,Update_Cycles, Temper_i,mesh_list,submesh_strech):
    # Use Sulzer's method: inplace = false
    Ratio_CeLi = Para_update["Ratio of Li-ion concentration change in electrolyte consider solvent consumption"]
    # print("Model is now using average EC Concentration of:",Para_update['Bulk solvent concentration [mol.m-3]'])
    # print("Ratio of electrolyte dry out in jelly roll is:",Para_update['Ratio of electrolyte dry out in jelly roll'])
    # print("Model is now using an electrode width of:",Para_update['Electrode width [m]'])
    # Important line: define new model based on previous solution
    list_short,dict_short = Get_Last_state(Model, Sol)
    dict_short["Negative electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Negative electrode porosity times concentration [mol.m-3]"] * Ratio_CeLi )# important: update sol here!
    dict_short["Separator porosity times concentration [mol.m-3]"] = (
        dict_short["Separator porosity times concentration [mol.m-3]"] * Ratio_CeLi )# important: update sol here!
    dict_short["Positive electrode porosity times concentration [mol.m-3]"] = (
        dict_short["Positive electrode porosity times concentration [mol.m-3]"] * Ratio_CeLi )# important: update sol here!
    Model_new = Model.set_initial_conditions_from(dict_short, inplace=False)
    Para_update.update(   {'Ambient temperature [K]':Temper_i });
    
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(mesh_list[0]),  
        var.x_s: int(mesh_list[1]),  
        var.x_p: int(mesh_list[2]),  
        var.r_n: int(mesh_list[3]),  
        var.r_p: int(mesh_list[4]),  }
    submesh_types = Model.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": int(submesh_strech)})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 
    Simnew = pb.Simulation(
        Model_new,
        experiment = ModelExperiment, 
        parameter_values=Para_update, 
        solver = pb.CasadiSolver(),
        var_pts = var_pts,
        submesh_types=submesh_types
    )
    Call_RPT = RioCallback()  # define callback
    # update 231208 - try 3 times until give up 
    i_run_try = 0
    while i_run_try<3:
        try:
            Sol_new = Simnew.solve(
                calc_esoh=False,
                # save_at_cycles = Update_Cycles,
                callbacks=Call_RPT) 
            if Call_RPT.success == False:
                raise Experiment_error_infeasible("Self detect")        
        except (
            pb.expression_tree.exceptions.ModelError,
            pb.expression_tree.exceptions.SolverError
            ) as e:
            i_run_try += 1
            Sol_new = "Model error or solver error"
            str_err = f"Fail to run RPT due to {Sol_new} for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        except Experiment_error_infeasible as custom_error:
            i_run_try += 1
            Sol_new = "Experiment error or infeasible"
            DeBug_List = [
                Model, Model_new, Call_RPT, Simnew,   Sol , Sol_new, Para_update, ModelExperiment, 
                Update_Cycles,Temper_i ,mesh_list,submesh_strech, var_pts,
                submesh_types,  list_short, dict_short, 
            ]
            str_err = f"{Sol_new}: {custom_error} for RPT for the {i_run_try}th time"
            print(str_err)
            DeBug_List = [ Para_update, Update_Cycles, dict_short, str_err ]
        else:
            # add 230221 - update 230317 try to access Throughput capacity more than once
            i_run_try += 1
            i_try = 0
            while i_try<3:
                try:
                    getSth2 = Sol_new['Throughput capacity [A.h]'].entries[-1]
                except:
                    i_try += 1
                    print(f"Fail to read Throughput capacity for the {i_try}th time")
                else:
                    break
            i_try = 0
            while i_try<3:
                try:
                    getSth = Sol['Throughput capacity [A.h]'].entries[-1]
                except:
                    i_try += 1
                    print(f"Fail to read Throughput capacity for the {i_try}th time")
                else:
                    break
            Sol_new['Throughput capacity [A.h]'].entries += getSth
            # update 23-11-17 change method to get throughput capacity to avioid problems of empty solution:
            thr_tot = Get_ThrCap(Sol_new)
            Sol_new['Throughput capacity [A.h]'].entries[-1] = getSth + thr_tot # only the last one is true
            DeBug_List = "Empty"
            print(f"Succeed to run RPT for the {i_run_try}th time")
            break # terminate the loop of trying to solve if you can get here. 
    # print("Solved this model in {}".format(ModelTimer.time()))
    Result_List_RPT = [Model_new, Sol_new,Call_RPT,DeBug_List]
    return Result_List_RPT



def recursive_scan(mylist,kvs, key_list, acc):
    # 递归终止条件
    if len(key_list) == 0:
        mylist.append(acc.copy())   # copy.deepcopy(acc) 如果value是个list，就要deep copy了
        return mylist
    # 继续递归
    k = key_list[0]
    for v in kvs[k]:
        acc[k] = v
        # print(v)
        recursive_scan(mylist,kvs, key_list[1:], acc)

# Update 2023-10-04 
def Overwrite_Initial_L_SEI_0_Neg_Porosity(Para_0,cap_loss):
    """ 
    This is to overwrite the initial negative electrode porosity 
    and initial SEI thickness (inner, outer) to be consistent 
    with the initial capacity loss 
    """
    delta_Q_SEI = cap_loss * 3600
    V_SEI = Para_0["Outer SEI partial molar volume [m3.mol-1]"] 
    # do this when finish updating
    F = 96485.3
    A = Para_0["Electrode width [m]"] * Para_0["Electrode height [m]"]
    z_SEI = Para_0["Ratio of lithium moles to SEI moles"] 
    L_neg = Para_0["Negative electrode thickness [m]"] 
    eps_act_neg = Para_0["Negative electrode active material volume fraction"]
    R_neg =   Para_0["Negative particle radius [m]"]
    l_cr_init = Para_0["Negative electrode initial crack length [m]"]
    w_cr = Para_0["Negative electrode initial crack width [m]"]
    rho_cr = Para_0["Negative electrode number of cracks per unit area [m-2]"]
    a_neg = (3 * eps_act_neg / R_neg)
    roughness = 1 + 2 * l_cr_init * w_cr * rho_cr
    L_SEI_init = delta_Q_SEI  * V_SEI / (
        z_SEI * F * A * L_neg * a_neg * roughness)

    delta_epi = (L_SEI_init ) * roughness * a_neg
    L_inner_init = L_SEI_init / 2
    epi = 0.25 - delta_epi
    # print(L_inner_init,epi)
    # important: update here!
    Para_0["Negative electrode porosity"] = epi
    Para_0["Initial outer SEI thickness [m]"] = L_inner_init + 2.5e-9
    Para_0["Initial inner SEI thickness [m]"] = L_inner_init + 2.5e-9
    print(f"Has Overwritten Initial outer SEI thickness [m] to be {(L_inner_init+2.5e-9):.2e} and Negative electrode porosity to be {epi:.3f} to account for initial capacity loss of {cap_loss:.3f} Ah")

    return Para_0












































# this function is to initialize the para with a known dict
def Para_init(Para_dict,cap_st):
    Para_dict_used = Para_dict.copy()
    






    Para_0=pb.ParameterValues(Para_dict_used["Para_Set"]  )
    Para_dict_used.pop("Para_Set")
    import ast,json

    if Para_dict_used.__contains__("Scan No"):
        Para_dict_used.pop("Scan No")
    if Para_dict_used.__contains__("Exp No."):
        Para_dict_used.pop("Exp No.")

    if Para_dict_used.__contains__("Total ageing cycles"):
        Total_Cycles = Para_dict_used["Total ageing cycles"]  
        Para_dict_used.pop("Total ageing cycles")
    if Para_dict_used.__contains__("Ageing cycles between RPT"):
        Cycle_bt_RPT = Para_dict_used["Ageing cycles between RPT"]  
        Para_dict_used.pop("Ageing cycles between RPT")
    if Para_dict_used.__contains__("Update cycles for ageing"):
        Update_Cycles = Para_dict_used["Update cycles for ageing"]  
        Para_dict_used.pop("Update cycles for ageing")
    if Para_dict_used.__contains__("Cycles within RPT"):
        RPT_Cycles = Para_dict_used["Cycles within RPT"]  
        Para_dict_used.pop("Cycles within RPT")
    if Para_dict_used.__contains__("Ageing temperature"):
        Temper_i = Para_dict_used["Ageing temperature"]  + 273.15 # update: change to K
        Para_dict_used.pop("Ageing temperature")
    if Para_dict_used.__contains__("RPT temperature"):
        Temper_RPT = Para_dict_used["RPT temperature"] + 273.15 # update: change to K 
        Para_dict_used.pop("RPT temperature")
    if Para_dict_used.__contains__("Mesh list"):
        mesh_list = Para_dict_used["Mesh list"]
        mesh_list = json.loads(mesh_list)  
        Para_dict_used.pop("Mesh list")
    if Para_dict_used.__contains__("Exponential mesh stretch"):
        submesh_strech = Para_dict_used["Exponential mesh stretch"]  
        Para_dict_used.pop("Exponential mesh stretch")
    else:
        submesh_strech = "nan";
    if Para_dict_used.__contains__("Model option"):
        model_options = Para_dict_used["Model option"] 
        # careful: need to change str to dict 
        model_options = ast.literal_eval(model_options)
        Para_dict_used.pop("Model option")
    
    
    
    
    
    
    cap_increase = 0
    
    




    if Para_dict_used.__contains__("Initial Neg SOC"):
        c_Neg1SOC_in = (
            Para_0["Maximum concentration in negative electrode [mol.m-3]"]
            *Para_dict_used["Initial Neg SOC"]  )
        Para_0.update(
            {"Initial concentration in negative electrode [mol.m-3]":
            c_Neg1SOC_in})
        Para_dict_used.pop("Initial Neg SOC")
    if Para_dict_used.__contains__("Initial Pos SOC"):    
        c_Pos1SOC_in = (
            Para_0["Maximum concentration in positive electrode [mol.m-3]"]
            *Para_dict_used["Initial Pos SOC"]  )
        Para_0.update(
            {"Initial concentration in positive electrode [mol.m-3]":
            c_Pos1SOC_in})
        Para_dict_used.pop("Initial Pos SOC")
    # 'Outer SEI partial molar volume [m3.mol-1]'
    if Para_dict_used.__contains__(
        "Outer SEI partial molar volume [m3.mol-1]"): 
        Para_0.update(
            {"Outer SEI partial molar volume [m3.mol-1]":
            Para_dict_used["Outer SEI partial molar volume [m3.mol-1]"]})
        Para_0.update(
            {"Inner SEI partial molar volume [m3.mol-1]":
            Para_dict_used["Outer SEI partial molar volume [m3.mol-1]"]})
        
        Para_dict_used.pop("Outer SEI partial molar volume [m3.mol-1]")


    CyclePack = [ 
        Total_Cycles,Cycle_bt_RPT,Update_Cycles,RPT_Cycles,
        Temper_i,Temper_RPT,mesh_list,submesh_strech,
        model_options,     cap_increase]
    # Mark Ruihe - updated 230222 - from P3
    for key, value in Para_dict_used.items():
        # risk: will update parameter that doesn't exist, 
        # so need to make sure the name is right 
        if isinstance(value, str):
            Para_0.update({key: eval(value)})
            #Para_dict_used.pop(key)
        else:
            Para_0.update({key: value},check_already_exists=False)




    cap_loss =   5.0 - 4.86491 # cap_st      # change for debug=   4.86491  
    Para_0 = Overwrite_Initial_L_SEI_0_Neg_Porosity(Para_0,cap_loss)



    return CyclePack,Para_0

# Add 220808 - to simplify the post-processing
def GetSol_dict (my_dict, keys_all, Sol, 
    cycle_no,step_CD, step_CC , step_RE, step_CV ):

    [keys_loc,keys_tim,keys_cyc]  = keys_all;
    # get time_based variables:
    if len(keys_tim): 
        for key in keys_tim:
            step_no = eval("step_{}".format(key[0:2]))
            if key[3:] == "Time [h]":
                my_dict[key].append  (  (
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries
                    -
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries[0]).tolist()  )
            elif key[3:] == "Anode potential [V]":
                sol_step = Sol.cycles[cycle_no].steps[step_no]
                mesh_sep = len(sol_step["Separator electrolyte potential [V]"].entries[:,0])
                if mesh_sep % 2 ==0:
                    phi_ref = (
                        sol_step["Separator electrolyte potential [V]"].entries[mesh_sep//2-1,:]
                        +
                        sol_step["Separator electrolyte potential [V]"].entries[mesh_sep//2+1,:]
                        ) / 2
                    #print(mesh_sep/2+0.5)
                else:
                    phi_ref = sol_step["Separator electrolyte potential [V]"].entries[mesh_sep//2,:]
                V_n = -phi_ref 
                my_dict[key].append(  V_n.tolist()  )
            elif key[3:] == "Cathode potential [V]":
                sol_step = Sol.cycles[cycle_no].steps[step_no]
                mesh_sep = len(sol_step["Separator electrolyte potential [V]"].entries[:,0])
                if mesh_sep % 2 ==0:
                    phi_ref = (
                        sol_step["Separator electrolyte potential [V]"].entries[mesh_sep//2-1,:]
                        +
                        sol_step["Separator electrolyte potential [V]"].entries[mesh_sep//2+1,:]
                        ) / 2
                    #print(mesh_sep/2+0.5)
                else:
                    phi_ref = sol_step["Separator electrolyte potential [V]"].entries[mesh_sep//2,:]
                V =  sol_step["Terminal voltage [V]"].entries
                V_p = V - phi_ref   
                my_dict[key].append(  V_p.tolist()  )
            else:
                my_dict[key].append(  (
                    Sol.cycles[cycle_no].steps[step_no][key[3:]].entries).tolist()  )
    # get cycle_step_based variables: # isn't an array
    if len(keys_cyc): 
        for key in keys_cyc:
            if key in ["Discharge capacity [A.h]"]:
                step_no = step_CD;
                my_dict[key].append  (
                    Sol.cycles[cycle_no].steps[step_no][key].entries[-1]
                    - 
                    Sol.cycles[cycle_no].steps[step_no][key].entries[0])
            elif key in ["Throughput capacity [A.h]"]: # 
                i_try = 0
                while i_try<3:
                    try:
                        getSth = Sol[key].entries[-1]
                    except:
                        i_try += 1
                        print(f"Fail to read Throughput capacity for the {i_try}th time")
                    else:
                        break
                my_dict[key].append(abs(getSth))
            elif key[0:5] in ["CDend","CCend","CVend","REend",
                "CDsta","CCsta","CVsta","REsta",]:
                step_no = eval("step_{}".format(key[0:2]))
                if key[2:5] == "sta":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[0])
                elif key[2:5] == "end":
                    my_dict[key].append  (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[-1])
                
    # get location_based variables:
    if len(keys_loc): 
        for key in keys_loc:
            if key in ["x_n [m]","x [m]","x_s [m]","x_p [m]"]:
                #print("These variables only add once")
                if not len(my_dict[key]):   # special: add only once
                    my_dict[key] = (Sol[key].entries[:,-1]).tolist()
            elif key[0:5] in ["CDend","CCend","CVend","REend",
                            "CDsta","CCsta","CVsta","REsta",]:      
                #print("These variables add multiple times")
                step_no = eval("step_{}".format(key[0:2]))
                if key[2:5] == "sta":
                    my_dict[key].append  ((
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[:,0]).tolist() )
                elif key[2:5] == "end":
                    my_dict[key].append ( (
                        Sol.cycles[cycle_no].steps[step_no][key[6:]].entries[:,-1]).tolist()  )
    return my_dict                              

def Get_SOH_LLI_LAM(my_dict_RPT,model_options,DryOut,mdic_dry,cap_0):
    my_dict_RPT['Throughput capacity [kA.h]'] = (
        np.array(my_dict_RPT['Throughput capacity [A.h]'])/1e3).tolist()
    my_dict_RPT['CDend SOH [%]'] = ((
        np.array(my_dict_RPT["Discharge capacity [A.h]"])
        /
        cap_0       # my_dict_RPT["Discharge capacity [A.h]"][0] # Mark: change to this so that every case has same standard for reservior paper
        )*100).tolist()
    my_dict_RPT["CDend LAM_ne [%]"] = ((1-
        np.array(my_dict_RPT['CDend Negative electrode capacity [A.h]'])
        /my_dict_RPT['CDend Negative electrode capacity [A.h]'][0])*100).tolist()
    my_dict_RPT["CDend LAM_pe [%]"] = ((1-
        np.array(my_dict_RPT['CDend Positive electrode capacity [A.h]'])
        /my_dict_RPT['CDend Positive electrode capacity [A.h]'][0])*100).tolist()
    my_dict_RPT["CDend LLI [%]"] = ((1-
        np.array(my_dict_RPT["CDend Total lithium capacity in particles [A.h]"])
        /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
    if model_options.__contains__("SEI"):
        my_dict_RPT["CDend LLI SEI [%]"] = ((
            np.array(
                my_dict_RPT["CDend Loss of capacity to negative SEI [A.h]"]-
                my_dict_RPT["CDend Loss of capacity to negative SEI [A.h]"][0]
                )
            /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
    else:
        my_dict_RPT["CDend LLI SEI [%]"] = (
            np.zeros(
                np.size(
                    my_dict_RPT["CDend LLI [%]"]))).tolist()
    if model_options.__contains__("SEI on cracks"):
        my_dict_RPT["CDend LLI SEI on cracks [%]"] = ((
            np.array(
                my_dict_RPT["CDend Loss of capacity to negative SEI on cracks [A.h]"]-
                my_dict_RPT["CDend Loss of capacity to negative SEI on cracks [A.h]"][0]
                )
            /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
    else:
        my_dict_RPT["CDend LLI SEI on cracks [%]"] = (
            np.zeros(
                np.size(
                    my_dict_RPT["CDend LLI [%]"]))).tolist()
    if model_options.__contains__("lithium plating"):
        my_dict_RPT["CDend LLI lithium plating [%]"] = ((
            np.array(
                my_dict_RPT["CDend Loss of capacity to negative lithium plating [A.h]"]-
                my_dict_RPT["CDend Loss of capacity to negative lithium plating [A.h]"][0]
                )
            /my_dict_RPT["CDend Total lithium capacity in particles [A.h]"][0])*100).tolist()
    else:
        my_dict_RPT["CDend LLI lithium plating [%]"] =  (
            np.zeros(
                np.size(
                    my_dict_RPT["CDend LLI [%]"]))).tolist()    
    my_dict_RPT["CDend LLI due to LAM [%]"] = (
        np.array(my_dict_RPT['CDend LLI [%]'])
        -my_dict_RPT["CDend LLI SEI [%]"]
        -my_dict_RPT["CDend LLI lithium plating [%]"]
        -my_dict_RPT["CDend LLI SEI on cracks [%]"] )
    if DryOut == "On":
        LAM_to_Dry   = 100-np.array(
            mdic_dry['Width_all']) /mdic_dry['Width_all'][0]*100
        LAM_to_Dry_end=LAM_to_Dry[-1]
    else:
        LAM_to_Dry_end = 0
    LAM_to_Crack_NE_end=my_dict_RPT['CDend LAM_ne [%]'][-1]-LAM_to_Dry_end
    LAM_to_Crack_PE_end=my_dict_RPT['CDend LAM_pe [%]'][-1]-LAM_to_Dry_end
    my_dict_RPT["LAM_to_Dry [%] end"] = LAM_to_Dry_end
    my_dict_RPT["LAM_to_Crack_NE [%] end"] = LAM_to_Crack_NE_end
    my_dict_RPT["LAM_to_Crack_PE [%] end"] = LAM_to_Crack_PE_end
    return my_dict_RPT

# define the model and run break-in cycle - 
# input parameter: model_options, Experiment_Breakin, Para_0, mesh_list, submesh_strech
# output: Sol_0 , Model_0, Call_Breakin
def Run_Breakin(
    model_options, Experiment_Breakin, 
    Para_0, mesh_list, submesh_strech,cap_increase):

    Model_0 = pb.lithium_ion.DFN(options=model_options)
    # update 220926 - add diffusivity and conductivity as variables:
    c_e = Model_0.variables["Electrolyte concentration [mol.m-3]"]
    T = Model_0.variables["Cell temperature [K]"]
    D_e = Para_0["Electrolyte diffusivity [m2.s-1]"]
    sigma_e = Para_0["Electrolyte conductivity [S.m-1]"]
    Model_0.variables["Electrolyte diffusivity [m2.s-1]"] = D_e(c_e, T)
    Model_0.variables["Electrolyte conductivity [S.m-1]"] = sigma_e(c_e, T)
    var = pb.standard_spatial_vars  
    var_pts = {
        var.x_n: int(mesh_list[0]),  
        var.x_s: int(mesh_list[1]),  
        var.x_p: int(mesh_list[2]),  
        var.r_n: int(mesh_list[3]),  
        var.r_p: int(mesh_list[4]),  }       
    submesh_types = Model_0.default_submesh_types
    if submesh_strech == "nan":
        pass
    else:
        particle_mesh = pb.MeshGenerator(
            pb.Exponential1DSubMesh, 
            submesh_params={"side": "right", "stretch": submesh_strech})
        submesh_types["negative particle"] = particle_mesh
        submesh_types["positive particle"] = particle_mesh 


    """ exp_breakin_text = [ (
        # refill
        "Rest for 10 s",  
        #f"Hold at 4.2V until C/100",
        #"Rest for 1 hours (20 minute period)", 
        ) ] 
    Experiment_Breakin= pb.Experiment( 
            exp_breakin_text * 1) 
    Sim_0    = pb.Simulation(
        Model_0,        experiment = Experiment_Breakin,
        parameter_values = Para_0,
        solver = pb.CasadiSolver(),
        var_pts=var_pts,
        submesh_types=submesh_types) 
    Call_Breakin = RioCallback()
    Sol_0    = Sim_0.solve(calc_esoh=False,callbacks=Call_Breakin) """
    # generate a list of capacity increase perturbation:
    import random
    Cap_in_perturbation = []
    try_no = 4
    for k in range(try_no):
        random_number = 0
        while abs(random_number) <= 1e-6:
            random_number = random.uniform(-1e-5, 1e-5)
        Cap_in_perturbation.append(random_number)
    Cap_in_perturbation[0] = 0.0 # first time must be zero
    Cap_in_perturbation[1] = 3.9116344182112036E-4
    c_s_neg_baseline = Para_0["Initial concentration in negative electrode [mol.m-3]"]
    # update 240603 - try shift neg soc 4 times until give up 
    i_run_try = 0
    while i_run_try<try_no:
        try:





            
            Sim_0    = pb.Simulation(
                Model_0,        experiment = Experiment_Breakin,
                parameter_values = Para_0,
                solver = pb.CasadiSolver(),
                var_pts=var_pts,
                submesh_types=submesh_types) 
            Call_Breakin = RioCallback()    
            Sol_0    = Sim_0.solve(calc_esoh=False,callbacks=Call_Breakin)
        except (
            pb.expression_tree.exceptions.ModelError,
            pb.expression_tree.exceptions.SolverError
            ) as e:
            Sol_0 = "Model error or solver error"
            str_err = (
                f"Fail to run break in due to {Sol_0} for the {i_run_try}th time"
                f" with capacity increase of {cap_increase}Ah and "
                f"perturbation of {Cap_in_perturbation[i_run_try]:.2e}Ah")
            print();print(str_err);print()
            i_run_try += 1
        else:
            str_err = (
                f"Succeed to run break in for the {i_run_try}th time "
                f"with capacity increase of {cap_increase}Ah and "
                f"perturbation of {Cap_in_perturbation[i_run_try]:.2e}Ah")
            print();print(str_err);print()
            break
    Result_list_breakin = [Model_0,Sol_0,Call_Breakin]

    return Result_list_breakin

# Input: Para_0
# Output: mdic_dry, Para_0
def Initialize_mdic_dry(Para_0,Int_ElelyExces_Ratio):
    mdic_dry = {
        "CeEC_All": [],
        "c_EC_r_new_All": [],
        "c_e_r_new_All": [],
        "Ratio_CeEC_All":[],
        "Ratio_CeLi_All":[],
        "Ratio_Dryout_All":[],

        "Vol_Elely_Tot_All": [],
        "Vol_Elely_JR_All":[],
        "Vol_Pore_tot_All": [],        
        "Vol_EC_consumed_All":[],
        "Vol_Elely_need_All":[],
        "Width_all":[],
        "Vol_Elely_add_All":[],
        "Vol_Pore_decrease_All":[],
        "Test_V_All":[],
        "Test_V2_All":[],
    }
    T_0                  =  Para_0['Initial temperature [K]']
    Porosity_Neg_0       =  Para_0['Negative electrode porosity']  
    Porosity_Pos_0       =  Para_0['Positive electrode porosity']  
    Porosity_Sep_0       =  Para_0['Separator porosity']  
    cs_Neg_Max           =  Para_0["Maximum concentration in negative electrode [mol.m-3]"];
    L_p                  =  Para_0["Positive electrode thickness [m]"]
    L_n                  =  Para_0["Negative electrode thickness [m]"]
    L_s                  =  Para_0["Separator thickness [m]"]
    L_y                  =  Para_0["Electrode width [m]"]
    Para_0.update({'Initial Electrode width [m]':L_y}, check_already_exists=False)
    L_y_0                =  Para_0["Initial Electrode width [m]"]
    L_z                  =  Para_0["Electrode height [m]"]
    Para_0.update({'Initial Electrode height [m]':L_z}, check_already_exists=False)
    L_z_0                =  Para_0["Initial Electrode height [m]"]
    
    Vol_Elely_Tot        = (
        ( L_n*Porosity_Neg_0 +  L_p*Porosity_Pos_0  +  L_s*Porosity_Sep_0  )  
        * L_y_0 * L_z_0 * Int_ElelyExces_Ratio
     ) # Set initial electrolyte amount [L] 
    Vol_Elely_JR         =  (
        ( L_n*Porosity_Neg_0 +  L_p*Porosity_Pos_0  +  L_s*Porosity_Sep_0  )  
        * L_y_0 * L_z_0 )
    Vol_Pore_tot         =  (
        ( L_n*Porosity_Neg_0 +  L_p*Porosity_Pos_0  +  L_s*Porosity_Sep_0  )  
        * L_y_0 * L_z_0 )
    Ratio_CeEC           =  1.0; Ratio_CeLi =  1.0  ;Ratio_Dryout         =  1.0
    Vol_EC_consumed      =  0;Vol_Elely_need       =  0;Vol_Elely_add        =  0
    Vol_Pore_decrease    =  0;Test_V2 = 0; Test_V = 0;
    print('Initial electrolyte amount is ', Vol_Elely_Tot*1e6, 'mL') 
    Para_0.update(
        {'Current total electrolyte volume in jelly roll [m3]':
        Vol_Elely_JR}, check_already_exists=False)
    Para_0.update(
        {'Current total electrolyte volume in whole cell [m3]':
        Vol_Elely_Tot}, check_already_exists=False)   
    
    mdic_dry["Vol_Elely_Tot_All"].append(Vol_Elely_Tot*1e6);            
    mdic_dry["Vol_Elely_JR_All"].append(Vol_Elely_JR*1e6);     
    mdic_dry["Vol_Pore_tot_All"].append(Vol_Pore_tot*1e6);           
    mdic_dry["Ratio_CeEC_All"].append(Ratio_CeEC);                      
    mdic_dry["Ratio_CeLi_All"].append(Ratio_CeLi);             
    mdic_dry["Ratio_Dryout_All"].append(Ratio_Dryout);
    mdic_dry["Vol_EC_consumed_All"].append(Vol_EC_consumed*1e6);        
    mdic_dry["Vol_Elely_need_All"].append(Vol_Elely_need*1e6);     
    mdic_dry["Width_all"].append(L_y_0);
    mdic_dry["Vol_Elely_add_All"].append(Vol_Elely_add*1e6);            
    mdic_dry["Vol_Pore_decrease_All"].append(Vol_Pore_decrease*1e6);
    mdic_dry["Test_V_All"].append(Test_V*1e6); 
    mdic_dry["Test_V2_All"].append(Test_V2*1e6); 
    mdic_dry["c_e_r_new_All"].append(Para_0["Initial concentration in electrolyte [mol.m-3]"]); 
    mdic_dry["c_EC_r_new_All"].append(Para_0['EC initial concentration in electrolyte [mol.m-3]'])

    return mdic_dry,Para_0
    
def Update_mdic_dry(Data_Pack,mdic_dry):
    [
        Vol_EC_consumed, 
        Vol_Elely_need, 
        Test_V, 
        Test_V2, 
        Vol_Elely_add, 
        Vol_Elely_Tot_new, 
        Vol_Elely_JR_new, 
        Vol_Pore_tot_new, 
        Vol_Pore_decrease, 
        c_e_r_new, c_EC_r_new,
        Ratio_Dryout, Ratio_CeEC_JR, 
        Ratio_CeLi_JR,
        Width_new, ]= Data_Pack;
    mdic_dry["Vol_Elely_Tot_All"].append(Vol_Elely_Tot_new*1e6);            
    mdic_dry["Vol_Elely_JR_All"].append(Vol_Elely_JR_new*1e6);     
    mdic_dry["Vol_Pore_tot_All"].append(Vol_Pore_tot_new*1e6);           
    mdic_dry["Ratio_CeEC_All"].append(Ratio_CeEC_JR);                      
    mdic_dry["Ratio_CeLi_All"].append(Ratio_CeLi_JR);             
    mdic_dry["Ratio_Dryout_All"].append(Ratio_Dryout);
    mdic_dry["Vol_EC_consumed_All"].append(Vol_EC_consumed*1e6);        
    mdic_dry["Vol_Elely_need_All"].append(Vol_Elely_need*1e6);     
    mdic_dry["Width_all"].append(Width_new);
    mdic_dry["Vol_Elely_add_All"].append(Vol_Elely_add*1e6);            
    mdic_dry["Vol_Pore_decrease_All"].append(Vol_Pore_decrease*1e6);
    mdic_dry["Test_V_All"].append(Test_V*1e6); 
    mdic_dry["Test_V2_All"].append(Test_V2*1e6); 
    mdic_dry["c_e_r_new_All"].append(c_e_r_new); 
    mdic_dry["c_EC_r_new_All"].append(c_EC_r_new)

    return mdic_dry
    

# Function to read exp
def Read_Exp(BasicPath,Exp_Any_Cell,Exp_Path,Exp_head,Exp_Any_Temp,i):
    Exp_Any_AllData  = {}
    for cell in Exp_Any_Cell:
        Exp_Any_AllData[cell] = {} # one cell
        # For extracted directly measured capacity, resistance, etc.
        Exp_Any_AllData[cell]["Extract Data"] = pd.read_csv(
            BasicPath+Exp_Path[i]+
            f"{Exp_head[i]} - cell {cell} ({Exp_Any_Temp[cell]}degC) - Extracted Data.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["Extract Data"].loc[     # update 230617- fill this with avg T
            0, 'Age set average temperature (degC)'] = float(Exp_Any_Temp[cell])
        # Read for DMA results, further a dictionary
        Exp_Any_AllData[cell]["DMA"] = {}
        Exp_Any_AllData[cell]["DMA"]["Cap_Offset"]=pd.read_csv(
            BasicPath+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - Capacity and offset data from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["DMA"]["LLI_LAM"]=pd.read_csv(
            BasicPath+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - DM data from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["DMA"]["Fit_SOC"]=pd.read_csv(
            BasicPath+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - fitting parameters from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        Exp_Any_AllData[cell]["DMA"]["RMSE"]=pd.read_csv(
            BasicPath+Exp_Path[i]+ "DMA Output/" + f"cell {cell}/" + 
            f"{Exp_head[i]} - RMSE data from OCV-fitting for cell {cell}.csv", 
            index_col=0)
        # update 230726: read 0.1C voltage curve, discharge only
        Exp_Any_AllData[cell]["0.1C voltage"] = {}
        for m in range(16):
            try:
                C_10_curve_temp = pd.read_csv(
                    BasicPath+Exp_Path[i]+ "0.1C Voltage Curves/"+ f"cell {cell}/" +
                    f"{Exp_head[i]} - cell {cell} - RPT{m} - 0.1C discharge data.csv", )    
            except:
                pass # print(f"Exp-{i+1} - Cell {cell} doesn't have RPT {m}")
            else:
                C_10_curve_temp["Time (h)"] = (
                    C_10_curve_temp["Time (s)"] - 
                    C_10_curve_temp["Time (s)"].iloc[0]) / 3600
                Exp_Any_AllData[cell]["0.1C voltage"][f"RPT{m}"] = C_10_curve_temp
                # print(f"Read Exp-{i+1} - Cell {cell} RPT {m}")
    print("Finish reading Experiment!")
    return Exp_Any_AllData
# judge ageing shape: - NOT Ready yet!
def check_concave_convex(x_values, y_values):
    # Check if the series has at least three points
    if len(x_values) < 3 or len(y_values) < 3:
        return ["None","None""None"]

    # Calculate the slopes between adjacent points based on "x" and "y" values
    slopes = []
    for i in range(1, len(y_values)):
        slope = (y_values[i] - y_values[i - 1]) / (x_values[i] - x_values[i - 1])
        slopes.append(slope)

    # Determine the indices for the start, middle, and end groups
    start_index = 0
    middle_index = len(y_values) // 2 - 1
    end_index = len(y_values) - 3

    # Check if the start group is concave, convex, or linear
    start_group = slopes[:3]
    start_avg = sum(start_group[:2]) / 2
    start_point = start_group[1]
    # TODO we still need to consider a proper index to judge linear
    if abs(start_avg - start_point) <= 0.005 * start_point:
        start_shape = "Linear"
    elif start_avg <= start_point:
        start_shape = "Sub"
    else:
        start_shape = "Super"

    # Check if the middle group is concave, convex, or linear
    middle_group = slopes[middle_index:middle_index + 3]
    middle_avg = sum(middle_group[:2]) / 2
    middle_point = middle_group[1]
    if abs(middle_avg - middle_point) <= 0.005 * middle_point:
        middle_shape = "Linear"
    elif middle_avg <= middle_point:
        middle_shape = "Sub"
    else:
        middle_shape = "Super"

    # Check if the end group is concave, convex, or linear
    end_group = slopes[end_index:]
    end_avg = sum(end_group[:2]) / 2
    end_point = end_group[1]
    if abs(end_avg - end_point) <= 0.005 * end_point:
        end_shape = "Linear"
    elif end_avg <= end_point:
        end_shape = "Sub"
    else:
        end_shape = "Super"

    # Return the shape results
    shape_results = [start_shape, middle_shape, end_shape]
    return shape_results

import shutil
def Copy_png(Root_Path,rows_per_file,Scan_end_end,purpose_i):
    Path_purpose = os.path.join(Root_Path, purpose_i)
    source_folders = []
    for i_bundle in range(int(Scan_end_end/rows_per_file)):
        Scan_start = (i_bundle)*rows_per_file+1;    
        Scan_end   = min(Scan_start + rows_per_file-1, Scan_end_end)    
        purpose = f"{purpose_i}_Case_{Scan_start}_{Scan_end}"
        source_folders.append(purpose)
        #print(purpose)

    # Create the Plot_Collect folder if it doesn't exist
    plot_collect_directory = os.path.join(Path_purpose, "Plot_Collect")
    os.makedirs(plot_collect_directory, exist_ok=True)

    # Move the .png files to the Plot_Collect folder
    for folder in source_folders:
        plots_directory = os.path.join(Path_purpose, folder, "Plots")
        if os.path.exists(plots_directory):
            for filename in os.listdir(plots_directory):
                if filename.startswith("0_") and filename.endswith("Summary.png"):
                    source_file = os.path.join(plots_directory, filename)
                    destination_file = os.path.join(plot_collect_directory, filename)
                    shutil.copy(source_file, destination_file)
    return 


# plot inside the function:
def Plot_Cyc_RPT_4(
        my_dict_RPT, Exp_Any_AllData,Temp_Cell_Exp,
        XY_pack,index_exp, Plot_Exp, R_from_GITT,
        Scan_i,Re_No,Temper_i,model_options,
        BasicPath, Target,fs,dpi):
    
    Num_subplot = 5;
    fig, axs = plt.subplots(2,3, figsize=(15,7.8),tight_layout=True)
    axs[0,0].plot(
        my_dict_RPT['Throughput capacity [kA.h]'], 
        my_dict_RPT['CDend SOH [%]'],     
        '-o', label="Scan=" + str(Scan_i) )
    axs[0,1].plot(
        my_dict_RPT['Throughput capacity [kA.h]'], 
        my_dict_RPT["CDend LLI [%]"],'-o', label="total LLI")
    if model_options.__contains__("lithium plating"):
        axs[0,1].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI lithium plating [%]"],'--o', label="LiP")
    if model_options.__contains__("SEI"):
        axs[0,1].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI SEI [%]"] ,'--o', label="SEI")
    if model_options.__contains__("SEI on cracks"):
        axs[0,1].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI SEI on cracks [%]"] ,
            '--o', label="SEI-on-cracks")
    axs[0,2].plot(
        my_dict_RPT["Throughput capacity [kA.h]"], 
        my_dict_RPT["CDend LAM_ne [%]"],     '-o', ) 
    axs[1,0].plot(
        my_dict_RPT["Throughput capacity [kA.h]"], 
        my_dict_RPT["CDend LAM_pe [%]"],     '-o',  ) 
    # axs[1,1].plot(
    #     my_dict_RPT["Throughput capacity [kA.h]"], 
    #     np.array(my_dict_RPT["Res_midSOC"]),     '-o', ) 
    axs[1,2].plot(
        my_dict_RPT["Throughput capacity [kA.h]"][1:], 
        np.array(my_dict_RPT["avg_Age_T"][1:]),     '-o', ) 
    # Plot Charge Throughput (A.h) vs SOH
    color_exp     = [0, 0, 0, 0.3]; marker_exp     = "v";
    color_exp_Avg = [0, 0, 0, 0.7]; marker_exp_Avg = "s";
    if index_exp in list(np.arange(1,6)) and int(Temper_i- 273.15) in [10,25,40]:
        Exp_temp_i_cell = Temp_Cell_Exp[str(int(Temper_i- 273.15))]
    else:
        Exp_temp_i_cell = "nan"
        Plot_Exp = False

    if Plot_Exp == True:
        for cell in Exp_temp_i_cell:
            df = Exp_Any_AllData[cell]["Extract Data"]
            chThr_temp = np.array(df["Charge Throughput (A.h)"])/1e3
            df_DMA = Exp_Any_AllData[cell]["DMA"]["LLI_LAM"]
            axs[0,0].plot(
                chThr_temp,np.array(df_DMA["SoH"])*100,
                color=color_exp,marker=marker_exp,label=f"Cell {cell}") 
            axs[0,1].plot(
                chThr_temp,np.array(df_DMA["LLI"])*100,
                color=color_exp,marker=marker_exp,label=f"Cell {cell}")  
            axs[0,2].plot(
                chThr_temp,np.array(df_DMA["LAM NE_tot"])*100,
                color=color_exp,marker=marker_exp, )
            axs[1,0].plot(
                chThr_temp,np.array(df_DMA["LAM PE"])*100,
                color=color_exp,marker=marker_exp,)
            # update 230312- plot resistance here
            # Exp_1_AllData["A"]["Extract Data"]["0.1s Resistance (Ohms)"]

            # index_Res = df[df['0.1s Resistance (Ohms)'].le(10)].index
            # axs[1,1].plot(
            #     #df["Days of degradation"][index_Res],
            #     np.array(df["Charge Throughput (A.h)"][index_Res])/1e3,
            #     np.array(df["0.1s Resistance (Ohms)"][index_Res])*1e3,
            #     color=color_exp,marker=marker_exp)
            # axs[1,2].plot(
            #     chThr_temp[1:],
            #     np.array(df["Age set average temperature (degC)"][1:]).astype(float),
            #     color=color_exp,marker=marker_exp,)

        # Update 230518: Plot Experiment Average - at 1 expeirment and 1 temperature
        [X_1_st,X_5_st,Y_1_st_avg,Y_2_st_avg,
            Y_3_st_avg,Y_4_st_avg,Y_5_st_avg,Y_6_st_avg]  = XY_pack
        axs[0,0].plot(
            X_1_st,Y_1_st_avg,color=color_exp_Avg,
            marker=marker_exp_Avg,label=f"Exp-Avg") 
        axs[0,1].plot(
            X_1_st,Y_2_st_avg,color=color_exp_Avg,
            marker=marker_exp_Avg,label=f"Exp-Avg")  
        axs[0,2].plot(
            X_1_st,Y_3_st_avg,color=color_exp_Avg,
            marker=marker_exp_Avg, )
        axs[1,0].plot(
            X_1_st,Y_4_st_avg,
            color=color_exp_Avg,marker=marker_exp_Avg,)
        axs[1,1].plot(
            X_5_st,Y_5_st_avg,
            color=color_exp_Avg,marker=marker_exp_Avg)
        axs[1,2].plot(
            X_1_st[1:],Y_6_st_avg[1:],
            color=color_exp_Avg,marker=marker_exp_Avg,)
    axs[0,0].set_ylabel("SOH %")
    axs[0,1].set_ylabel("LLI %")
    axs[0,2].set_ylabel("LAM NE %")
    axs[1,0].set_ylabel("LAM PE %")
    axs[1,1].set_ylabel(r"Lump resistance [m$\Omega$]")
    axs[1,2].set_ylabel(r"Avg age T [$^\circ$C]")
    axs[0,2].set_xlabel("Charge Throughput (kA.h)")
    axs[1,2].set_xlabel("Charge Throughput (kA.h)")
    axf = axs.flatten()
    for i in range(0,6):
        labels = axf[i].get_xticklabels() + axf[i].get_yticklabels(); 
        [label.set_fontname('DejaVu Sans') for label in labels]
        axf[i].tick_params(labelcolor='k', labelsize=fs, width=1);del labels
    axs[1,1].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
    axs[0,0].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
    axs[0,1].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
    fig.suptitle(
        f"Scan_{Scan_i}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}"
        +r"$^\circ$C - Summary", fontsize=fs+2)
    plt.savefig(
        BasicPath + Target+    "Plots/" +  
        f"0_Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Summary.png", dpi=dpi)
    plt.close()  # close the figure to save RAM

    if model_options.__contains__("SEI on cracks"):
        """ Num_subplot = 2;
        fig, axs = plt.subplots(1,Num_subplot, figsize=(12,4.8),tight_layout=True)
        axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend X-averaged total SEI on cracks thickness [m]"],     '-o', label="Scan=" + str(Scan_i) )
        axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend X-averaged negative electrode roughness ratio"],'-o', label="Scan=" + str(Scan_i) )
        axs[0].set_ylabel("SEI on cracks thickness [m]",   fontdict={'family':'DejaVu Sans','size':fs})
        axs[1].set_ylabel("Roughness ratio",   fontdict={'family':'DejaVu Sans','size':fs})
        for i in range(0,Num_subplot):
            axs[i].set_xlabel("Charge Throughput (kA.h)",   fontdict={'family':'DejaVu Sans','size':fs})
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
        axs[0].set_title("X-avg tot Neg SEI on cracks thickness",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[1].set_title("X-avg Neg roughness ratio",   fontdict={'family':'DejaVu Sans','size':fs+1})
        plt.savefig(BasicPath + Target+"Plots/" +
            f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC - Cracks related_Scan.png", dpi=dpi)
        plt.close()  # close the figure to save RAM """

        Num_subplot = 2;
        fig, axs = plt.subplots(1,Num_subplot, figsize=(12,4.8),tight_layout=True)
        axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend Negative electrode capacity [A.h]"][0]
            -
            my_dict_RPT["CDend Negative electrode capacity [A.h]"],'-o',label="Neg Scan=" + str(Scan_i))
        axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend Positive electrode capacity [A.h]"][0]
            -
            my_dict_RPT["CDend Positive electrode capacity [A.h]"],'-^',label="Pos Scan=" + str(Scan_i))
        """ axs[0].plot(
            my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend X-averaged total SEI on cracks thickness [m]"],                  
            '-o',label="Scan="+ str(Scan_i)) """
        for i in range(0,2):
            axs[i].set_xlabel("Charge Throughput (kA.h)",   fontdict={'family':'DejaVu Sans','size':fs})
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
        #axs[0].set_ylabel("SEI on cracks thickness [m]",   fontdict={'family':'DejaVu Sans','size':fs})
        #axs[0].set_title("CDend X-avg tot SEI on cracks thickness",   fontdict={'family':'DejaVu Sans','size':fs+1})
        for i in range(0,2):
            axs[i].set_xlabel("Charge Throughput (kA.h)",   fontdict={'family':'DejaVu Sans','size':fs})
            axs[i].set_ylabel("Capacity [A.h]",   fontdict={'family':'DejaVu Sans','size':fs})
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)
            axs[i].set_title("LAM of Neg and Pos",   fontdict={'family':'DejaVu Sans','size':fs+1})
        plt.savefig(BasicPath + Target+"Plots/" +
            f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC LAM-IR.png", dpi=dpi)
        plt.close()  # close the figure to save RAM
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDsta Positive electrode stoichiometry"] ,'-o',label="Start" )
    axs[0].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend Positive electrode stoichiometry"] ,'-^',label="End" )
    axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDsta Negative electrode stoichiometry"],'-o',label="Start" )
    axs[1].plot(my_dict_RPT['Throughput capacity [kA.h]'], my_dict_RPT["CDend Negative electrode stoichiometry"],'-^',label="End" )
    for i in range(0,2):
        axs[i].set_xlabel("Charge Throughput (kA.h)",   fontdict={'family':'DejaVu Sans','size':fs})
        axs[i].set_ylabel("Stoichiometry",   fontdict={'family':'DejaVu Sans','size':fs})
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
    axs[0].set_title("Neg Sto. range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Pos Sto. range (Dis)",   fontdict={'family':'DejaVu Sans','size':fs+1})
    plt.savefig(BasicPath + Target+"Plots/"+
        f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC SOC_RPT_dis.png", dpi=dpi) 
    plt.close()  # close the figure to save RAM

    # update 230518: plot resistance in C/2 discharge:
    # N_RPT = len(my_dict_RPT["Res_full"])
    # colormap_i = mpl.cm.get_cmap("gray", 14) 
    # fig, axs = plt.subplots(figsize=(4,3.2),tight_layout=True)
    # for i in range(N_RPT):
    #     axs.plot(
    #         my_dict_RPT["SOC_Res"][i], my_dict_RPT["Res_full"][i] ,
    #         color=colormap_i(i),marker="o",  label=f"RPT {i}" )
    # if R_from_GITT: 
    #     axs.set_xlabel("SOC-GITT %",   fontdict={'family':'DejaVu Sans','size':fs})
    #     axs.set_ylabel(r'Res GITT (m$\Omega$)',   fontdict={'family':'DejaVu Sans','size':fs})
    #     # axs.legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
    #     axs.set_title("Res during GITT Dis",   fontdict={'family':'DejaVu Sans','size':fs+1})
    # else:
    #     axs.set_xlabel("SOC-C/2 %",   fontdict={'family':'DejaVu Sans','size':fs})
    #     axs.set_ylabel(r'Res C/2 (m$\Omega$)',   fontdict={'family':'DejaVu Sans','size':fs})
    #     # axs.legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)     
    #     axs.set_title("Res during C/2 Dis",   fontdict={'family':'DejaVu Sans','size':fs+1})
    # plt.savefig(BasicPath + Target+ "Plots/"+
    #     f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Res_full.png", dpi=dpi) 
    # plt.close()  # close the figure to save RAM

    return

def adjust_category_names(categories):
    return ['\n'.join(category.split()) if len(category) > 12 else category for category in categories]


def Plot_DMA_Dec(my_dict_RPT,Scan_i,Re_No,Temper_i,
                 model_options,BasicPath, Target,fs,dpi):
    fig, ax = plt.subplots(figsize=(6,5),tight_layout=True) 
    categories = ['SEI', 'SEI on cracks', 'Li plating', 'LAM (cracking and dry-out)']
    adjusted_categories = adjust_category_names(categories)
    values = [
        my_dict_RPT["CDend LLI SEI [%]"][-1], 
        my_dict_RPT["CDend LLI SEI on cracks [%]"][-1], 
        my_dict_RPT["CDend LLI lithium plating [%]"][-1] , 
        my_dict_RPT["CDend LLI due to LAM [%]"][-1]   ]
    plt.bar(adjusted_categories, values ,width=0.45 )
    plt.ylabel('LLI %')
    fig.suptitle(
        f"Scan_{Scan_i}_Re_{Re_No}-{str(int(Temper_i- 273.15))}"
        +r"$^\circ$C - LLI break down", fontsize=fs+2)
    plt.savefig(
        BasicPath + Target+    "Plots/" +  
        f"0_Scan_{Scan_i}_Re_{Re_No}-{str(int(Temper_i- 273.15))}degC LLI break down.png", dpi=dpi)
    plt.close()  # close the figure to save RAM

    fig, axs = plt.subplots(1,2, figsize=(12,4.5),tight_layout=True) 
    width = 0.45
    categories = ['Dry out', 'Stress',]
    values_ne = [my_dict_RPT["LAM_to_Dry [%] end"],my_dict_RPT["LAM_to_Crack_NE [%] end"] ]
    values_pe = [my_dict_RPT["LAM_to_Dry [%] end"],my_dict_RPT["LAM_to_Crack_PE [%] end"] ]
    axs[0].bar(categories, values_ne, width )
    axs[1].bar(categories, values_pe, width )
    axs[0].set_ylabel("LAM NE %")
    axs[1].set_ylabel("LAM PE %")
    fig.suptitle(
        f"Scan_{Scan_i}_Re_{Re_No}-{str(int(Temper_i- 273.15))}"
        +r"$^\circ$C - LAM break down", fontsize=fs+2)
    plt.savefig(
        BasicPath + Target+    "Plots/" +  
        f"0_Scan_{Scan_i}_Re_{Re_No}-{str(int(Temper_i- 273.15))}degC LAM break down.png", dpi=dpi)
    plt.close()  # close the figure to save RAM
    return

#
"""    
    "CD Time [h]",
    "CD Terminal voltage [V]",
    "CD Anode potential [V]",    # self defined
    "CD Cathode potential [V]",  # self defined
    "CC Time [h]",
    "CC Terminal voltage [V]",
    "CC Anode potential [V]",    # self defined
    "CC Cathode potential [V]",  # self defined """

# def Plot_HalfCell_V(
#         my_dict_RPT,my_dict_AGE,Scan_i,Re_No,index_exp,colormap,
#         Temper_i,model_options,BasicPath, Target,fs,dpi):
#     #~~~~~~~~~~~~~~~~~ plot RPT 
#     def inFun_Plot(my_dict,str_jj):
#         fig, axs = plt.subplots(3,2, figsize=(12,10),tight_layout=True)
#         Str_Front= ["CD ","CC ",]
#         Str_Back = ["Cathode potential [V]","Terminal voltage [V]","Anode potential [V]",]
#         for j in range(2): 
#             for i in range(3):
#                 Time = my_dict[Str_Front[j]+"Time [h]"]
#                 Num_Lines = len(Time)
#                 cmap = mpl.cm.get_cmap(colormap, Num_Lines) # cmap(i)
#                 for k in range(Num_Lines):
#                     Y = my_dict[Str_Front[j]+Str_Back[i]] [k]
#                     axs[i,j].plot( Time[k] , Y, color = cmap(k)   )
#                 if j == 0 :
#                     axs[i,j].set_ylabel(
#                         Str_Back[i],   fontdict={'family':'DejaVu Sans','size':fs})
#             axs[2,j].set_xlabel("Time [h]",   fontdict={'family':'DejaVu Sans','size':fs})
#         axs[0,0].set_title("During Discharge",   fontdict={'family':'DejaVu Sans','size':fs+1})
#         axs[0,1].set_title("During Charge",   fontdict={'family':'DejaVu Sans','size':fs+1})
#         fig.suptitle(
#             f"Scan {str(Scan_i)}-Exp-{index_exp}-{str(int(Temper_i-273.15))}"
#             +r"$^\circ$C"+f" - Half cell Potential ({str_jj})", fontsize=fs+2)
#         plt.savefig(BasicPath + Target+"Plots/" +
#             f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC"
#             f" Half cell Potential ({str_jj}).png", dpi=dpi) 
#         plt.close() 
#     inFun_Plot(my_dict_RPT,"RPT")
#     inFun_Plot(my_dict_AGE,"AGE")
#     return 

# refined 241021
def Plot_HalfCell_V(
        my_dict_RPT, my_dict_AGE, Scan_i, Re_No, index_exp, colormap,
        Temper_i, model_options, BasicPath, Target, fs, dpi):
    #~~~~~~~~~~~~~~~~~ plot RPT 
    def inFun_Plot(my_dict, str_jj):
        fig, axs = plt.subplots(3, 2, figsize=(12, 10), tight_layout=True)
        Str_Front = ["CD ", "CC "]
        Str_Back = ["Cathode potential [V]", "Terminal voltage [V]", "Anode potential [V]"]
        
        for j in range(2): 
            for i in range(3):
                Time = my_dict[Str_Front[j] + "Time [h]"]
                Num_Lines = len(Time)
                
                # 使用 plt.get_cmap 获取颜色映射
                cmap = plt.get_cmap(colormap)  # 移除 Num_Lines 参数
                
                for k in range(Num_Lines):
                    Y = my_dict[Str_Front[j] + Str_Back[i]][k]
                    # 使用颜色映射为每条线赋颜色
                    axs[i, j].plot(Time[k], Y, color=cmap(k / Num_Lines))  # 使用归一化的 k 值
                    
                if j == 0:
                    axs[i, j].set_ylabel(Str_Back[i], fontdict={'family':'DejaVu Sans', 'size':fs})
                    
            axs[2, j].set_xlabel("Time [h]", fontdict={'family':'DejaVu Sans', 'size':fs})
        
        axs[0, 0].set_title("During Discharge", fontdict={'family':'DejaVu Sans', 'size':fs+1})
        axs[0, 1].set_title("During Charge", fontdict={'family':'DejaVu Sans', 'size':fs+1})
        
        fig.suptitle(
            f"Scan {str(Scan_i)}-Exp-{index_exp}-{str(int(Temper_i-273.15))}"
            + r"$^\circ$C" + f" - Half cell Potential ({str_jj})", fontsize=fs+2)
        
        # 保存图片
        plt.savefig(BasicPath + Target + "Plots/" +
                    f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i - 273.15))}degC"
                    f" Half cell Potential ({str_jj}).png", dpi=dpi)
        plt.close()

    # 绘制 RPT 和 AGE
    inFun_Plot(my_dict_RPT, "RPT")
    inFun_Plot(my_dict_AGE, "AGE")




def Plot_Loc_AGE_4(my_dict_AGE,Scan_i,Re_No,index_exp,Temper_i,model_options,BasicPath, Target,fs,dpi):
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Porosity"][0],'-o',label="First")
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Porosity"][-1],'-^',label="Last"  )
    axs[1].plot(
        my_dict_AGE["x_n [m]"],
        my_dict_AGE["CDend Negative electrode reaction overpotential [V]"][0],'-o',label="First" )
    axs[1].plot(
        my_dict_AGE["x_n [m]"],
        my_dict_AGE["CDend Negative electrode reaction overpotential [V]"][-1],'-^',label="Last" )

    axs[0].set_xlabel("Dimensional Cell thickness",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1].set_xlabel("Dimensional Neg thickness",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[0].set_title("Porosity",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Neg electrode reaction overpotential",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[0].set_ylabel("Porosity",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1].set_ylabel("Overpotential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
    for i in range(0,2):
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath + Target+"Plots/" +
        f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Por Neg_S_eta.png", dpi=dpi) 
    plt.close()  # close the figure to save RAM

    """ update: disable cracking related things just because when sei-on-crack is false it will easily give errors; 
    if model_options.__contains__("SEI on cracks"):
        Num_subplot = 2;
        fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
        axs[0].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Negative electrode roughness ratio"][0],'-o',label="First")
        axs[0].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Negative electrode roughness ratio"][-1],'-^',label="Last"  )
        axs[1].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Total SEI on cracks thickness [m]"][0],'-o',label="First" )
        axs[1].plot(my_dict_AGE["x_n [m]"], my_dict_AGE["CDend Total SEI on cracks thickness [m]"][-1],'-^',label="Last" )
        axs[0].set_title("Tot Neg SEI on cracks thickness",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[1].set_title("Neg roughness ratio",   fontdict={'family':'DejaVu Sans','size':fs+1})
        axs[0].set_ylabel("SEI on cracks thickness [m]",   fontdict={'family':'DejaVu Sans','size':fs})
        axs[1].set_ylabel("Roughness ratio",   fontdict={'family':'DejaVu Sans','size':fs})
        for i in range(0,2):
            axs[i].set_xlabel("Dimensional Neg thickness",   fontdict={'family':'DejaVu Sans','size':fs})
            labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
            axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
            axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
            axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
        plt.savefig(BasicPath + Target+"Plots/" +
            f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Cracks related spatial.png", dpi=dpi) 
        plt.close()  # close the figure to save RAM """
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte concentration [mol.m-3]"][0],'-o',label="First")
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte concentration [mol.m-3]"][-1],'-^',label="Last"  )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte potential [V]"][0],'-o',label="First" )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte potential [V]"][-1],'-^',label="Last" )
    axs[0].set_title("Electrolyte concentration",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Electrolyte potential",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[0].set_ylabel("Concentration [mol.m-3]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1].set_ylabel("Potential [V]",   fontdict={'family':'DejaVu Sans','size':fs})
    for i in range(0,2):
        axs[i].set_xlabel("Dimensional Cell thickness",   fontdict={'family':'DejaVu Sans','size':fs})
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath + Target+"Plots/" +
        f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Electrolyte concentration and potential.png", dpi=dpi)
    plt.close()  # close the figure to save RAM
    Num_subplot = 2;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(8,3.2),tight_layout=True)
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte diffusivity [m2.s-1]"][0],'-o',label="First")
    axs[0].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte diffusivity [m2.s-1]"][-1],'-^',label="Last"  )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte conductivity [S.m-1]"][0],'-o',label="First" )
    axs[1].plot(my_dict_AGE["x [m]"], my_dict_AGE["CDend Electrolyte conductivity [S.m-1]"][-1],'-^',label="Last" )
    axs[0].set_title("Electrolyte diffusivity",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[1].set_title("Electrolyte conductivity",   fontdict={'family':'DejaVu Sans','size':fs+1})
    axs[0].set_ylabel("Diffusivity [m2.s-1]",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[1].set_ylabel("Conductivity [S.m-1]",   fontdict={'family':'DejaVu Sans','size':fs})
    for i in range(0,2):
        axs[i].set_xlabel("Dimensional Cell thickness",   fontdict={'family':'DejaVu Sans','size':fs})
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].ticklabel_format(style='sci', axis='x', scilimits=(-1e-2,1e-2))
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath + Target+"Plots/" +
        f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Electrolyte diffusivity and conductivity.png", dpi=dpi)
    plt.close()  # close the figure to save RAM
    return

def Plot_Dryout(
    Cyc_Update_Index,mdic_dry,ce_EC_0,index_exp,Temper_i,
    Scan_i,Re_No,BasicPath, Target,fs,dpi):
    
    CeEC_All =np.full(np.size(mdic_dry["Ratio_CeEC_All"]),ce_EC_0); 
    for i in range(1,np.size(mdic_dry["Ratio_CeEC_All"])):
        for k in range(0,i):
            CeEC_All[i] *= mdic_dry["Ratio_CeEC_All"][k]
    mdic_dry["CeEC_All"]= CeEC_All.tolist()

    Num_subplot = 3;
    fig, axs = plt.subplots(1,Num_subplot, figsize=(12,3.2),tight_layout=True)
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_EC_consumed_All"],'-o',label="EC consumed")
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_need_All"],'-.',label="Elely needed")
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_add_All"],'-s',label="Elely added")
    axs[0].plot(Cyc_Update_Index, mdic_dry["Vol_Pore_decrease_All"],'--',label="Pore decreased")
    axs[1].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_Tot_All"],     '-o', label="Total electrolyte in cell" )
    axs[1].plot(Cyc_Update_Index, mdic_dry["Vol_Elely_JR_All"],     '--', label="Total electrolyte in JR" )
    axs[1].plot(Cyc_Update_Index, mdic_dry["Vol_Pore_tot_All"],     '-s', label="Total pore in JR" )
    axs[2].plot(Cyc_Update_Index, mdic_dry["Ratio_Dryout_All"],     '-s', label="Dry out ratio" )
    axs[2].set_ylabel("Ratio",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[2].set_xlabel("Cycle number",   fontdict={'family':'DejaVu Sans','size':fs})
    axs[2].set_title("Dry out ratio",   fontdict={'family':'DejaVu Sans','size':fs+1})
    for i in range(0,2):
        axs[i].set_xlabel("Cycle number",   fontdict={'family':'DejaVu Sans','size':fs})
        axs[i].set_ylabel("Volume [mL]",   fontdict={'family':'DejaVu Sans','size':fs})
        axs[i].set_title("Volume",   fontdict={'family':'DejaVu Sans','size':fs+1})
        labels = axs[i].get_xticklabels() + axs[i].get_yticklabels(); [label.set_fontname('DejaVu Sans') for label in labels]
        axs[i].tick_params(labelcolor='k', labelsize=fs, width=1) ;  del labels;
        axs[i].legend(prop={'family':'DejaVu Sans','size':fs-2},loc='best',frameon=False)    
    plt.savefig(BasicPath + Target+"Plots/" +
        f"Scan_{Scan_i}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i- 273.15))}degC Volume_total.png", 
        dpi=dpi)
    plt.close()  # close the figure to save RAM
    return

# update 230312: add a function to get the discharge capacity and resistance
def Get_0p1s_R0(sol_RPT,Index,cap_full):
    Res_0p1s = []; SOC = [100,];
    for i,index in enumerate(Index):
        cycle = sol_RPT.cycles[index]
        Res_0p1s.append(   (
            np.mean(cycle.steps[1]["Terminal voltage [V]"].entries[-10:-1])
            - cycle.steps[2]["Terminal voltage [V]"].entries[0]
        ) / cycle.steps[2]["Current [A]"].entries[0] * 1000)
        if i > 0:
            Dis_Cap = abs(
                cycle.steps[2]["Discharge capacity [A.h]"].entries[0] 
                - cycle.steps[2]["Discharge capacity [A.h]"].entries[-1] )
            SOC.append(SOC[-1]-Dis_Cap/cap_full*100)
    return Res_0p1s[12],Res_0p1s,SOC

# update 230517 add a function to get R_50%SOC from C/2 discharge
def Get_R_from_0P5C_CD(step_0P5C_CD,cap_full):
    # print("Total data points: ",len(step_0P5C_CD["Time [h]"].entries))
    Dis_Cap = abs(
        step_0P5C_CD["Discharge capacity [A.h]"].entries[0] 
        - step_0P5C_CD["Discharge capacity [A.h]"].entries )
    SOC_0p5C = (1-Dis_Cap/cap_full)*100
    #print(SOC_0p5C)
    V_ohmic = (
    step_0P5C_CD['Battery open-circuit voltage [V]'].entries 
    - step_0P5C_CD["Terminal voltage [V]"].entries
    + step_0P5C_CD["Battery particle concentration overpotential [V]"].entries 
    + step_0P5C_CD["X-averaged battery concentration overpotential [V]" ].entries
    )
    # print("Applied current [A]:",step_0P5C_CD["Current [A]"].entries[0])
    Res_0p5C = V_ohmic/step_0P5C_CD["Current [A]"].entries[0] * 1e3
    Res_0p5C_50SOC = np.interp(50,np.flip(SOC_0p5C),np.flip(Res_0p5C),)
    SOC_0p5C = SOC_0p5C.tolist()
    Res_0p5C = Res_0p5C.tolist()
    Res_0p5C_50SOC = Res_0p5C_50SOC.tolist()
    return Res_0p5C_50SOC,Res_0p5C,SOC_0p5C  #,Rohmic_CD_2

# Update 23-05-18
# function to get cell average index from several cells in one T and one Exp
def Get_Cell_Mean_1T_1Exp(Exp_Any_AllData,Exp_temp_i_cell):
    # Get the cell with smallest X "Charge Throughput (A.h)" - to interpolate
    X_1_last = [];    X_5_last = [];
    for cell in Exp_temp_i_cell:
        df = Exp_Any_AllData[cell]["Extract Data"]
        xtemp = np.array(df["Charge Throughput (A.h)"])/1e3
        X_1_last.append(xtemp[-1])
        index_Res = df[df['0.1s Resistance (Ohms)'].le(10)].index
        xtemp2 = np.array(df["Charge Throughput (A.h)"][index_Res])/1e3
        X_5_last.append(xtemp2[-1])
    index_X1 = np.argmin(X_1_last) # find index of the smallest value
    index_X5 = np.argmin(X_5_last)
    # Get the interpolate ones
    Y_1_st = []; Y_2_st = []; Y_3_st = []; Y_4_st = []; Y_5_st = [];Y_6_st = [];
    df_t1 = Exp_Any_AllData[Exp_temp_i_cell[index_X1]]["Extract Data"]
    X_1_st = np.array(df_t1["Charge Throughput (A.h)"])/1e3
    df_t5 = Exp_Any_AllData[Exp_temp_i_cell[index_X5]]["Extract Data"]
    index_Res_5 = df_t5[df_t5['0.1s Resistance (Ohms)'].le(10)].index
    X_5_st = np.array(df_t5["Charge Throughput (A.h)"][index_Res_5])/1e3
    for cell in Exp_temp_i_cell:
        df = Exp_Any_AllData[cell]["Extract Data"]
        chThr_temp = np.array(df["Charge Throughput (A.h)"])/1e3
        df_DMA = Exp_Any_AllData[cell]["DMA"]["LLI_LAM"]
        Y_1_st.append( list(np.interp(X_1_st,chThr_temp,
                        np.array(df_DMA["SoH"])*100)) )
        Y_2_st.append( list(np.interp(X_1_st,chThr_temp,
                        np.array(df_DMA["LLI"])*100)) )
        Y_3_st.append( list(np.interp(X_1_st,chThr_temp,
                        np.array(df_DMA["LAM NE_tot"])*100)) )
        Y_4_st.append( list(np.interp(X_1_st,chThr_temp,
                        np.array(df_DMA["LAM PE"])*100)) )
        Y_6_st.append(list(
            np.interp(X_1_st, chThr_temp, 
            df["Age set average temperature (degC)"].astype(float))))

        index_Res = df[df['0.1s Resistance (Ohms)'].le(10)].index
        Y_5_st.append( list(np.interp(
            X_5_st,np.array(df["Charge Throughput (A.h)"][index_Res])/1e3,
            np.array(df["0.1s Resistance (Ohms)"][index_Res])*1e3, )) )
    # Get the mean values of the interpolate ones:
    def Get_Mean(X_1_st,Y_1_st):
        Y_1_st_avg = []
        for i in range(len(X_1_st)):
            sum=0
            for j in range(len(Y_1_st)):
                sum += Y_1_st[j][i]
            Y_1_st_avg.append(sum/len(Y_1_st))
        return Y_1_st_avg   
    Y_1_st_avg = Get_Mean(X_1_st,Y_1_st)
    Y_2_st_avg = Get_Mean(X_1_st,Y_2_st)
    Y_3_st_avg = Get_Mean(X_1_st,Y_3_st)
    Y_4_st_avg = Get_Mean(X_1_st,Y_4_st)
    Y_5_st_avg = Get_Mean(X_5_st,Y_5_st)
    Y_6_st_avg = Get_Mean(X_1_st,Y_6_st)
    XY_pack = [
        X_1_st,X_5_st,Y_1_st_avg,Y_2_st_avg,
        Y_3_st_avg,Y_4_st_avg,Y_5_st_avg,Y_6_st_avg]
    return XY_pack

# Update 23-05-18 Compare MPE - code created by ChatGPT
def mean_percentage_error(A, B):
    if 0 in A:
        indices = np.where(A == 0)[0]
        A = np.delete(A, indices)
        B = np.delete(B, indices)
        print(A,B)
    errors = np.abs(A - B) / np.abs(A)
    mpe = np.mean(errors) * 100
    return mpe

# Update 23-05-18
# compare modelling result with modelling:
# idea: do interpolation following TODO this is where weighting works
# initial:: X_1_st,X_5_st,Y_1_st_avg,Y_2_st_avg,Y_3_st_avg,Y_4_st_avg,Y_5_st_avg
# to compare: my_dict_RPT
def Compare_Exp_Model(
        my_dict_RPT, XY_pack, Scan_i, Re_No,
        index_exp, Temper_i,BasicPath, Target,fs,dpi, PlotCheck):
    
    [X_1_st,X_5_st,Y_1_st_avg,Y_2_st_avg,
        Y_3_st_avg,Y_4_st_avg,Y_5_st_avg,Y_6_st_avg] = XY_pack
    mX_1 = my_dict_RPT['Throughput capacity [kA.h]']
    if mX_1[-1] > X_1_st[-1]:
        punish = 1; 
        mX_1_st = X_1_st   # do interpolation on modelling result
        mY_1_st = np.interp(mX_1_st,my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT['CDend SOH [%]'])
        mY_2_st = np.interp(mX_1_st,my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LLI [%]"])
        mY_3_st = np.interp(mX_1_st,my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LAM_ne [%]"])
        mY_4_st = np.interp(mX_1_st,my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["CDend LAM_pe [%]"])
        mY_6_st = np.interp(mX_1_st,my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["avg_Age_T"])
        # experiment result remain unchanged
        Y_1_st_avgm = Y_1_st_avg
        Y_2_st_avgm = Y_2_st_avg
        Y_3_st_avgm = Y_3_st_avg
        Y_4_st_avgm = Y_4_st_avg
        Y_6_st_avgm = Y_6_st_avg
    else:                # do interpolation on expeirment results
        punish = X_1_st[-1] / mX_1[-1]  # punishment error, add when simulation end early
        mX_1_st = mX_1 #  standard for experiment following modelling
        Y_1_st_avgm = np.interp(mX_1_st,X_1_st,Y_1_st_avg)
        Y_2_st_avgm = np.interp(mX_1_st,X_1_st,Y_2_st_avg)
        Y_3_st_avgm = np.interp(mX_1_st,X_1_st,Y_3_st_avg)
        Y_4_st_avgm = np.interp(mX_1_st,X_1_st,Y_4_st_avg)
        Y_6_st_avgm = np.interp(mX_1_st,X_1_st,Y_6_st_avg)
        mY_1_st = my_dict_RPT['CDend SOH [%]']
        mY_2_st = my_dict_RPT["CDend LLI [%]"]
        mY_3_st = my_dict_RPT["CDend LAM_ne [%]"]
        mY_4_st = my_dict_RPT["CDend LAM_pe [%]"]
        mY_6_st = my_dict_RPT["avg_Age_T"]
    if mX_1[-1] > X_5_st[-1]:
        mX_5_st = X_5_st   
        mY_5_st = np.interp(mX_5_st,my_dict_RPT['Throughput capacity [kA.h]'], 
            my_dict_RPT["Res_midSOC"])
        Y_5_st_avgm = Y_5_st_avg
    else:
        mX_5_st = mX_1 #  standard for experiment following modelling
        mY_5_st = my_dict_RPT["Res_midSOC"]
        Y_5_st_avgm = np.interp(mX_5_st,X_5_st,Y_5_st_avg)
    # Now we can calculate MPE! mean_percentage_error
    mpe_1 = np.sum(abs(np.array(Y_1_st_avgm)-np.array(mY_1_st)))/len(Y_1_st_avgm) # SOH [%]
    mpe_2 = np.sum(abs(np.array(Y_2_st_avgm)-np.array(mY_2_st)))/len(Y_2_st_avgm) # LLI [%]
    mpe_3 = np.sum(abs(np.array(Y_3_st_avgm)-np.array(mY_3_st)))/len(Y_3_st_avgm) # LAM_ne [%]
    mpe_4 = np.sum(abs(np.array(Y_4_st_avgm)-np.array(mY_4_st)))/len(Y_4_st_avgm) # LAM_pe [%]
    mpe_5 = mean_percentage_error(Y_5_st_avgm, mY_5_st) # Res_midSOC
    mpe_6 = mean_percentage_error(Y_6_st_avgm, mY_6_st) # Age set average temperature (degC)
    # total MPE: TODO this is where weighting works
    # SOH and Resistance are directly measured so give more weight; 
    # DMA result is derived from pOCV and come with certain errors
    mpe_tot = (
        0.55*mpe_1 + 0.1*(mpe_2+mpe_3+mpe_4) 
        + 0.15*mpe_5 + 0.15* mpe_6 +
        (punish>1.0)* punish * 2 )
    # plot and check:
    if PlotCheck == True:
        fig, axs = plt.subplots(5,1, figsize=(6,13),tight_layout=True)
        axs[0].plot(mX_1_st, mY_1_st,'-o', label="Model" )
        axs[0].plot(mX_1_st, Y_1_st_avgm,'-o', label="Exp"  )
        axs[1].plot(mX_1_st, mY_2_st,'-o', label="Model" )
        axs[1].plot(mX_1_st, Y_2_st_avgm,'-o', label="Exp"  )
        axs[2].plot(mX_1_st, mY_3_st,'-o', label="Model" )
        axs[2].plot(mX_1_st, Y_3_st_avgm,'-o', label="Exp"  )
        axs[3].plot(mX_1_st, mY_4_st,'-o', label="Model" )
        axs[3].plot(mX_1_st, Y_4_st_avgm,'-o', label="Exp"  )
        axs[4].plot(mX_5_st, mY_5_st,'-o', label="Model" )
        axs[4].plot(mX_5_st, Y_5_st_avgm,'-o', label="Exp"  )
        axs[0].legend()
        axs[0].set_ylabel("SOH %")
        axs[1].set_ylabel("LLI %")
        axs[2].set_ylabel("LAM NE %")
        axs[3].set_ylabel("LAM PE %")
        axs[4].set_ylabel(r"Lump resistance [m$\Omega$]")
        axs[4].set_xlabel("Charge Throughput (kA.h)")
        fig.suptitle(
            f"Scan {str(Scan_i)}-Exp-{index_exp}-{str(int(Temper_i-273.15))}"
            +r"$^\circ$C - Calculate Error", fontsize=fs+2)
        plt.savefig(BasicPath + Target
            +"Plots/"+ f"Scan {str(Scan_i)}_Re_{Re_No}-Exp-{index_exp}-{str(int(Temper_i-273.15))}degC"
            +r"- Calculate Error.png", dpi=dpi)
        plt.close()  # close the figure to save RAM
    mpe_all = np.around(
        [
            mpe_tot,mpe_1,mpe_2,
            mpe_3,mpe_4,mpe_5,mpe_6,
            punish],2)
    return mpe_all

# read scan files:
# def load_combinations_from_csv(Para_file):
#     dataframe = pd.read_csv(Para_file)
#     parameter_names = dataframe.columns.tolist()
#     combinations = dataframe.values.tolist()
#     # Get all para
#     Para_dict_list = []
#     # get all dictionaries
#     for combination in combinations:
#         input_dict = {}
#         for parameter_name,para_value in zip(parameter_names,combination ):
#             input_dict[parameter_name] = para_value
#         Para_dict_list.append(input_dict)
#     print(f"Total scan case is {len(Para_dict_list)}")

#     return Para_dict_list

def load_combinations_from_csv(Para_file):
    dataframe = pd.read_csv(Para_file)
    # 清理数据，确保格式正确
    dataframe.reset_index(inplace=True)
    dataframe.columns = ['Scan No', 'Exp No.'] + dataframe.columns[2:].tolist()
    
    parameter_names = dataframe.columns.tolist()
    combinations = dataframe.values.tolist()

    # 获取所有参数字典
    Para_dict_list = []
    for combination in combinations:
        input_dict = {}
        for parameter_name, para_value in zip(parameter_names, combination):
            input_dict[parameter_name] = para_value
        Para_dict_list.append(input_dict)

    print(f"总扫描案例数为 {len(Para_dict_list)}")
    return Para_dict_list


def Initialize_exp_text( index_exp, V_max, V_min, Add_Rest):
    """ TODO: if at the end of ageing cycles the cell SOC is NOT 100%SOC, then
    need to top it up before RPT or at the start of RPT! That should involve 
    "0.3C charge to 4.2V" if the SOC after ageing cycles is far away from 100%
    like Exp-1 and 2. Doing a CV hold at 4.2V at 85%SOC can sometimes work but 
    prone to fail
    This applied to Exp-1 and 2, but Exps 3,5-9 are fine!  """
    if index_exp ==1:
        charge_time_mins = 1 * 60 
        exp_AGE_text = [(
            f"Discharge at 1C until {V_min}V", 
            f"Charge at 0.3C for {charge_time_mins} minutes or until {V_max}V",
            ),  ]  # *  setting on cycler is 516, rather than 514 in wiki
        step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
    
    elif index_exp ==2:
        discharge_time_mins = 0.15* 60 * 4.86491/5
        charge_time_mins = 0.5* 60 * 4.86491/5
        exp_AGE_text = [(
            f"Discharge at 1C for {discharge_time_mins} minutes or until {V_min}V", 
            f"Charge at 0.3C for {charge_time_mins} minutes or until {V_max}V",
            ),  ]  # *  setting on cycler is 516, rather than 514 in wiki
        step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
    elif index_exp ==3:
        discharge_time_mins = 0.15* 60 * 4.86491/5
        charge_time_mins = 0.5* 60 * 4.86491/5
        exp_AGE_text = [(
            f"Discharge at 1C for {discharge_time_mins} minutes or until {V_min}V", 
            f"Charge at 0.3C until {V_max}V",
            f"Hold at {V_max} V until C/100",
            ),  ]   # *  setting on cycler is 515, rather than 514 in wiki
        step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
    elif index_exp ==5:
        if Add_Rest:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                "Rest for 1 second", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]  # *  78
            step_AGE_CD =0;   step_AGE_CC =2;   step_AGE_CV =3;
        else:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]  # *  78
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;
    elif index_exp ==6:
        if Add_Rest:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                "Rest for 1 second", 
                f"Charge at 1.2C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =2;   step_AGE_CV =3;
        else:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 1.2C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2;

    elif index_exp ==7:
        if Add_Rest:
            exp_AGE_text = [(
                f"Discharge at 0.5C until {V_min}V", 
                "Rest for 1 second", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =2;   step_AGE_CV =3
        else:
            exp_AGE_text = [(
                f"Discharge at 0.5C until {V_min}V",  
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
    elif index_exp ==8:
        if Add_Rest:
            exp_AGE_text = [(
                f"Discharge at 2C until {V_min}V", 
                "Rest for 1 second", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]
            step_AGE_CD =0;   step_AGE_CC =2;   step_AGE_CV =3
        else:
            exp_AGE_text = [(
                f"Discharge at 2C until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ]
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
    elif index_exp ==9:
        if Add_Rest:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                "Rest for 1 second", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =2;   step_AGE_CV =3
        else:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at 0.3C until {V_max}V",
                f"Hold at {V_max} V until C/100",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2

    elif index_exp ==10: # Added in 241020
        if Add_Rest:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                "Rest for 1 second", 
                f"Charge at C/3 until {V_max}V",
                f"Hold at {V_max} V until C/20",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =2;   step_AGE_CV =3
        else:
            exp_AGE_text = [(
                f"Discharge at 1C until {V_min}V", 
                f"Charge at C/3 until {V_max}V",
                f"Hold at {V_max} V until C/20",
                ),  ] 
            step_AGE_CD =0;   step_AGE_CC =1;   step_AGE_CV =2
    else:
        print("Not yet implemented!")

    # now for RPT: 

    # Added 241020 Xinlei
    exp_RPT_Warwick = [(
        # 1/3C discharge cycle 
        f"Discharge at C/3 until {V_min} V",  
        "Rest for 1 hours",  
        f"Charge at C/3 until {V_max} V",
        f"Hold at {V_max}V until C/20",
        "Rest for 1 hours",
        # 1C discharge capacity 
        f"Discharge at 1C until {V_min} V",  
        "Rest for 1 hours",
        f"Charge at 1C until {V_max} V",
        f"Hold at {V_max}V until C/100",
        "Rest for 1 hours",

        # # 0.5C top up 
        # f"Charge at 0.5C until {V_max} V",
        # f"Hold at {V_max}V until C/100",
        )
    ]

    step_0p1C_CD = 5; # 1C discharge
    step_0p1C_CC = 7; # 1C charge
    step_0p1C_RE =6;     # Rest after 1C discharge
    step_0p5C_CD = 5;  # Meaningless,choose to be 0.1C discharge following 1C discharge

    exp_RPT_Need_TopUp = [ (
        f"Charge at 0.3C until {V_max}V",
        f"Hold at {V_max}V until C/100",
        "Rest for 1 hours",
        "Rest for 10 s",   # add here to correct values of step_0p1C_CD 
        # 0.1C cycle 
        f"Discharge at 0.1C until {V_min} V",  
        "Rest for 3 hours",  
        f"Charge at 0.1C until {V_max} V",
        f"Hold at {V_max}V until C/100",
        "Rest for 1 hours",
        # 0.5C cycle 
        f"Discharge at 0.5C until {V_min} V",  
        "Rest for 3 hours",
        f"Charge at 0.5C until {V_max} V",
        f"Hold at {V_max}V until C/100",
        # Update 23-11-17: add one more 0.5C cycle to increase throughput capacity
        f"Discharge at 0.5C until {V_min} V",  
        "Rest for 3 hours",
        f"Charge at 0.5C until {V_max} V",
        f"Hold at {V_max}V until C/100",   
        "Rest for 3 hours",  
        ) ] 
    exp_RPT_No_TopUp = [ (
        "Rest for 10 s",   # add here to correct values of step_0p1C_CD
        "Rest for 10 s",   # add here to correct values of step_0p1C_CD
        "Rest for 1 hours", 
        "Rest for 10 s",   # add here to correct values of step_0p1C_CD
        # 0.1C cycle 
        f"Discharge at 0.1C until {V_min} V",  
        "Rest for 3 hours",  
        f"Charge at 0.1C until {V_max} V",
        f"Hold at {V_max}V until C/100",
        "Rest for 1 hours",
        # 0.5C cycle 
        f"Discharge at 0.5C until {V_min} V",  
        "Rest for 3 hours",
        f"Charge at 0.5C until {V_max} V",
        f"Hold at {V_max}V until C/100",
        # Update 23-11-17: add one more 0.5C cycle to increase throughput capacity
        f"Discharge at 0.5C until {V_min} V",  
        "Rest for 3 hours",
        f"Charge at 0.5C until {V_max} V",
        f"Hold at {V_max}V until C/100",   
        "Rest for 3 hours",  
        ) ] 
    exp_breakin_text = [ (
        # refill
        #"Rest for 10 s",   # add here to correct values of step_0p1C_CD
        #f"Hold at {V_max}V until C/100",
        #"Rest for 10 s", # Mark Ruihe change ad hoc setting for LFP 

        # Ruihe's old break-in
        # f"Discharge at 0.5C until {V_max-0.2}V", # start from discharge as it is easier for unbalanced cells
        # f"Charge at 0.3C until {V_max}V",
        # f"Hold at {V_max}V until C/100",
        # "Rest for 1 hours", 
        # # 0.1C cycle 
        # f"Discharge at 0.1C until {V_min} V",  
        # "Rest for 3 hours",  
        # f"Charge at 0.1C until {V_max} V",
        # f"Hold at {V_max}V until C/100",
        # "Rest for 1 hours",
        # # 0.5C cycle 
        # f"Discharge at 0.5C until {V_min} V",  
        # "Rest for 3 hours",
        # f"Charge at 0.5C until {V_max} V",
        # f"Hold at {V_max}V until C/100",
        # # Update 23-11-17: add one more 0.5C cycle to increase throughput capacity
        # f"Discharge at 0.5C until {V_min} V",  
        # "Rest for 3 hours",
        # f"Charge at 0.5C until {V_max} V",
        # f"Hold at {V_max}V until C/100",   
        # "Rest for 3 hours",  

        ## Here break-in is set as RPT, 241021
    
        # 1/3C discharge cycle 
        f"Discharge at C/3 until {V_min} V",  
        "Rest for 1 hours",  
        f"Charge at C/3 until {V_max} V",
        f"Hold at {V_max}V until C/20",
        "Rest for 1 hours",
        # 1C discharge capacity 
        f"Discharge at 1C until {V_min} V",  
        "Rest for 1 hours",
        f"Charge at 1C until {V_max} V",
        f"Hold at {V_max}V until C/100",
        "Rest for 1 hours",

        # # 0.5C top up 
        # f"Charge at 0.5C until {V_max} V",
        # f"Hold at {V_max}V until C/100",
        ) ] 
    exp_RPT_GITT_text = [ (
        "Rest for 5 minutes (1 minute period)",  
        "Rest for 1.2 seconds (0.1 second period)",  
        f"Discharge at C/2 for 4.8 minutes or until {V_min}V (0.1 second period)",
        "Rest for 1 hour", # (5 minute period)  
        ) ]
    exp_refill = [ (
        f"Charge at 0.3C until {V_max}V",
        f"Hold at {V_max}V until C/100",
        "Rest for 1 hours", 
        ) ] 
    if index_exp == 1: 
        # adjust to 30% SOC, SOC before this step must be 100%
        charge_time_mins = 1 * 60 
        exp_adjust_before_age = [ (
            f"Discharge at 1C until {V_min} V",
            f"Hold at {V_min}V until C/100",
            "Rest for 4 hours", 
            f"Charge at 0.3C for {charge_time_mins} minutes or until {V_max}V",
            ) ] 
        # Update: 231219 for Exp-2, need to top up the cell after ageing and before C/10
        exp_RPT_text = exp_RPT_Need_TopUp 
    elif index_exp == 2:
        # adjust to 85% SOC, SOC before this step must be 100%
        discharge_time_mins = 0.15* 60 * 4.86491/5
        exp_adjust_before_age = [ (
            f"Discharge at 1C for {discharge_time_mins} minutes or until {V_min}V", # discharge for 15%SOC
            "Rest for 3 hours", 
            ) ] 
        # Update: 231219 for Exp-2, need to top up the cell after ageing and before C/10
        exp_RPT_text = exp_RPT_Need_TopUp 
    elif index_exp ==3:
        exp_adjust_before_age = [ (
            # just a place holder for now TODO 
            "Rest for 1 hours", 
            ) ] 
        exp_RPT_text = exp_RPT_No_TopUp 

    elif index_exp ==10: ## Warick 241021
        exp_adjust_before_age = [ (
            # just a place holder for now TODO 
            "Rest for 1 hours", 
            ) ] 
        exp_RPT_text = exp_RPT_Warwick 

    else:
        exp_adjust_before_age = [ (
            # just a place holder for now TODO 
            "Rest for 1 hours", 
            ) ] 
        exp_RPT_text = exp_RPT_No_TopUp 
    # step index for RPT # revised Xinlei 241021 to match new RPT
    # step_0p1C_CD = 4; step_0p1C_CC = 6;   step_0p1C_RE =5;    
    # step_0p5C_CD = 9;  
    Pack_return = [
        exp_AGE_text, step_AGE_CD, step_AGE_CC, step_AGE_CV,
        exp_breakin_text, exp_RPT_text, exp_RPT_GITT_text, 
        exp_refill,exp_adjust_before_age,
        step_0p1C_CD, step_0p1C_CC, step_0p1C_RE, step_0p5C_CD
        ] 
    return Pack_return

# update 231205: write a function to get tot_cyc,cyc_age,update
def Get_tot_cyc(Runshort,index_exp,Temp_K,Scan_i):
    print(f"Aging T is {Temp_K}degC")
    if Runshort == "GEM-2":
        if index_exp == 1: # 0~30% SOC cycling
            tot_cyc = 3084; cyc_age = 257; update = 257; 
        if index_exp == 2: # 70~85% SOC cycling                 
            tot_cyc = 6192; cyc_age = 516; update = 516;  # actually 516 on cycler
        if index_exp == 3: # 70~85% SOC cycling
            tot_cyc = 6192; cyc_age = 514; update = 514;  # actually 516 on cycler
        if index_exp == 5: # 0~100% SOC cycling
            tot_cyc = 1170; cyc_age = 78; update = 78;    # actually 77 on cycler
        if index_exp == 6:
            tot_cyc = 1170*10; cyc_age = 78; update = 78;
        if index_exp == 7:
            tot_cyc = 1170*10; cyc_age = 78; update = 78;
        if index_exp == 8:
            tot_cyc = 1170*10; cyc_age = 78; update = 78;
        if index_exp == 9:
            tot_cyc = 1170*10; cyc_age = 78; update = 78;
    # Xinlei add 100 cycles    
    elif Runshort == "War":
        if index_exp == 10:
            tot_cyc = 10; cyc_age = 10; update = 10; ## Changed 241021
    elif Runshort == "Reservoir":
        pass

    else:  # really running short    
        if index_exp in list(np.arange(1,10)):
            tot_cyc = 2; cyc_age = 1; update =1           # 2 1 1 
        else:
            print(f"Exp {index_exp} Not yet implemented!")
    return tot_cyc,cyc_age,update

# Get information for Experiment - LG M50 Ageing dataset
def Get_Exp_Pack():
    Exp_Path = [
        "Expt 1 - Si-based Degradation/",
        "Expt 2,2 - C-based Degradation 2/",
        "Expt 3 - Cathode Degradation and Li-Plating/",
        "Expt 4 - Drive Cycle Aging (Control)/",
        "Expt 5 - Standard Cycle Aging (Control)/",]
    Exp_head = [
        "Expt 1",
        "Expt 2,2",
        "Expt 3",
        "Expt 4",
        "Expt 5",]
    Exp_1_Cell = ["A","B","J","D","E","F","K","L","M"];
    Exp_1_Temp = {
        "A":"10","B":"10","J":"10",
        "D":"25","E":"25","F":"25",
        "K":"40","L":"40","M":"40",}
    Temp_Cell_Exp_1 = {
        "10":["A","B","J"],
        "25":["D","E","F"],
        "40":["K","L","M"],}
    Exp_2_Cell = ["A","B","C","D","E","F"];
    Exp_2_Temp = {
        "A":"10","B":"10",
        "C":"25","D":"25",
        "E":"40","F":"40",}
    Temp_Cell_Exp_2 = {
        "10":["A","B"],
        "25":["C","D"],
        "40":["E","F"],}
    Exp_3_Cell = ["A","B","C","D","E","F","G","H","I"];
    Exp_3_Temp = {
        "A":"10","B":"10","C":"10",
        "D":"25","E":"25","F":"25",
        "G":"40","H":"40","I":"40"}
    Temp_Cell_Exp_3 = {
        "10":["A","B","C"],
        "25":["D","E","F"],
        "40":["G","H","I"],}
    Exp_4_Cell = ["A","B","C","D","E","F","G","H"];
    Exp_4_Temp = {
        "A":"10","B":"10","C":"10",
        "D":"25","E":"25",
        "F":"40","G":"40","H":"40",}
    Temp_Cell_Exp_4 = {
        "10":["A","B","C"],
        "25":["D","E",],
        "40":["F","G","H"],}
    Exp_5_Cell = ["A","B","C","D","E","F","G","H"];
    Exp_5_Temp = {
        "A":"10","B":"10","C":"10",
        "D":"25","E":"25",
        "F":"40","G":"40","H":"40",}
    Temp_Cell_Exp_5 = {
        "10":["A","B","C"],
        "25":["D","E",],
        "40":["F","G","H"],}
    Exp_All_Cell  = [Exp_1_Cell,Exp_2_Cell,Exp_3_Cell,Exp_4_Cell,Exp_5_Cell]
    Exp_Temp_Cell = [Exp_1_Temp,Exp_2_Temp,Exp_3_Temp,Exp_4_Temp,Exp_5_Temp]
    Temp_Cell_Exp_All = [
        Temp_Cell_Exp_1,Temp_Cell_Exp_2,Temp_Cell_Exp_3,
        Temp_Cell_Exp_4,Temp_Cell_Exp_5]
    Mark_Cell_All = [
        {
        "A":"o","B":">","J":"v",
        "D":"o","E":">","F":"v",
        "K":"o","L":">","M":"v",},
        {
        "A":"10","B":"10",
        "C":"25","D":"25",
        "E":"40","F":"40",},
        {
        "A":"o","B":">","C":"v",
        "D":"o","E":">","F":"v",
        "G":"o","H":">","I":"v",},
        {
        "A":"o","B":">","C":"v",
        "D":"o","E":">",
        "F":"o","G":">","H":"v",},
        {
        "A":"o","B":">","C":"v",
        "D":"o","E":">",
        "F":"o","G":">","H":"v",}]
    Color_Cell_All = [
        {
        "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
        "J":[2/255, 3/255, 226/255,0.7],
        "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],"F":[0, 0, 0,0.7],
        "K":[1,0,0,0.4],"L":[1,0,0,0.4],"M":[1,0,0,0.4],},
        {
        "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
        "D":[0, 0, 0,0.7],"C":[0, 0, 0,0.7],
        "E":[1,0,0,0.4],"F":[1,0,0,0.4],},
        {
        "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
        "C":[2/255, 3/255, 226/255,0.7],
        "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],"F":[0, 0, 0,0.7],
        "G":[1,0,0,0.4],"H":[1,0,0,0.4],"I":[1,0,0,0.4],},
        {
        "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
        "C":[2/255, 3/255, 226/255,0.7],
        "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],
        "F":[1,0,0,0.4],"G":[1,0,0,0.4],"H":[1,0,0,0.4],},
        {
        "A":[2/255, 3/255, 226/255,0.7],"B":[2/255, 3/255, 226/255,0.7],
        "C":[2/255, 3/255, 226/255,0.7],
        "D":[0, 0, 0,0.7],"E":[0, 0, 0,0.7],
        "F":[1,0,0,0.4],"G":[1,0,0,0.4],"H":[1,0,0,0.4],}]
    Exp_pack = [
        Exp_All_Cell,Temp_Cell_Exp_All,
        Exp_Path,Exp_head,Exp_Temp_Cell] 
    return Exp_pack
# Get keys for output
def Get_Output_Keys(Para_dict_i):
    model_options = Para_dict_i["Model option"] 
    keys_loc_RPT = [ # MAY WANT TO SELECT AGEING CYCLE later
        # Default output:
        "x [m]",
        "x_n [m]",
        "x_s [m]",
        "x_p [m]",
        # default: end; 
        "CCend Porosity",
        "CCend Negative electrode interfacial current density [A.m-2]",
        "CCend Electrolyte potential [V]",
        "CCend Electrolyte concentration [mol.m-3]",
        "CCend Negative electrode reaction overpotential [V]",
        "CCend Negative particle surface concentration [mol.m-3]",
        #"CCend Negative electrode roughness ratio",
        #"CCend Total SEI on cracks thickness [m]",

        "CDend Porosity",
        "CDend Negative electrode interfacial current density [A.m-2]",
        "CDend Electrolyte potential [V]",
        "CDend Electrolyte concentration [mol.m-3]",
        "CDend Negative electrode reaction overpotential [V]",
        "CDend Negative particle surface concentration [mol.m-3]",
        #"CDend Negative electrode roughness ratio",
        #"CDend Total SEI on cracks thickness [m]",
        #"REend Total SEI on cracks thickness [m]",
    ]
    keys_tim_RPT = [
        # default: CD
        "CD Time [h]",
        "CD Terminal voltage [V]",
        "CD Anode potential [V]",    # self defined
        "CD Cathode potential [V]",  # self defined
        "CC Time [h]",
        "CC Terminal voltage [V]",
        "CC Anode potential [V]",    # self defined
        "CC Cathode potential [V]",  # self defined
    ]
    keys_cyc_RPT = [   # default: CDend
        "Discharge capacity [A.h]",
        "Throughput capacity [A.h]",
        "CDend Total lithium capacity in particles [A.h]",
        "CDend Loss of capacity to negative lithium plating [A.h]", ## revised negative 241021
        "CDend Loss of capacity to negative SEI [A.h]",## revised negative
        "CDend Loss of capacity to negative SEI on cracks [A.h]",## revised negative
        #"CDend X-averaged total SEI on cracks thickness [m]",
        #"CDend X-averaged negative electrode roughness ratio",
        "CDend Local ECM resistance [Ohm]",
        "CDsta Negative electrode stoichiometry", 
        "CDend Negative electrode stoichiometry",
        "CDsta Positive electrode stoichiometry", 
        "CDend Positive electrode stoichiometry",
        "CDend Negative electrode capacity [A.h]",
        "CDend Positive electrode capacity [A.h]",
    ]

    keys_loc_AGE = [ # MAY WANT TO SELECT AGEING CYCLE later
        # Default output:
        "x [m]",
        "x_n [m]",
        "x_s [m]",
        "x_p [m]",
        # default: end; 
        "CCend Porosity",
        "CCend Negative electrode interfacial current density [A.m-2]",
        "CCend Electrolyte potential [V]",
        "CCend Electrolyte concentration [mol.m-3]",
        "CCend Negative electrode reaction overpotential [V]",
        "CCend Negative particle surface concentration [mol.m-3]",
        "CCend Negative electrode surface potential difference [V]",
        "CCend Negative electrode SEI film overpotential [V]", ## revised negative 241021
        #"CCend Negative electrode roughness ratio",
        #"CCend Total SEI on cracks thickness [m]",

        "CDend Porosity",
        "CDend Negative electrode interfacial current density [A.m-2]",
        "CDend Electrolyte potential [V]",
        "CDend Electrolyte concentration [mol.m-3]",
        "CDend Negative electrode reaction overpotential [V]",
        "CDend Negative particle surface concentration [mol.m-3]",
        #"CDend Negative electrode roughness ratio",
        #"CDend Total SEI on cracks thickness [m]",
        "CDend Negative electrode surface potential difference [V]",
        "CDend Negative electrode SEI film overpotential [V]",
        "CDend Electrolyte diffusivity [m2.s-1]",
        "CDend Electrolyte conductivity [S.m-1]",
    ]
    keys_tim_AGE = [
        # default: CD
        "CD Time [h]",
        "CD Terminal voltage [V]",
        "CD Anode potential [V]",    # self defined
        "CD Cathode potential [V]",  # self defined
        
        "CC Time [h]",
        "CC Terminal voltage [V]",
        "CC Anode potential [V]",    # self defined
        "CC Cathode potential [V]",  # self defined
    ]
    keys_cyc_AGE = [];
    keys_all_RPT = [keys_loc_RPT,keys_tim_RPT,keys_cyc_RPT]
    keys_all_AGE = [keys_loc_AGE,keys_tim_AGE,keys_cyc_AGE]
    keys_all = [keys_all_RPT,keys_all_AGE]
    return keys_all

# Update 240429: Define a new dictionary to write summary result into Excel file
def Write_Dict_to_Excel(
        Flag_Breakin,Para_dict_i,cap_0,
        BasicPath,Target,book_name_xlsx,
        index_exp,Pass_Fail,
        mpe_all,
        model_options,my_dict_RPT,mdic_dry,
        DryOut,Scan_i,str_exp_AGE_text,
        str_exp_RPT_text,str_error_AGE_final,
        str_error_RPT,log_messages,
        ): 
    # *mpe_all = [mpe_tot,mpe_1,mpe_2,mpe_3,mpe_4,mpe_5,mpe_6,punish]
    Dict_Excel = {}
    # 1 - start from basic summary information: pre
    keys_pre = [
        "Scan No","Exp No.","Y or N",
        "Error tot %","Error SOH %","Error LLI %",
        "Error LAM NE %","Error LAM PE %",
        "Error Res %","Error ageT %","Punish",
        "Dry out",]
    value_Pre = [
        Scan_i,index_exp,Pass_Fail,
        *mpe_all,DryOut]
    for key, value in zip(keys_pre,value_Pre):
        Dict_Excel[key] = value
    # 2 - middle is input parameters
    Dict_Excel.update(Para_dict_i)  # include "Total ageing cycles", "Ageing cycles between RPT", 
    # 3 - Finally post-processing data, may be "nan"
    Dict_Excel["exp_AGE_text"] =  str_exp_AGE_text
    Dict_Excel["exp_RPT_text"] =  str_exp_RPT_text
    if Flag_Breakin: # break in cycle succeed:
        tot_char_throughput = my_dict_RPT['Throughput capacity [kA.h]'][-1]
        Dict_Excel["Throughput capacity [kA.h]"] =  tot_char_throughput
        Dict_Excel["CDend SOH [%]"] =  my_dict_RPT['CDend SOH [%]'][-1]
        Porosity = my_dict_RPT["CDend Porosity"][-1]
        x_n = my_dict_RPT["x_n [m]"]
        Dict_Excel["CDend Neg Porosity Sep side"] =  Porosity[len(x_n)-1]
        Dict_Excel["CDend LLI [%]"] =  my_dict_RPT["CDend LLI [%]"][-1]

        if model_options.__contains__("SEI on cracks"):
            Dict_Excel["LLI to sei-on-cracks [%]"] = my_dict_RPT["CDend LLI SEI on cracks [%]"][-1]
        else:
            Dict_Excel["LLI to sei-on-cracks [%]"] = 0.0
        if model_options.__contains__("lithium plating"):
            Dict_Excel["LLI to LiP [%]"] = my_dict_RPT["CDend LLI lithium plating [%]"][-1]
        else:
            Dict_Excel["LLI to LiP [%]"] = 0.0
        if model_options.__contains__("SEI"):
            Dict_Excel["LLI to SEI [%]"] = my_dict_RPT["CDend LLI SEI [%]"][-1]
        else:
            Dict_Excel["LLI to SEI [%]"] = 0.0
        Dict_Excel["CDend LLI due to LAM [%]"] = my_dict_RPT["CDend LLI due to LAM [%]"][-1]
        Dict_Excel["LAM_to_Crack_NE [%]"] = my_dict_RPT["LAM_to_Crack_NE [%] end"]
        Dict_Excel["LAM_to_Crack_PE [%]"] = my_dict_RPT["LAM_to_Crack_PE [%] end"]
        Dict_Excel["LAM_to_Dry [%]"] = my_dict_RPT["LAM_to_Dry [%] end"]
        Dict_Excel["CDend LAM_ne_tot [%]"] = my_dict_RPT["CDend LAM_ne [%]"][-1]
        Dict_Excel["CDend LAM_pe_tot [%]"] = my_dict_RPT["CDend LAM_pe [%]"][-1]
        
        # per kAh:
        Dict_Excel["Cap per kAh"]=(1-Dict_Excel["CDend SOH [%]"])*cap_0/tot_char_throughput
        Dict_Excel["LLI% per kAh"]=Dict_Excel["CDend LLI [%]"]/tot_char_throughput
        Dict_Excel["LAM_ne% per kAh"]=Dict_Excel["CDend LAM_ne_tot [%]"]/tot_char_throughput
        Dict_Excel["LAM_pe% per kAh"]=Dict_Excel["CDend LAM_pe_tot [%]"]/tot_char_throughput
        if DryOut == "On":
            Dict_Excel["Vol_Elely_Tot_All_final"] = mdic_dry["Vol_Elely_Tot_All"][-1]
            Dict_Excel["Vol_Elely_JR_All_final"]  = mdic_dry["Vol_Elely_JR_All"][-1]
            Dict_Excel["Width_all_final"]         = mdic_dry["Width_all"][-1]
        else:
            Dict_Excel["Vol_Elely_Tot_All_final"] = np.nan
            Dict_Excel["Vol_Elely_JR_All_final"]  = np.nan
            Dict_Excel["Width_all_final"]         = np.nan
    else: # break in cycle fail:
        Keys_pos = [
            "Throughput capacity [kA.h]","CDend SOH [%]",
            "CDend Neg Porosity Sep side",
            "CDend LLI [%]","LLI to sei-on-cracks [%]",
            "LLI to LiP [%]","LLI to SEI [%]",
            "CDend LLI due to LAM [%]","LAM_to_Crack_NE [%]",
            "LAM_to_Crack_PE [%]","LAM_to_Dry [%]",
            "CDend LAM_ne_tot [%]","CDend LAM_pe_tot [%]",
            "Cap per kAh","LLI% per kAh","LAM_ne% per kAh",
            "LAM_pe% per kAh","Vol_Elely_Tot_All_final",
            "Vol_Elely_JR_All_final","Width_all_final",
        ]
        for key in Keys_pos:
            Dict_Excel[key] = np.nan

    Dict_Excel["Error_AGE"] = str_error_AGE_final
    Dict_Excel["Error_RPT"] = str_error_RPT
    Dict_Excel["Log"] = log_messages
    # df_excel = pd.DataFrame(Dict_Excel)
    df_excel = pd.DataFrame([Dict_Excel])
    book_name_xlsx_seperate =   str(Scan_i)+ '_' + book_name_xlsx
    df_excel.to_excel(
        BasicPath + Target + "Excel/" 
        +book_name_xlsx_seperate,
        index=False, engine='openpyxl')
    return df_excel


def Run_P2_Excel(
    Para_dict_i,  Path_List,  Re_No,        
    Timelimit,    Options): 
    ##########################################################
    ##############    Part-0: Log of the scripts    ##########
    ##########################################################
    # add 221205: if Timeout=='True', use Patrick's version, disable pool
    #             else, use pool to accelerate 
    # add Return_Sol, on HPC, always set to False, as it is useless, 
    # add 230221: do sol_new['Throughput capacity [A.h]'].entries += sol_old['Throughput capacity [A.h]'].entries 
    #             and for "Throughput energy [W.h]", when use Model.set_initial_conditions_from
    #             this is done inside the two functions Run_Model_Base_On_Last_Solution(_RPT)
    # change 230621: add ability to scan one parameter set at different temperature 
    #                and at Exp-2,3,5
    # 240429: tidy up the script
    # idea to tidy up: create a class: 

    ##########################################################
    ##############    Part-1: Initialization    ##############
    ##########################################################
    # Unpack options:
    [ 
        On_HPC,Runshort,Add_Rest,
        Plot_Exp,Timeout,Return_Sol,
        Check_Small_Time,R_from_GITT,
        dpi,fs] = Options
    ModelTimer = pb.Timer() # start counting time
    if Check_Small_Time == True:
        SmallTimer = pb.Timer()
    import io
    log_buffer = io.StringIO()

    
    [BasicPath, Path_to_ExpData,Target,purpose] = Path_List
    # Gey output keys:
    keys_all = Get_Output_Keys(Para_dict_i)

    # Create folders
    if not os.path.exists(BasicPath +Target):
        os.mkdir(BasicPath +Target)
    if not os.path.exists(BasicPath +Target+"Mats"):
        os.mkdir(BasicPath +Target +"Mats");
    if not os.path.exists(BasicPath +Target+"Plots"):
        os.mkdir(BasicPath +Target+"Plots");
    if not os.path.exists(BasicPath +Target+"Excel"):
        os.mkdir(BasicPath +Target+"Excel");
    
    
    font = {'family' : 'DejaVu Sans','size'   : fs}
    mpl.rc('font', **font)

    # new to define: exp_text_list,  Path_pack
    # exp_index_pack , 
    # define here:
    index_i   = Para_dict_i["Scan No"]  
    Scan_i = int(index_i)
    print(f'Start Now! Scan {Scan_i} Re {Re_No}')  
    index_exp = Para_dict_i["Exp No."] # index for experiment set, can now go for 2,3,5
    print("Exp:", index_exp)
    Temp_K = Para_dict_i["Ageing temperature"]  
    print("T:", Temp_K)
    Round_No = f"Case_{Scan_i}_Exp_{index_exp}_{Temp_K}oC"  # index to identify different rounds of running 
    # Load Niall's data
    [
        Exp_All_Cell,Temp_Cell_Exp_All,
        Exp_Path,Exp_head,Exp_Temp_Cell
        ]  = Get_Exp_Pack()
    book_name_xlsx = f'Re_{Re_No}_{purpose}.xlsx'
    # index_exp should be 1~5 to really have experimental data
    if index_exp in list(np.arange(1,6)) and int(Temp_K) in [10,25,40]:
        Temp_Cell_Exp = Temp_Cell_Exp_All[index_exp-1] 
        Exp_Any_AllData = Read_Exp(
            Path_to_ExpData,Exp_All_Cell[index_exp-1],
            Exp_Path,Exp_head,Exp_Temp_Cell[index_exp-1],
            index_exp-1)
    else:
        Temp_Cell_Exp = "nan"
        Exp_Any_AllData = "nan"
    # update 231205: write a function to get tot_cyc,cyc_age,update
    tot_cyc,cyc_age,update = Get_tot_cyc(Runshort,index_exp,Temp_K,Scan_i)
    Para_dict_i["Total ageing cycles"]       = int(tot_cyc)
    Para_dict_i["Ageing cycles between RPT"] = int(cyc_age)
    Para_dict_i["Update cycles for ageing"]  = int(update) # keys
    
    # set up experiment
    # update 231205: get a new function to initialize exp_AGE_text and exp_RPT_text
    if Para_dict_i["Para_Set"] == "OKane2023_Sunil":
        V_max = 3.65
        cap_0 = 72  #  Mark: change to 72 for LFP  
    else:
        V_max = 4.2     
        cap_0 = 4.86491   # set initial capacity to standardize SOH and get initial SEI thickness  
    V_min = 2.5; 
    [
        exp_AGE_text, step_AGE_CD, step_AGE_CC, step_AGE_CV,
        exp_breakin_text, exp_RPT_text, exp_RPT_GITT_text, 
        exp_refill,exp_adjust_before_age,
        step_0p1C_CD, step_0p1C_CC, step_0p1C_RE, step_0p5C_CD
        ] = Initialize_exp_text(
        index_exp, V_max, V_min, Add_Rest)
    cycle_no = -1; 


    Sol_RPT = [];  Sol_AGE = [];   
    
    # pb.set_logging_level('INFO') # show more information!
    # set_start_method('fork') # from Patrick

    # Un-pack data:
    CyclePack,Para_0 = Para_init(Para_dict_i,4.86491) # initialize the parameter
    [
        Total_Cycles,Cycle_bt_RPT,Update_Cycles,
        RPT_Cycles,Temper_i,Temper_RPT,mesh_list,
        submesh_strech,model_options,
        cap_increase] = CyclePack
    [keys_all_RPT,keys_all_AGE] = keys_all
    str_exp_AGE_text  = str(exp_AGE_text)
    str_exp_RPT_text  = str(exp_RPT_text)
    str_exp_RPT_GITT_text  = str(exp_RPT_GITT_text)

    # define experiment
    Experiment_Long   = pb.Experiment( exp_AGE_text * Update_Cycles  )  
    # update 24-04-2023: delete GITT
    # Update 01-11-2023 add GITT back but with an option 
    # update 231210: refine experiment to avoid charge or hold at 4.2V
    if R_from_GITT: 
        Experiment_Breakin= pb.Experiment( 
            exp_breakin_text * 1
            + exp_RPT_GITT_text*24  + exp_refill * 1    # only do refil if have GITT
            + exp_adjust_before_age*1) 
        Experiment_RPT    = pb.Experiment( 
            exp_RPT_text * 1
            + exp_RPT_GITT_text*24  + exp_refill * 1    # only do refil if have GITT
            + exp_adjust_before_age*1) 
        Cyc_Index_Res = np.arange(1,25,1) 
    else:   # then get resistance from C/2
        Experiment_Breakin= pb.Experiment( 
            exp_breakin_text * 1
            + exp_adjust_before_age*1) 
        Experiment_RPT    = pb.Experiment( 
            exp_RPT_text * 1
            + exp_adjust_before_age*1) 
    

    #####  index definition ######################
    Small_Loop =  int(Cycle_bt_RPT/Update_Cycles);   
    SaveTimes = int(Total_Cycles/Cycle_bt_RPT);   

    # initialize my_dict for outputs
    my_dict_RPT = {}
    for keys in keys_all_RPT:
        for key in keys:
            my_dict_RPT[key]=[];
    my_dict_AGE = {}; 
    for keys in keys_all_AGE:
        for key in keys:
            my_dict_AGE[key]=[]
    my_dict_RPT["Cycle_RPT"] = []
    my_dict_RPT["Res_full"] = []
    my_dict_RPT["Res_midSOC"] = []
    my_dict_RPT["SOC_Res"] = []
    my_dict_AGE["Cycle_AGE"] = []
    my_dict_RPT["avg_Age_T"] = [] # Update add 230617 
    Cyc_Update_Index     =[]

    # update 220924: merge DryOut and Int_ElelyExces_Ratio
    temp_Int_ElelyExces_Ratio =  Para_0["Initial electrolyte excessive amount ratio"] 
    ce_EC_0 = Para_0['EC initial concentration in electrolyte [mol.m-3]'] # used to calculate ce_EC_All
    if temp_Int_ElelyExces_Ratio < 1:
        Int_ElelyExces_Ratio = -1;
        DryOut = "Off";
    else:
        Int_ElelyExces_Ratio = temp_Int_ElelyExces_Ratio;
        DryOut = "On";
    print(f"Scan {Scan_i} Re {Re_No}: DryOut = {DryOut}")
    if DryOut == "On":  
        mdic_dry,Para_0 = Initialize_mdic_dry(Para_0,Int_ElelyExces_Ratio)
    else:
        mdic_dry ={}
    if Check_Small_Time == True:
        print(f'Scan {Scan_i} Re {Re_No}: Spent {SmallTimer.time()} on Initialization')
        SmallTimer.reset()
    ##########################################################
    ##############    Part-2: Run model         ##############
    ##########################################################
    ##########################################################
    Timeout_text = 'I timed out'
    ##########    2-1: Define model and run break-in cycle
    try:  
        # Timelimit = int(3600*2)
        # the following turns on for HPC only!
        if Timeout == True:
            timeout_RPT = TimeoutFunc(
                Run_Breakin, 
                timeout=Timelimit, 
                timeout_val=Timeout_text)
            Result_list_breakin  = timeout_RPT(
                model_options, Experiment_Breakin, 
                Para_0, mesh_list, submesh_strech,
                cap_increase)
        else:
            Result_list_breakin  = Run_Breakin(
                model_options, Experiment_Breakin, 
                Para_0, mesh_list, submesh_strech,
                cap_increase)
        [Model_0,Sol_0,Call_Breakin] = Result_list_breakin
        if Return_Sol == True:
            Sol_RPT.append(Sol_0)
        if Call_Breakin.success == False:
            print("Fail due to Experiment error or infeasible")
            1/0
        if Sol_0 == Timeout_text: # to do: distinguish different failure cases
            print("Fail due to Timeout")
            1/0
        if Sol_0 == "Model error or solver error":
            print("Fail due to Model error or solver error")
            1/0
    except ZeroDivisionError as e:
        str_error_Breakin = str(e)
        if Check_Small_Time == True:
            str_error_Breakin = f"Scan {Scan_i} Re {Re_No}: Fail break-in cycle within {SmallTimer.time()}, need to exit the whole scan now due to {str_error_Breakin} but do not know how!"
            print(str_error_Breakin)
            SmallTimer.reset()
        else:
            str_error_Breakin = f"Scan {Scan_i} Re {Re_No}: Fail break-in cycle, need to exit the whole scan now due to {str_error_Breakin} but do not know how!"
            print(str_error_Breakin)
        Flag_Breakin = False 
    else:
        if Check_Small_Time == True:    
            print(f"Scan {Scan_i} Re {Re_No}: Finish break-in cycle within {SmallTimer.time()}")
            SmallTimer.reset()
        else:
            print(f"Scan {Scan_i} Re {Re_No}: Finish break-in cycle")
        # post-process for break-in cycle - 0.1C only
        my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_0, 
            0, step_0p1C_CD, step_0p1C_CC,step_0p1C_RE , step_AGE_CV   )
        # update 230517 - Get R from C/2 discharge only, discard GITT

        # cap_full = 5; 
        # if R_from_GITT: 
        #     Res_midSOC,Res_full,SOC_Res = Get_0p1s_R0(Sol_0,Cyc_Index_Res,cap_full)
        # else: 
        #     step_0P5C_CD = Sol_0.cycles[0].steps[step_0p5C_CD]
        #     Res_midSOC,Res_full,SOC_Res = Get_R_from_0P5C_CD(step_0P5C_CD,cap_full)
        # my_dict_RPT["SOC_Res"].append(SOC_Res)
        # my_dict_RPT["Res_full"].append(Res_full)
        # my_dict_RPT["Res_midSOC"].append(Res_midSOC)    

        my_dict_RPT["avg_Age_T"].append(Temper_i-273.15)  # Update add 230617              
        #del SOC_Res,Res_full,Res_midSOC

        cycle_count =0
        my_dict_RPT["Cycle_RPT"].append(cycle_count)
        Cyc_Update_Index.append(cycle_count)
        Flag_Breakin = True
        if Check_Small_Time == True:    
            print(f"Scan {Scan_i} Re {Re_No}: Finish post-process for break-in cycle within {SmallTimer.time()}")
            SmallTimer.reset()
        else:
            print(f"Scan {Scan_i} Re {Re_No}: Finish post-process for break-in cycle")
        
    Flag_AGE = True; Flag_partial_AGE = False
    str_error_AGE_final = "Empty";   str_error_RPT = "Empty"; 
    DeBug_List_RPT = "Break in fail"; DeBug_List_AGE = "Break in fail"
    #############################################################
    #######   2-2: Write a big loop to finish the long experiment    
    if Flag_Breakin == True: 
        k=0
        # Para_All.append(Para_0);Model_All.append(Model_0);Sol_All_i.append(Sol_0); 
        Para_0_Dry_old = Para_0;     Model_Dry_old = Model_0  ; Sol_Dry_old = Sol_0;   del Model_0,Sol_0
        while k < SaveTimes:    
            i=0  
            avg_Age_T = []  
            while i < Small_Loop:
                if DryOut == "On":
                    Data_Pack,Paraupdate   = Cal_new_con_Update (  Sol_Dry_old,   Para_0_Dry_old )
                if DryOut == "Off":
                    Paraupdate = Para_0
                # Run aging cycle:
                try:
                    #Timelimit = int(60*60*2)
                    if Timeout == True:
                        timeout_AGE = TimeoutFunc(
                            Run_Model_Base_On_Last_Solution, 
                            timeout=Timelimit, 
                            timeout_val=Timeout_text)
                        Result_list_AGE = timeout_AGE( 
                            Model_Dry_old  , Sol_Dry_old , Paraupdate ,Experiment_Long, 
                            Update_Cycles,Temper_i,mesh_list,submesh_strech )
                    else:
                        Result_list_AGE = Run_Model_Base_On_Last_Solution( 
                            Model_Dry_old  , Sol_Dry_old , Paraupdate ,Experiment_Long, 
                            Update_Cycles,Temper_i,mesh_list,submesh_strech )
                    [Model_Dry_i, Sol_Dry_i , Call_Age,DeBug_List_AGE ] = Result_list_AGE
                    
                    if Return_Sol == True:
                        Sol_AGE.append(Sol_Dry_i)
                    if "Partially" in DeBug_List_AGE[-1]:
                        Flag_partial_AGE = True
                        succeed_cycs = len(Sol_Dry_i.cycles) 
                        if succeed_cycs < Update_Cycles:
                            print(f"Instead of {Update_Cycles}, succeed only {succeed_cycs} cycles")
                            Flag_partial_AGE = True
                    else:
                        if Call_Age.success == False:
                            print("Fail due to Experiment error or infeasible")
                            str_error_AGE = "Experiment error or infeasible"
                            1/0
                        elif Sol_Dry_i == "Model error or solver error":
                            print("Fail due to Model error or solver error")
                            str_error_AGE = "Model error or solver error"
                            1/0
                        else:
                            pass
                    if Sol_Dry_i == Timeout_text: # fail due to timeout
                        print("Fail due to Timeout")
                        str_error_AGE = "Timeout"
                        1/0
                except ZeroDivisionError as e: # ageing cycle fails
                    if Check_Small_Time == True:    
                        str_error_AGE_final = f"Scan {Scan_i} Re {Re_No}: Fail during No.{Cyc_Update_Index[-1]} ageing cycles within {SmallTimer.time()} due to {str_error_AGE}"
                        print(str_error_AGE_final)
                        SmallTimer.reset()
                    else:
                        str_error_AGE_final = f"Scan {Scan_i} Re {Re_No}: Fail during No.{Cyc_Update_Index[-1]} ageing cycles due to {str_error_AGE}"
                        print(str_error_AGE_final)
                    Flag_AGE = False
                    break
                else:                           # ageing cycle SUCCEED
                    succeed_cycs = len(Sol_Dry_i.cycles) 
                    Para_0_Dry_old = Paraupdate; Model_Dry_old = Model_Dry_i; Sol_Dry_old = Sol_Dry_i;   
                    del Paraupdate,Model_Dry_i,Sol_Dry_i
                    
                    if Check_Small_Time == True:    
                        print(f"Scan {Scan_i} Re {Re_No}: Finish for No.{Cyc_Update_Index[-1]} ageing cycles within {SmallTimer.time()}")
                        SmallTimer.reset()
                    else:
                        print(f"Scan {Scan_i} Re {Re_No}: Finish for No.{Cyc_Update_Index[-1]} ageing cycles")

                    # post-process for first ageing cycle and every -1 ageing cycle
                    if k==0 and i==0:    
                        my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                            0, step_AGE_CD , step_AGE_CC , step_0p1C_RE, step_AGE_CV   )     
                        my_dict_AGE["Cycle_AGE"].append(1)
                    # update 240111
                    if Flag_partial_AGE == True:
                        try:
                            my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                                cycle_no, step_AGE_CD , step_AGE_CC , step_0p1C_RE, step_AGE_CV   )    
                        except IndexError:
                            print("The last cycle is incomplete, try [-2] cycle")
                            try:
                                my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                                    -2, step_AGE_CD , step_AGE_CC , step_0p1C_RE, step_AGE_CV   )   
                            except IndexError:
                                print("[-2] cycle also does not work, try first one")
                                try:
                                    my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                                        0, step_AGE_CD , step_AGE_CC , step_0p1C_RE, step_AGE_CV   )  
                                except:
                                    print("Still does not work, less than one cycle, we are in trouble")
                    else:
                        try:
                            my_dict_AGE = GetSol_dict (my_dict_AGE,keys_all_AGE, Sol_Dry_old, 
                                cycle_no, step_AGE_CD , step_AGE_CC , step_0p1C_RE, step_AGE_CV   )    
                        except:
                            print("GetSol_dict fail for a complete ageing set for unknown reasons!!!")
                    cycle_count +=  succeed_cycs 
                    avg_Age_T.append(np.mean(
                        Sol_Dry_old["Volume-averaged cell temperature [C]"].entries))
                    my_dict_AGE["Cycle_AGE"].append(cycle_count)           
                    Cyc_Update_Index.append(cycle_count)
                    
                    if DryOut == "On":
                        mdic_dry = Update_mdic_dry(Data_Pack,mdic_dry)
                    
                    if Check_Small_Time == True:    
                        print(f"Scan {Scan_i} Re {Re_No}: Finish post-process for No.{Cyc_Update_Index[-1]} ageing cycles within {SmallTimer.time()}")
                        SmallTimer.reset()
                    else:
                        pass
                    i += 1;   ##################### Finish small loop and add 1 to i 
            
            # run RPT, and also update parameters (otherwise will have problems)
            if DryOut == "On":
                Data_Pack , Paraupdate  = Cal_new_con_Update (  
                    Sol_Dry_old,   Para_0_Dry_old   )
            if DryOut == "Off":
                Paraupdate = Para_0     
            try:
                # Timelimit = int(60*60*2)
                if Timeout == True:
                    timeout_RPT = TimeoutFunc(
                        Run_Model_Base_On_Last_Solution_RPT, 
                        timeout=Timelimit, 
                        timeout_val=Timeout_text)
                    Result_list_RPT = timeout_RPT(
                        Model_Dry_old  , Sol_Dry_old ,   
                        Paraupdate,      Experiment_RPT, RPT_Cycles, 
                        Temper_RPT ,mesh_list ,submesh_strech
                    )
                else:
                    Result_list_RPT = Run_Model_Base_On_Last_Solution_RPT(
                        Model_Dry_old  , Sol_Dry_old ,   
                        Paraupdate,      Experiment_RPT, RPT_Cycles, 
                        Temper_RPT ,mesh_list ,submesh_strech
                    )
                [Model_Dry_i, Sol_Dry_i,Call_RPT,DeBug_List_RPT]  = Result_list_RPT
                if Return_Sol == True:
                    Sol_RPT.append(Sol_Dry_i)
                #print(f"Temperature for RPT is now: {Temper_RPT}")  
                if Call_RPT.success == False:
                    print("Fail due to Experiment error or infeasible")
                    str_error_RPT = "Experiment error or infeasible"
                    1/0 
                if Sol_Dry_i == Timeout_text:
                    #print("Fail due to Timeout")
                    str_error_RPT = "Timeout"
                    1/0
                if Sol_Dry_i == "Model error or solver error":
                    print("Fail due to Model error or solver error")
                    str_error_RPT = "Model error or solver error"
                    1/0
            except ZeroDivisionError as e:
                if Check_Small_Time == True:    
                    str_error_RPT = f"Scan {Scan_i} Re {Re_No}: Fail during No.{Cyc_Update_Index[-1]} RPT cycles within {SmallTimer.time()}, due to {str_error_RPT}"
                    print(str_error_RPT)
                    SmallTimer.reset()
                else:
                    str_error_RPT = f"Scan {Scan_i} Re {Re_No}: Fail during No.{Cyc_Update_Index[-1]} RPT cycles, due to {str_error_RPT}"
                    print(str_error_RPT)
                break
            else:
                # post-process for RPT
                Cyc_Update_Index.append(cycle_count)
                if Check_Small_Time == True:    
                    print(f"Scan {Scan_i} Re {Re_No}: Finish for No.{Cyc_Update_Index[-1]} RPT cycles within {SmallTimer.time()}")
                    SmallTimer.reset()
                else:
                    print(f"Scan {Scan_i} Re {Re_No}: Finish for No.{Cyc_Update_Index[-1]} RPT cycles")
                # update 231210: delete the first hold at 4.2V for later RPT
                my_dict_RPT = GetSol_dict (my_dict_RPT,keys_all_RPT, Sol_Dry_i, 
                    0,step_0p1C_CD, step_0p1C_CC,step_0p1C_RE , step_AGE_CV   ) 
                my_dict_RPT["Cycle_RPT"].append(cycle_count)
                my_dict_RPT["avg_Age_T"].append(np.mean(avg_Age_T))  # Make sure avg_Age_T and 
                
                # update 230517 - Get R from C/2 discharge only, discard GITT
                cap_full = Paraupdate["Nominal cell capacity [A.h]"] # 5

                # if R_from_GITT: 
                #     Res_midSOC,Res_full,SOC_Res = Get_0p1s_R0(Sol_Dry_i,Cyc_Index_Res,cap_full)
                # else: 
                #     step_0P5C_CD = Sol_Dry_i.cycles[0].steps[step_0p5C_CD]
                #     Res_midSOC,Res_full,SOC_Res = Get_R_from_0P5C_CD(step_0P5C_CD,cap_full)
                # my_dict_RPT["SOC_Res"].append(SOC_Res)
                # my_dict_RPT["Res_full"].append(Res_full)
                # my_dict_RPT["Res_midSOC"].append(Res_midSOC)             
                # del SOC_Res,Res_full,Res_midSOC

                if DryOut == "On":
                    mdic_dry = Update_mdic_dry(Data_Pack,mdic_dry)
                Para_0_Dry_old = Paraupdate;    Model_Dry_old = Model_Dry_i  ;     Sol_Dry_old = Sol_Dry_i    ;   
                del Paraupdate,Model_Dry_i,Sol_Dry_i
                if Check_Small_Time == True:    
                    print(f"Scan {Scan_i} Re {Re_No}: Finish post-process for No.{Cyc_Update_Index[-1]} RPT cycles within {SmallTimer.time()}")
                    SmallTimer.reset()
                else:
                    pass
                if Flag_AGE == False or Flag_partial_AGE == True:
                    break
            k += 1 
    DeBug_Lists = [DeBug_List_RPT,DeBug_List_AGE]
    Keys_error = ["Error tot %","Error SOH %","Error LLI %",
        "Error LAM NE %","Error LAM PE %",
        "Error Res %","Error ageT %","Punish"]
    ############################################################# 
    #########   An extremely bad case: cannot even finish breakin
    if Flag_Breakin == False: 
        log_messages = log_buffer.getvalue()
        log_messages = log_messages.strip()
        value_list_temp = list(Para_dict_i.values())
        values_para = []
        for value_list_temp_i in value_list_temp:
            values_para.append(str(value_list_temp_i))
        # sequence: scan no, exp, pass or fail, mpe, dry-out, 
        mpe_all = [np.nan] * 8
        Pass_Fail = "Fail at break-in cycle"
        # Update 240430 new script: 
        df_excel = Write_Dict_to_Excel(
            Flag_Breakin,Para_dict_i,cap_0,
            BasicPath,Target,book_name_xlsx,
            index_exp,Pass_Fail,mpe_all,
            model_options,my_dict_RPT,mdic_dry,
            DryOut,Scan_i,str_exp_AGE_text,
            str_exp_RPT_text,str_error_AGE_final,
            str_error_RPT,log_messages,)
        for key in Keys_error:
            my_dict_RPT[key] =  np.nan
        
        midc_merge = [my_dict_RPT, my_dict_AGE,mdic_dry]
        import pickle
        with open(
            BasicPath + Target+"Mats/" 
            + str(Scan_i)+ f'-DeBug_Lists_Re_{Re_No}.pkl', 'wb') as file:
            pickle.dump(DeBug_Lists, file)

        return midc_merge,Sol_RPT,Sol_AGE,DeBug_Lists
    ##########################################################
    ##############   Part-3: Post-prosessing    ##############
    ##########################################################
    # Newly add (220517): save plots, not just a single line in excel file:     
    # Newly add (221114): make plotting as functions
    # Niall_data = loadmat( 'Extracted_all_cell.mat'
    else:
        #if not os.path.exists(BasicPath + Target + str(Scan_i)):
        #os.mkdir(BasicPath + Target + str(Scan_i) );
        # Update 230221 - Add model LLI, LAM manually 
        my_dict_RPT = Get_SOH_LLI_LAM(my_dict_RPT,model_options,DryOut,mdic_dry,cap_0)

        if Check_Small_Time == True:    
            print(f"Scan {Scan_i} Re {Re_No}: Getting extra variables within {SmallTimer.time()}")
            SmallTimer.reset()
        else:
            pass
        # Newly add: update 23-05-18: evaluate errors systematically:
        if index_exp in list(np.arange(1,6)) and int(Temper_i- 273.15) in [10,25,40]:
            Exp_temp_i_cell = Temp_Cell_Exp[str(int(Temper_i- 273.15))]
            XY_pack = Get_Cell_Mean_1T_1Exp(Exp_Any_AllData,Exp_temp_i_cell) # interpolate for exp only
            mpe_all = Compare_Exp_Model( my_dict_RPT, XY_pack, Scan_i, Re_No,
                index_exp, Temper_i,BasicPath, Target,fs,dpi, PlotCheck=True)
        else:
            Exp_temp_i_cell = "nan"
            XY_pack         = "nan"
            mpe_all         = [np.nan]*8
        for mpe_i,key in zip(mpe_all,Keys_error):
            my_dict_RPT[key] =  mpe_i
        # [mpe_tot,mpe_1,mpe_2,mpe_3,mpe_4,mpe_5,mpe_6,punish] = mpe_all
        # set pass or fail TODO figure out how much should be appropriate:
        if isinstance(mpe_all[0],float): # is mpe_tot
            if mpe_all[0] < 3:
                Pass_Fail = "Pass"
            else:
                Pass_Fail = "Fail"
        else:
            Pass_Fail = f"Not Exp"

        ##########################################################
        #########      3-1: Plot cycle,location, Dryout related 
        # update 23-05-25 there is a bug in Cyc_Update_Index, need to slide a bit:
        Cyc_Update_Index.insert(0,0); del Cyc_Update_Index[-1]

        Plot_Cyc_RPT_4(
            my_dict_RPT, Exp_Any_AllData,Temp_Cell_Exp,
            XY_pack,index_exp, Plot_Exp,   R_from_GITT,
            Scan_i,Re_No,Temper_i,model_options,
            BasicPath, Target,fs,dpi)
        Plot_DMA_Dec(my_dict_RPT,Scan_i,Re_No,Temper_i,model_options,
            BasicPath, Target,fs,dpi)
        if len(my_dict_AGE["CDend Porosity"])>1:
            Plot_Loc_AGE_4(
                my_dict_AGE,Scan_i,Re_No,index_exp,Temper_i,
                model_options,BasicPath, Target,fs,dpi)
            colormap = "cool"
            Plot_HalfCell_V(
                my_dict_RPT,my_dict_AGE,Scan_i,Re_No,index_exp,colormap,
                Temper_i,model_options,BasicPath, Target,fs,dpi)
        if DryOut == "On":
            Plot_Dryout(
                Cyc_Update_Index,mdic_dry,ce_EC_0,index_exp,Temper_i,
                Scan_i,Re_No,BasicPath, Target,fs,dpi)
        if Check_Small_Time == True:    
            print(f"Scan {Scan_i} Re {Re_No}: Finish all plots within {SmallTimer.time()}")
            SmallTimer.reset()
        else:
            pass
        # update 23-05-26 judge ageing shape: contain 3 index
        """ [Start_shape,Middle_shape,End_shape] = check_concave_convex(
        my_dict_RPT['Throughput capacity [kA.h]'], 
        my_dict_RPT['CDend SOH [%]']) """
        ##########################################################
        ##########################################################
        #########      3-3: Save summary to excel 
        log_messages = log_buffer.getvalue()
        log_messages = log_messages.strip()
        df_excel = Write_Dict_to_Excel(
            Flag_Breakin,Para_dict_i,cap_0,
            BasicPath,Target,book_name_xlsx,
            index_exp,Pass_Fail,
            mpe_all,
            model_options,my_dict_RPT,mdic_dry,
            DryOut,Scan_i,str_exp_AGE_text,
            str_exp_RPT_text,str_error_AGE_final,
            str_error_RPT,log_messages,
            )
        # pick only the cyc - most concern
        Keys_cyc_mat =[
            'Throughput capacity [kA.h]',
            'CDend SOH [%]',
            "CDend LLI [%]",
            "CDend LAM_ne [%]",
            "CDend LAM_pe [%]",
            "Res_midSOC",
        ]
        my_dict_mat = {}
        for key in Keys_cyc_mat:
            my_dict_mat[key]=my_dict_RPT[key]
        #########      3-2: Save data as .mat or .json
        my_dict_RPT["Cyc_Update_Index"] = Cyc_Update_Index
        my_dict_RPT["SaveTimes"]    = SaveTimes
        # new_dict = {}
        # add "age" after age keys to avoid overwritten
        #for key in my_dict_AGE: 
        #    new_key = key + " Age"
        #    new_dict[new_key] = my_dict_AGE[key]
        # midc_merge = {**my_dict_RPT, **my_dict_AGE,**mdic_dry}
        midc_merge = [my_dict_RPT, my_dict_AGE,mdic_dry]
        if isinstance(Sol_RPT[-1],pb.solvers.solution.Solution):
            _,dict_short = Get_Last_state(Model_Dry_old, Sol_RPT[-1])
            sol_last = Sol_RPT[-1]
        else:
            _,dict_short = Get_Last_state(Model_Dry_old, Sol_AGE[-1])
            sol_last = Sol_AGE[-1]
            print("!!!!!!!Big problem! The last RPT fails, need to restart from RPT next time")
        # calculate dry-out parameter first
        if DryOut == "On":
            Data_Pack,Paraupdate   = Cal_new_con_Update (  sol_last,   Para_0_Dry_old )
        else: 
            Paraupdate = Para_0; Data_Pack = "nan"
        i_try = 0
        while i_try<3:
            try:
                getSth = sol_last['Throughput capacity [A.h]'].entries[-1] 
                # Sol_new['Throughput capacity [A.h]'].entries[-1]
            except:
                i_try += 1
                print(f"Fail to read Throughput capacity for the {i_try}th time")
            else:
                break
        Save_for_Reload = [ midc_merge, dict_short, Paraupdate, Data_Pack, getSth]
        import pickle,json

        with open(
            BasicPath + Target+"Mats/" 
            + str(Scan_i)+ f'_Re_{Re_No}-midc_merge.pkl', 'wb') as file:
            pickle.dump(midc_merge, file)
        
        with open(
            BasicPath + Target+"Mats/" 
            + str(Scan_i)+ f'_Re_{Re_No}-Save_for_Reload.pkl', 'wb') as file:
            pickle.dump(Save_for_Reload, file)


        try:
            savemat(
                BasicPath + Target+"Mats/" 
                + str(Scan_i)+ f'_Re_{Re_No}-Ageing_summary_only.mat',
                my_dict_mat)  
        except:
            print(f"Scan {Scan_i} Re {Re_No}: Encounter problems when saving mat file!")
        else: 
            print(f"Scan {Scan_i} Re {Re_No}: Successfully save mat file!")

        if Check_Small_Time == True:    
            print(f"Scan {Scan_i} Re {Re_No}: Try saving within {SmallTimer.time()}")
            SmallTimer.reset()
        else:
            pass

        
        with open(
            BasicPath + Target+"Mats/" 
            + str(Scan_i)+ f'_Re_{Re_No}-DeBug_Lists.pkl', 'wb') as file:
            pickle.dump(DeBug_Lists, file)
        

        # update 231217: save ageing solution if partially succeed in ageing set
        if Flag_partial_AGE == True:
            try:
                Sol_partial_AGE_list = [
                    Sol_AGE[-1].cycles[0],    Sol_AGE[-1].cycles[-1],
                    len(Sol_AGE[-1].cycles)]
            except IndexError:
                try:
                    Sol_partial_AGE_list = [
                        Sol_AGE[-1].cycles[0],    Sol_AGE[-1].cycles[-2],
                        len(Sol_AGE[-1].cycles)]
                except IndexError:
                    print("Problems in saving last two cycles in Sol_AGE[-1],"
                            " now saving the whole Sol_AGE[-1]")
                    Sol_partial_AGE_list = [
                        Sol_AGE[-1],    "nan",
                        len(Sol_AGE[-1].cycles)]
            with open(
                BasicPath + Target+"Mats/" 
                + str(Scan_i)+ f'_Re_{Re_No}-Sol_partial_AGE_list.pkl', 'wb') as file:
                pickle.dump(Sol_partial_AGE_list, file)
            print(f"Last AGE succeed partially, save Sol_partial_AGE_list.pkl for Scan {Scan_i} Re {Re_No}")
        else:
            pass
        print("Succeed doing something in {}".format(ModelTimer.time()))
        print(f'This is the end of No. {Scan_i} scan, Re {Re_No}')
        return midc_merge,Sol_RPT,Sol_AGE,DeBug_Lists


