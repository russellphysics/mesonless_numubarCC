################################################################################
##                                                                            ##
##    CONTAINS: Script to create a JSON dictionary and plots describing muons ##
##              in dirt background events. ONLY TRUTH INFORMATION IS USED.    ##
##                                                                            ##
################################################################################

import h5py
import argparse
import numpy as np
import glob
import os
import dirt_backgrounds
import sys
sys.path.append('../../common')
import file_parsing
import geometry_methods as geo_methods
import particlePDG_defs as pdg_defs
import truth_methods as truth
sys.path.append('../plotting_scripts')
from plot_dirt_background import plot_dirt_backgrounds

    
def main(sim_dir, input_type):
 
    #sim_h5 = h5py.File(sim_file,'r')
    dirt_muon_dict = dict()
    file_count =0
    
    #print('start')
    files = glob.glob(sim_dir+'*h5')
    files.sort(key=os.path.getctime)
    #print(files)
    
    for sim_file in files:
        file_count+=1
        if file_count ==10:
            print('search aborted, File count: ', file_count)
            break
        
        sim_h5 = h5py.File(sim_file, 'r')
        print('Openning new file: ', sim_file)
        #file_parsing.print_keys_attributes(sim_h5)    
        print('file count: ', file_count)
        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories']['eventID']) #spillID
        for spill_id in unique_spill:

            ghdr, gstack, traj, vert, seg = file_parsing.get_spill_data(sim_h5, spill_id, input_type)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):
                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
                vert_in_active_LAr = geo_methods.fiducialized_vertex( vert_pos )
                #print(vert_pos)
                
                ##### REQUIRE neutrino vertex out of LAr active volume #####
                if vert_in_active_LAr==True: continue

                vert_id = vert['vertexID'][v_i]

                #nu_mu_bar = pdg_defs.signal_nu_pdg(ghdr, vert_id)
                #is_cc = truth.signal_cc(ghdr, vert_id)
                #pionless = truth.signal_meson_status(gstack, vert_id)
                
                fv_particle_origin=geo_methods.fiducialized_particle_origin(traj, vert_id)

                ##### Dirt Muon BACKGROUNDS #####
                ##### REQUIRE: (A) nu_mu_bar, (B) CC, (C) pions present, (D) final state particle start point in FV
                #if nu_mu_bar==True and is_cc==True and pionless==False and fv_particle_origin==True:

                if fv_particle_origin==True: ## change the criteria
                    dirt_backgrounds.dirt_muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, dirt_muon_dict)

            # end vertex loop
        #end spill loop
    #end of file loop

    plot_dirt_backgrounds(dirt_muon_dict,1) #plotting

    file_parsing.save_dict_to_json(dirt_muon_dict, "dirt_muon_dict", True)

    print (dirt_muon_dict)
    
    #end
    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, type = str, help = ''' string path to the directory of simulation files ''' )
    #parser.add_argument('-f', '--sim_file', default=None, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
