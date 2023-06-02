import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import auxiliary
import signal_characterization as sig_char
from plot_signal_muons import plot_muons
from plot_signal_hadrons import plot_hadrons

def main(sim_dir, input_type, n_files_processed):

    test_count = 0

    ### NOTE: Current POT scaling is based on MiniRun3 larnd file situation
    if int(n_files_processed) < 1023.: 
        scale_factor = (1./(int(n_files_processed)/1023.))*2.5
    else:
        scale_factor = 2.5

    muon_dict = dict() # Initialize muon dictionary
    hadron_dict = dict() # Initialize hadron dictionary

    # Dictionaries for combining with other background explorations
    signal_dict = dict() # Initialize dictionary for signal muons for full comparison
    
    file_ext = '' ## Changes based on input type

    if input_type == 'larnd': 
        file_ext = '.LARNDSIM.h5'
    elif input_type == 'edep':
        file_ext = '.EDEPSIM.h5'

    for sim_file in glob.glob(sim_dir+'/*'+file_ext): # Loop over simulation files

        if test_count ==int(n_files_processed) : break
        test_count+=1

        if (test_count % 5 == 0):
            print("Processing file: ", str(test_count), "/", str(n_files_processed))

        sim_h5 = h5py.File(sim_file,'r')

        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories']['eventID'])
        for spill_id in unique_spill:

            ghdr, gstack, traj, vert, seg = auxiliary.get_spill_data(sim_h5, spill_id, input_type)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):

                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]] 
                vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos ) # Check vertex location relative to FV

                ##### REQUIRE: neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: continue

                vert_id = vert['vertexID'][v_i]

                nu_mu = auxiliary.signal_nu_pdg(ghdr, vert_id) # nu_mu OR nu_mu_bar
                is_cc = auxiliary.signal_cc(ghdr, vert_id)
                mesonless = auxiliary.signal_meson_status(gstack, vert_id)
                fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)

                ### REQUIRE: (A) nu_mu(_bar), (B) CC interaction, (C) NO final state mesons, (D) final state particle start point in FV
                if nu_mu==True and is_cc==True and mesonless==True and fv_particle_origin==True:
                    sig_char.muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict)
                    sig_char.hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict)
                    sig_char.get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict)

    # Save all Python dictionaries to JSON files
    auxiliary.save_dict_to_json(signal_dict, "signal_dict", True)
    auxiliary.save_dict_to_json(muon_dict, "muon_dict", True)
    auxiliary.save_dict_to_json(hadron_dict, "hadron_dict", True)

    # Save full signal and w.s. bkg counts to TXT file
    signal_count = len(signal_dict)*scale_factor
    outfile = open('signal_event_counts.txt', "w")
    outfile.writelines(["Signal Events (scaled to 2.5e19 POT): "+str(signal_count)+"\n", \
                        "Number of files used to get count: "+str(n_files_processed)+"\n", \
                        "Scale factor:"+str(scale_factor)+"\n"])
    outfile.close()

    # PLOT: Signal Event Info      
    plot_muons(muon_dict, scale_factor, sig_bkg = 0)
    plot_hadrons(hadron_dict, scale_factor, sig_bkg = 0)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='larnd', choices=['edep', 'larnd'], type=str, \
                        help='''string corresponding to the output file type: edep or larnd''')
    parser.add_argument('-n', '--n_files_processed', default=1, required=True, type=int, \
                        help='''File count of number of files processed in production sample''')
    args = parser.parse_args()
    main(**vars(args))
