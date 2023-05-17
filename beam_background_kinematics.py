import matplotlib
import matplotlib.pyplot as plt
import h5py
import argparse
import numpy as np
import twoBytwo_defs
import threshold_backgrounds
import auxiliary
import glob

nu_signal_pdg=-14
pion_pdg={111,211,-211}




def signal_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if ghdr_nu_interaction[0]==nu_signal_pdg: return True
    else: return False
    


def signal_cc(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    return ghdr[ghdr_vert_mask]['isCC'][0]

    

def signal_pion_status(gstack, vert_id):
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    if len(pion_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False


    
def main(sim_file, input_type):
    cc_dict, nc_dict, cc_primaries_dict, nc_primaries_dict = [dict() for i in range(4)]
    file_ctr=0

    file_ext = '' ### modified by commandline argument
    if input_type=='larnd': file_ext='.LARNDSIM.h5'
    elif input_typw=='edep': file_ext='.EDEPSIM.h5'
    
    for sim_file in glob.glob(sim_dir+'/*'+file_ext):
        file_ctr+=1
        if file_ctr>30: break
        if file_ctr%10==0: print('FILE #: ',file_ctr)
        sim_h5 = h5py.File(sim_file,'r')
    
        ### partition file by spill
        unique_spill = np.unique(sim_h5['trajectories']['spillID'])
        for spill_id in unique_spill:
            ghdr, gstack, traj, vert, seg = auxiliary.get_spill_data(sim_h5, spill_id)

            ### partition by vertex ID within beam spill
            for v_i in range(len(vert['vertexID'])):
                vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
                vert_in_active_LAr = twoBytwo_defs.fiducialized_vertex( vert_pos )

                ##### REQUIRE neutrino vertex in LAr active volume #####
                if vert_in_active_LAr==False: continue

                vert_id = vert['vertexID'][v_i]
            
                nu_mu_bar = signal_nu_pdg(ghdr, vert_id)
                is_cc = signal_cc(ghdr, vert_id)
                pionless = signal_pion_status(gstack, vert_id)
                fv_particle_origin=twoBytwo_defs.fiducialized_particle_origin(traj, vert_id)
                        
                ##### THRESHOLD BACKGROUNDS #####
                if is_cc==True and pionless==False and fv_particle_origin==True:
                    threshold_backgrounds.pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, cc_dict)
                    threshold_backgrounds.primaries(spill_id, vert_id, ghdr, gstack, traj, vert, seg, cc_primaries_dict)

                ##### PID BACKGROUNDS #####
                if is_cc==False and pionless==False and fv_particle_origin==True:
                    threshold_backgrounds.pion_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, nc_dict)
                    threshold_backgrounds.primaries(spill_id, vert_id, ghdr, gstack, traj, vert, seg, nc_primaries_dict)

    auxiliary.save_dict_to_json(cc_dict, 'cc_pion_backgrounds', True)
    auxiliary.save_dict_to_json(cc_primaries_dict, 'cc_primaries', True)
    auxiliary.save_dict_to_json(nc_dict, 'nc_pion_backgrounds', True)
    auxiliary.save_dict_to_json(nc_primaries_dict, 'nc_primaries', True)
    


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--sim_dir', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the directory containing edep-sim or larnd ouput simulation file(s)''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, \
                        help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
