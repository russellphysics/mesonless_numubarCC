################################################################################
##                                                                            ##
##    CONTAINS: Methods to parse files or objects and/or assist in new file/  ##
##              object creation e.g. hdf5 file parsing and json dictionary    ##
##              creation.                                                     ##
##                                                                            ##
################################################################################
import json

####--------------------------- HDF5 FILE PARSING --------------------------####

def print_keys_attributes(sim_h5, input_type):
    print(sim_h5.keys(),'\n')
    if input_type=='edep':
        print('GENIE HDR: ',sim_h5['genie_hdr'].dtype,'\n')
        print('GENIE STACK: ',sim_h5['genie_stack'].dtype,'\n')
        print('SEGMENTS: ', sim_h5['segments'].dtype,'\n')
        print('TRAJECTORIES', sim_h5['trajectories'].dtype,'\n')
        print('VERTICES', sim_h5['vertices'].dtype)
    elif input_type=='larnd':
        print('GENIE HDR: ',sim_h5['genie_hdr'].dtype,'\n')
        print('GENIE STACK: ',sim_h5['genie_stack'].dtype,'\n')
        print('LIGHT DAT: ',sim_h5['light_dat'].dtype,'\n')
        print('LIGHT TRIGGER: ',sim_h5['light_trig'].dtype,'\n')
        print('LIGHT WAVEFORM: ',sim_h5['light_wvfm'].dtype,'\n')
        print('MC PACKETS ASSN: ',sim_h5['mc_packets_assn'].dtype,'\n')
        print('MESSAGES: ',sim_h5['messages'].dtype,'\n')
        print('PACKETS: ',sim_h5['packets'].dtype,'\n')
        print('TRACKS: ',sim_h5['tracks'].dtype,'\n')
        print('TRAJECTORIES', sim_h5['trajectories'].dtype,'\n')
        print('VERTICES', sim_h5['vertices'].dtype)
        

def get_spill_data(sim_h5, spill_id, input_type):
    parse_var=''; seg_var=''
    if input_type=='edep': parse_var='spillID'; seg_var='segments'
    elif input_type=='larnd': parse_var='eventID'; seg_var='tracks'        

    ghdr_spill_mask = sim_h5['genie_hdr'][:][parse_var]==spill_id
    gstack_spill_mask = sim_h5['genie_stack'][:][parse_var]==spill_id
    traj_spill_mask = sim_h5['trajectories'][:][parse_var]==spill_id
    vert_spill_mask = sim_h5['vertices'][:][parse_var]==spill_id
    seg_spill_mask = sim_h5[seg_var][:][parse_var]==spill_id

    ghdr = sim_h5['genie_hdr'][ghdr_spill_mask]
    gstack = sim_h5['genie_stack'][gstack_spill_mask]
    traj = sim_h5['trajectories'][traj_spill_mask]
    vert = sim_h5['vertices'][vert_spill_mask]
    seg = sim_h5[seg_var][seg_spill_mask]

    return ghdr, gstack, traj, vert, seg


####----------------------- OUTPUT DICTIONARY TO JSON ----------------------####

def tuple_key_to_string(d):
    out={}
    for key in d.keys():
        string_key=""
        max_length=len(key)
        for i in range(max_length):
            if i<len(key)-1: string_key+=str(key[i])+"-"
            else: string_key+=str(key[i])
        out[string_key]=d[key]
    return out                    
            

def save_dict_to_json(d, name, if_tuple):
    with open(name+".json", "w") as outfile:
        if if_tuple==True:
            updated_d = tuple_key_to_string(d)
            json.dump(updated_d, outfile, indent=4)
        else:
            json.dump(d, outfile, indent=4)    


####-------------------- NUMPY/PYTHON OBJECT CONVERSIONS -------------------####

def np_array_of_array_to_flat_list(a):
    b = list(a)
    return [list(c)[0] for c in b]
