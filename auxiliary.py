import twoBytwo_defs
import numpy as np
import json

neutral_pdg=[111] #, 22] #, 2112] # add K0, rho0, eta0?
meson_pdg={111,211,-211,130,310,311,321,-321,221,331,421,-421,411,-411, 431,-431}
nu_mu_pdg=14

hadron_pdg_dict ={2112:'n',
                  2212:'p',
                 -2212:r'$\bar{p}$', 
                  3112:r'$\Sigma^-$',
                  3122:r'$\Lambda^0$',
                 -3122:r'$\bar{\Lambda}^0$',
                  3212:r'$\Sigma^0$',
                  3222:r'$\Sigma^+$', 
                  4212:r'$\Sigma_c^+$',
                  4222:r'$\Sigma_c^{++}$',
                  4112:r'$\Sigma_c^0$', 
                  4122:r'$\Lambda_c^+$'} 

neutral_hadron_pdg_dict ={2112:'n',
                          3122:r'$\Lambda^0$',
                         -3122:r'$\bar{\Lambda}^0$',
                          3212:r'$\Sigma^0$',
                          4112:r'$\Sigma_c^0$'} 

threshold = 1.2 # stand-in threshold in cm (~3 pixels)
rest_mass_dict ={2212: 938.27208816,
                   13: 105.6583755  }  # [MeV], from PDG

##### HDF5 FILE PARSING-------------------------------------


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
    


##### OUTPUT DICTIONARY TO JSON------------------------------


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


##### SIGNAL BASIC CHARACTERISTICS--------------------------
def signal_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if abs(ghdr_nu_interaction[0])==nu_mu_pdg: return True
    else: return False


def signal_cc(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    return ghdr[ghdr_vert_mask]['isCC'][0]

def signal_meson_status(gstack, vert_id):
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    if len(meson_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False

def nu_int_type(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertexID']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]
    int_type = ''
    if ghdr_nu_interaction['isQES'] == True:
        int_type = 'QES'
    elif ghdr_nu_interaction['isMEC'] == True:
        int_type = 'MEC'
    elif ghdr_nu_interaction['isRES'] == True:
        int_type = 'RES'
    elif ghdr_nu_interaction['isDIS'] == True:
        int_type = 'DIS'
    elif ghdr_nu_interaction['isCOH'] == True:
        int_type = 'COH'
    else:
        int_type = 'UND'
    return int_type

    
##### FIND PARENT PDG --------------------------------------
def find_parent_pdg(parent_id, vertex_id, traj, ghdr):
    if parent_id==-1:
        ghdr_mask = ghdr['vertexID']==vertex_id
        parent_pdg=ghdr[ghdr_mask]['nu_pdg']
    else:
        parent_mask = traj['trackID']==parent_id
        parent_pdg = traj[parent_mask]['pdgId']
    if parent_pdg==[]: parent_pdg=[0]
    return parent_pdg

            

##### SAME PDG CONNECTED TRAJECTORIES -----------------------

# Method to connect trajectories of a single particle which have been assigned 
# different trackIDs. If a particle scatters or reinteracts in a way such that 
# multiple child particles are produced, the lineage stops, even if one of those 
# child particles has the same PDG ID as the original particle. We do this by 
# ensuring that the original particle/trackID (and any other trackID accepted as 
# being the same particle) does not have the same parent as any other 
# particles/trackIDs (is an "only child") and also does not have multiple children.
def same_pdg_connected_trajectories(track_pdg, track_id, vertex_assoc_traj,\
                                    traj, ghdr):
    trackid_set = {track_id} # initialize a set of trackIDs associated with a single particle
    
    ## WALK UP THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg: # stop if a member of an older generation has a different PDG ID than the original track/particle
        particle_mask = vertex_assoc_traj['trackID'] == this_track_id # mask to find trajectories assoc. w/ current track
        parent_track_id = vertex_assoc_traj[particle_mask]['parentID'] # find parent ID of current track
        parent_mask = vertex_assoc_traj['parentID'] == parent_track_id # mask to find trajectories w/ current track AS PARENT TRACK
        this_generation =  vertex_assoc_traj[parent_mask] # find all tracks with the same parent as the current track
        if len(this_generation) == 1: # only move forward if current track is an "only child"
            this_pdg = find_parent_pdg(vertex_assoc_traj[particle_mask]['parentID'],
                                   vertex_assoc_traj[particle_mask]['vertexID'],
                                   traj, ghdr) # get parent PDG ID of current track
            if this_pdg==track_pdg: # if parent PDG ID of track matches current/original track's PDG ID, add parent to track id set
                this_track_id = vertex_assoc_traj[particle_mask]['parentID'].tolist()[0] # also makes parent track the new "current" track
                trackid_set.add(this_track_id)
        else: break # break while loop if current track/particle is not an "only child"

    ## WALK DOWN THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg: # stop if/when a child has a different PDG ID than the original track
        particle_mask = vertex_assoc_traj['parentID'] == this_track_id  # mask to find trajectories w/ current track AS PARENT TRACK
        child_particle=vertex_assoc_traj[particle_mask] # find all tracks which are children of the current track
        if len(child_particle)==1: # only move forward if current track has only one child
            this_pdg = child_particle['pdgId']
            if child_particle['pdgId']==track_pdg: # if the child's PDG ID matches the original track PDG ID, add child track id to set
                this_track_id = child_particle['trackID'].tolist()[0] # also makes child track the new "current" track
                trackid_set.add(this_track_id)
        else: break # break while loop if current track/particle has more than one child particle

    return trackid_set
    

def is_primary_particle(trackid_set, vertex_assoc_traj,traj, ghdr):
    is_prim = False

    for tid in trackid_set:

        particle_mask = vertex_assoc_traj['trackID'] == tid
        parent_pdg = find_parent_pdg(vertex_assoc_traj[particle_mask]['parentID'],
                                     vertex_assoc_traj[particle_mask]['vertexID'],
                                     traj, ghdr)
        #print("Parent PDG:", parent_pdg)
        if abs(parent_pdg)==14:
            is_prim = True
            break
        else: continue

    return is_prim


def find_trajectory_at_vertex(trackid_set, vertex_assoc_traj,traj, ghdr):
    
    is_prim = is_primary_particle(trackid_set, vertex_assoc_traj,traj, ghdr)

    if is_prim == False:
        print("Track ID set does not represent a primary particle. Therefore, there does not exist a trajectory coming from the vertex.")
        return 
    else:
        trackid_at_vertex = 0 # initialize trackid_at_vertex variable
        for tid in trackid_set: # loop through trajectories in set to find trajectory with muon (anti)neutrino parent
            particle_mask = vertex_assoc_traj['trackID'] == tid
            parent_pdg = find_parent_pdg(vertex_assoc_traj[particle_mask]['parentID'],
                                         vertex_assoc_traj[particle_mask]['vertexID'],
                                         traj, ghdr) 
            
            if abs(parent_pdg)==14:
                trackid_at_vertex = tid
                break
            else: continue

    return trackid_at_vertex


def angle_wrt_beam_direction(trackid_set, vertex_assoc_traj,traj, ghdr):

    trackid_at_vertex = find_trajectory_at_vertex(trackid_set, vertex_assoc_traj,traj, ghdr) # find track id for trajectory at vertex

    trackid_at_vertex_mask = vertex_assoc_traj['trackID']==trackid_at_vertex
    track_at_vertex = vertex_assoc_traj[trackid_at_vertex_mask] # primary particle trajectory from vertex

    particle_start = track_at_vertex['xyz_start'][0]
    particle_end = track_at_vertex['xyz_end'][0]
    #print("Particle Start: ", particle_start)

    # If the particle goes along the beam axis, we expect x_start == x_end and y_start == y_end.
    # This means that the length would be z_start - z_end. Comparing this value to the true trajectory length
    # allows us to find the angle between the particle trajectory and the beam (z) axis.
    particle_vector = particle_end - particle_start # get vector describing particle path
    #print("Particle Vector: ", particle_vector)
    particle_vector_length = np.sqrt(np.sum(particle_vector**2))
    beam_direction_length = particle_vector[2]
    angle_wrt_beam = np.arccos(beam_direction_length / particle_vector_length) # angle is returned in radians

    #print("Angle w/ resp. to beam:", angle_wrt_beam)
    return angle_wrt_beam 


def angle_between_two_trajectories(traj_trackid_one, traj_trackid_two, vertex_assoc_traj):

    traj_trackid_one_mask = vertex_assoc_traj['trackID']==traj_trackid_one
    traj_one = vertex_assoc_traj[traj_trackid_one_mask] # trajectory #1
    #print("Trajectory track ID One:", traj_trackid_one)

    traj_trackid_two_mask = vertex_assoc_traj['trackID']==traj_trackid_two
    traj_two = vertex_assoc_traj[traj_trackid_two_mask] # trajectory #2
    #print("Trajectory track ID Two:", traj_trackid_two)

    particle_one_start = traj_one['xyz_start'][0]
    particle_one_end = traj_one['xyz_end'][0]
    #print("Particle One Start:", particle_one_start)

    particle_two_start = traj_two['xyz_start'][0]
    particle_two_end = traj_two['xyz_end'][0]
    #print("Particle Two Start:", particle_two_start)

    if not (particle_one_start[0] == particle_two_start[0] and \
        particle_one_start[1] == particle_two_start[1] and \
        particle_one_start[2] == particle_two_start[2]):
        print("WARNING: You seem to be comparing two trajectories originating from different locations ...")

    particle_one_vector = particle_one_end - particle_one_start # get vector describing particle one path
    particle_two_vector = particle_two_end - particle_two_start # get vector describing particle one path
    
    particle_one_vector_length = np.sqrt(np.sum(particle_one_vector**2))
    particle_two_vector_length = np.sqrt(np.sum(particle_two_vector**2))
    #print("Particle One Vector Length:", particle_one_vector_length)
    #print("Particle Two Vector Length:", particle_two_vector_length)

    particle_vector_dot_product = np.sum(abs(particle_one_vector * particle_two_vector))

    angle = np.arccos(particle_vector_dot_product / \
                      (particle_one_vector_length * particle_two_vector_length)) # angle is returned in radians using a \dot b = |a||b|cos\theta

    #print("Angle b/w Two Trajectories:", angle)
    return angle

            
##### FIDUCIAL VOLUME/ TOTAL VOLUME ENERGY DEPOSITION -------


def tpc_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_e={}
    for i in range(8): contained_e[i]=0.
    for sg in seg[seg_id_mask]:
        tpc_fv=twoBytwo_defs.tpc_vertex([(sg['x_start']+sg['x_end'])/2.,
                                         (sg['y_start']+sg['y_end'])/2.,
                                         (sg['z_start']+sg['z_end'])/2.])
        for key in tpc_fv.keys():
            if tpc_fv[key]==True:
                contained_e[key]+=sg['dE']
    return contained_e


            
def fv_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_e=0.
    for sg in seg[seg_id_mask]:
        if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                              (sg['y_start']+sg['y_end'])/2.,
                                              (sg['z_start']+sg['z_end'])/2.]):
            contained_e+=sg['dE']
    return contained_e



def total_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    total_e=0.
    for sg in seg[seg_id_mask]: total_e+=sg['dE']
    return total_e



def fv_edep_neutral_e(trackID, traj, seg):
    flag=True; daughter_track_ids=set()
    temp_track_id=[trackID]
    while flag==True:
        second_temp_track_id=[]
        for ttid in temp_track_id:
            parent_mask = traj['parentID']==trackID
            daughter_id=traj[parent_mask]['trackID']
            if len(daughter_id)!=0:
                for did in daughter_id:
                    daughter_track_ids.add(did)
                    second_temp_track_id.append(did)
            else: flag=False
        temp_track_id=second_temp_track_id
    contained_e=0
    print('[CONTAINED EDEP] NEUTRAL daughter track IDs: ', daughter_track_ids)
    for dti in daughter_track_ids:
        contained_e+=fv_edep_charged_e(dti, traj, seg)
    return contained_e



def total_edep_neutral_e(trackID, traj, seg):
    flag=True; daughter_track_ids=set()
    temp_track_id=[trackID]
    while flag==True:
        second_temp_track_id=[]
        for ttid in temp_track_id:
            parent_mask = traj['parentID']==trackID
            daughter_id=traj[parent_mask]['trackID']
            if len(daughter_id)!=0:
                for did in daughter_id:
                    daughter_track_ids.add(did)
                    second_temp_track_id.append(did)
            else: flag=False
        temp_track_id=second_temp_track_id
    total_e=0
    print('[TOTAL EDEP] NEUTRAL daughter track IDs: ', daughter_track_ids)
    for dti in daughter_track_ids:
        total_e+=total_edep_charged_e(dti, traj, seg)
    return total_e


def tpc_contained_energy(pdg, trackID, traj, seg):
    if pdg in neutral_pdg:
        temp={}
        for i in range(8): temp[i]=0.
        return temp
    else: return tpc_edep_charged_e(trackID, traj, seg)

            

def fv_contained_energy(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.#fv_edep_neutral_e(trackID, traj, seg)
    else: return fv_edep_charged_e(trackID, traj, seg)


    
def total_energy(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.#total_edep_neutral_e(trackID, traj, seg)
    else: return total_edep_charged_e(trackID, traj, seg)

def truth_primary_particle_kinetic_energy(pdg, track_id_set, vertex_assoc_traj, traj, ghdr):

    traj_id_at_vertex = find_trajectory_at_vertex(track_id_set, vertex_assoc_traj, traj, ghdr)
    trackid_at_vertex_mask = vertex_assoc_traj['trackID']==traj_id_at_vertex
    track_at_vertex = vertex_assoc_traj[trackid_at_vertex_mask] 
    energy = track_at_vertex['E_start']
    ke = energy - rest_mass_dict[pdg]

    return ke
    
def truth_primary_particle_momentum(track_id_set, vertex_assoc_traj, traj, ghdr):

    traj_id_at_vertex = find_trajectory_at_vertex(track_id_set, vertex_assoc_traj, traj, ghdr)
    trackid_at_vertex_mask = vertex_assoc_traj['trackID']==traj_id_at_vertex
    track_at_vertex = vertex_assoc_traj[trackid_at_vertex_mask] 
    mom = np.sqrt(np.sum(track_at_vertex['pxyz_start']**2))

    return mom


##### FIDUCIAL VOLUME/ TOTAL VOLUME LENGTH ------------------------
    


def fv_edep_charged_length(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_length=0.
    for sg in seg[seg_id_mask]:
        if twoBytwo_defs.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                              (sg['y_start']+sg['y_end'])/2.,
                                              (sg['z_start']+sg['z_end'])/2.]):
            contained_length+=np.sqrt( (sg['x_start']-sg['x_end'])**2+
                                       (sg['y_start']-sg['y_end'])**2.+
                                       (sg['z_start']-sg['z_end'])**2. )
    return contained_length



def total_edep_charged_length(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_length=0.
    for sg in seg[seg_id_mask]:
            contained_length+=np.sqrt( (sg['x_start']-sg['x_end'])**2+
                                       (sg['y_start']-sg['y_end'])**2.+
                                       (sg['z_start']-sg['z_end'])**2. )
    return contained_length



def fv_contained_length(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.
    else: return fv_edep_charged_length(trackID, traj, seg)

    

def total_length(pdg, trackID, traj, seg):
    if pdg in neutral_pdg: return 0.
    else: return total_edep_charged_length(trackID, traj, seg)



##### MISCELLANEOUS ----------------------------------------------------------------------------

def np_array_of_array_to_flat_list(a):
    b = list(a)
    return [list(c)[0] for c in b]
