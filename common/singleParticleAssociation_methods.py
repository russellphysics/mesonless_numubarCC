################################################################################
##                                                                            ##
##    CONTAINS: Methods related to associating segments and/or trajectories   ##
##              corresponding to the same particle. ALL METHODS CURRENTLY     ##
##              RELY ON TRUTH INFORMATION (07/11/2023).                       ##
##    NOTES: "ghdr" refers to a genie_hdr dataset extracted from an hdf5      ##
##           file as seen in common/file_parsing.py                           ##
##           "traj" refers to trajectories dataset extracted from an hdf5     ##
##           file as seen in common/file_parsing.py                           ##
##                                                                            ##
################################################################################

import truth_methods


####------------------- CONNECT TRAJECTORIES USING PDG ---------------------####

# Method below used to connect trajectories of a single particle which have been  
# assigned different traj_ids. If a particle scatters or reinteracts in a way such 
# that multiple child particles are produced, the lineage stops, even if one of 
# those child particles has the same PDG ID as the original particle. We do this by 
# ensuring that the original particle/traj_id (and any other traj_id accepted as 
# being the same particle) does not have the same parent as any other 
# particles/traj_ids (is an "only child") and also does not have multiple children.
def same_pdg_connected_trajectories(track_pdg, track_id, vertex_assoc_traj,\
                                    traj, ghdr):
    traj_id_set = {track_id} # initialize a set of traj_ids associated with a single particle
    
    ## WALK UP THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg: # stop if a member of an older generation has a different PDG ID than the original track/particle
        particle_mask = vertex_assoc_traj['traj_id'] == this_track_id # mask to find trajectories assoc. w/ current track
        parent_track_id = vertex_assoc_traj[particle_mask]['parent_id'] # find parent ID of current track
        parent_mask = vertex_assoc_traj['parent_id'] == parent_track_id # mask to find trajectories w/ current track AS PARENT TRACK
        this_generation =  vertex_assoc_traj[parent_mask] # find all tracks with the same parent as the current track
        if len(this_generation) == 1: # only move forward if current track is an "only child"
            this_pdg = truth_methods.find_parent_pdg(vertex_assoc_traj[particle_mask]['parent_id'],
                                   vertex_assoc_traj[particle_mask]['vertex_id'],
                                   traj, ghdr) # get parent PDG ID of current track
            if this_pdg==track_pdg: # if parent PDG ID of track matches current/original track's PDG ID, add parent to track id set
                this_track_id = vertex_assoc_traj[particle_mask]['parent_id'].tolist()[0] # also makes parent track the new "current" track
                traj_id_set.add(this_track_id)
        else: break # break while loop if current track/particle is not an "only child"

    ## WALK DOWN THE FAMILY TREE
    this_pdg = track_pdg
    this_track_id = track_id

    while this_pdg==track_pdg: # stop if/when a child has a different PDG ID than the original track
        particle_mask = vertex_assoc_traj['parent_id'] == this_track_id  # mask to find trajectories w/ current track AS PARENT TRACK
        child_particle=vertex_assoc_traj[particle_mask] # find all tracks which are children of the current track
        if len(child_particle)==1: # only move forward if current track has only one child
            this_pdg = child_particle['pdg_id']
            if child_particle['pdg_id']==track_pdg: # if the child's PDG ID matches the original track PDG ID, add child track id to set
                this_track_id = child_particle['traj_id'].tolist()[0] # also makes child track the new "current" track
                traj_id_set.add(this_track_id)
        else: break # break while loop if current track/particle has more than one child particle

    return traj_id_set
    

####-------- FIND TRAJECTORY AT END OR BEGINNING OF PARTICLE TRACK ---------####

def find_trajectory_at_vertex(traj_id_set, vertex_assoc_traj,traj, ghdr):
    
    is_prim = truth_methods.is_primary_particle(traj_id_set, vertex_assoc_traj,traj, ghdr)

    if is_prim == False:
        print("Track ID set does not represent a primary particle. Therefore, \
              there does not exist a trajectory coming from the vertex.")
        return 
    else:
        traj_id_at_vertex = 0 # initialize traj_id_at_vertex variable
        for tid in traj_id_set: # loop through trajectories in set to find trajectory with muon (anti)neutrino parent
            particle_mask = vertex_assoc_traj['traj_id'] == tid
            parent_pdg = truth_methods.find_parent_pdg(vertex_assoc_traj[particle_mask]['parent_id'],
                                         vertex_assoc_traj[particle_mask]['vertex_id'],
                                         traj, ghdr) 
            
            if abs(parent_pdg)==14:
                traj_id_at_vertex = tid
                break
            else: continue

    return traj_id_at_vertex


def find_forward_primary_particle_end_trajectory(traj_id_set, vertex_assoc_traj,traj, ghdr):
    
    tid_at_vertex = find_trajectory_at_vertex(traj_id_set, vertex_assoc_traj,traj, ghdr)
    particle_mask = vertex_assoc_traj['traj_id'] == tid_at_vertex
    start_pt = vertex_assoc_traj[particle_mask]['xyz_start']

    traj_id_at_end = tid_at_vertex # initialize traj_id_at_vertex variable
    end_z = start_pt[2] # initialize end z value
    for tid in traj_id_set: # loop through trajectories in set to find trajectory with largest z value
        traj_mask = vertex_assoc_traj['traj_id'] == tid
        tid_end = vertex_assoc_traj[traj_mask]['xyz_end']
        if tid_end[2]>end_z:
            end_z = tid_end[2]
            traj_id_at_end = tid
        else: continue

    return traj_id_at_end