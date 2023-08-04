################################################################################
##                                                                            ##
##    CONTAINS: Methods and definitions related to determining particle       ##
##              kinematic variables (e.g. angle wrt beam, momentum, energy).  ##
##              ALL METHODS CURRENTLY RELY ON TRUTH INFORMATION (07/11/2023). ##
##    NOTES: "ghdr" refers to a genie_hdr dataset extracted from an hdf5      ##
##           file as seen in common/file_parsing.py                           ##
##           "traj" refers to trajectories dataset extracted from an hdf5     ##
##           file as seen in common/file_parsing.py                           ##
##                                                                            ##
################################################################################

import numpy as np
import geometry_methods
import particlePDG_defs
import singleParticleAssociation_methods as particle_assoc

####--------- DEFINITIONS RELATED TO FINDING KINEMATIC VARIABLES -----------####

threshold = 1.2 # stand-in threshold in cm (~3 pixels)


####----------------- TRUTH PRIMARY PARTICLE KINEMATICS --------------------####

def truth_primary_particle_kinetic_energy(pdg, track_id_set, vertex_assoc_traj, traj, ghdr):

    traj_id_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, vertex_assoc_traj, traj, ghdr)
    trackid_at_vertex_mask = vertex_assoc_traj['trackID']==traj_id_at_vertex
    track_at_vertex = vertex_assoc_traj[trackid_at_vertex_mask] 
    energy = track_at_vertex['E_start']
    ke = energy - particlePDG_defs.rest_mass_dict[pdg]

    return ke
    
    
def truth_primary_particle_momentum(track_id_set, vertex_assoc_traj, traj, ghdr):

    traj_id_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, vertex_assoc_traj, traj, ghdr)
    trackid_at_vertex_mask = vertex_assoc_traj['trackID']==traj_id_at_vertex
    track_at_vertex = vertex_assoc_traj[trackid_at_vertex_mask] 
    mom = np.sqrt(np.sum(track_at_vertex['pxyz_start']**2))

    return mom


####--------------------- ENERGY DEPOSITION (IN TPC) -----------------------####

def tpc_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_e={}
    for i in range(8): contained_e[i]=0.
    for sg in seg[seg_id_mask]:
        tpc_fv=geometry_methods.tpc_vertex([(sg['x_start']+sg['x_end'])/2.,
                                         (sg['y_start']+sg['y_end'])/2.,
                                         (sg['z_start']+sg['z_end'])/2.])
        for key in tpc_fv.keys():
            if tpc_fv[key]==True:
                contained_e[key]+=sg['dE']
    return contained_e


def tpc_contained_energy(pdg, trackID, traj, seg):
    if pdg in particlePDG_defs.neutral_pdg:
        temp={}
        for i in range(8): temp[i]=0.
        return temp
    else: return tpc_edep_charged_e(trackID, traj, seg)


####---------------- ENERGY DEPOSITION (FIDUCIAL VOLUME) -------------------####

def fv_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_e=0.
    for sg in seg[seg_id_mask]:
        if geometry_methods.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                              (sg['y_start']+sg['y_end'])/2.,
                                              (sg['z_start']+sg['z_end'])/2.]):
            contained_e+=sg['dE']
    return contained_e


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


def fv_contained_energy(pdg, trackID, traj, seg):
    if pdg in particlePDG_defs.neutral_pdg: return 0.#fv_edep_neutral_e(trackID, traj, seg)
    else: return fv_edep_charged_e(trackID, traj, seg)


####------------------ ENERGY DEPOSITION (TOTAL VOLUME) --------------------####

def total_edep_charged_e(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    total_e=0.
    for sg in seg[seg_id_mask]: total_e+=sg['dE']
    return total_e


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


def total_energy(pdg, trackID, traj, seg):
    if pdg in particlePDG_defs.neutral_pdg: return 0.#total_edep_neutral_e(trackID, traj, seg)
    else: return total_edep_charged_e(trackID, traj, seg)


####------------------- TRACK LENGTH (FIDUCIAL VOLUME) ---------------------####

def fv_edep_charged_length(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_length=0.
    for sg in seg[seg_id_mask]:
        if geometry_methods.fiducialized_vertex([(sg['x_start']+sg['x_end'])/2.,
                                              (sg['y_start']+sg['y_end'])/2.,
                                              (sg['z_start']+sg['z_end'])/2.]):
            contained_length+=np.sqrt( (sg['x_start']-sg['x_end'])**2+
                                       (sg['y_start']-sg['y_end'])**2.+
                                       (sg['z_start']-sg['z_end'])**2. )
    return contained_length


def fv_contained_length(pdg, trackID, traj, seg):
    if pdg in particlePDG_defs.neutral_pdg: return 0.
    else: return fv_edep_charged_length(trackID, traj, seg)


####--------------------- TRACK LENGTH (TOTAL VOLUME) ----------------------####


def total_edep_charged_length(trackID, traj, seg):
    seg_id_mask=seg['trackID']==trackID
    contained_length=0.
    for sg in seg[seg_id_mask]:
            contained_length+=np.sqrt( (sg['x_start']-sg['x_end'])**2+
                                       (sg['y_start']-sg['y_end'])**2.+
                                       (sg['z_start']-sg['z_end'])**2. )
    return contained_length


def total_length(pdg, trackID, traj, seg):
    if pdg in particlePDG_defs.neutral_pdg: return 0.
    else: return total_edep_charged_length(trackID, traj, seg)


####---------------------------- TRACK ANGLES ------------------------------####

def angle_wrt_beam_direction(trackid_set, vertex_assoc_traj,traj, ghdr):

    trackid_at_vertex = particle_assoc.find_trajectory_at_vertex(trackid_set, vertex_assoc_traj,traj, ghdr) # find track id for trajectory at vertex

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
