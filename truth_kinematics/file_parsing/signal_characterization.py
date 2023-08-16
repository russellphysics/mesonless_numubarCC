################################################################################
##                                                                            ##
##    CONTAINS: Methods to create python dictionaries with information to     ##
##              describe general signal events (get_signal_dict, muons in     ##
##              signal events (muon_characterization), and hadrons in signal  ##
##              events (hadron_characterization).                             ##
##                                                                            ##
################################################################################

import sys
sys.path.append('../../common')
import geometry_methods as geo_methods
import particlePDG_defs as pdg_defs
import truth_methods as truth
import singleParticleAssociation_methods as particle_assoc
import kinematicVariable_methods as kinematics
import numpy as np

''' TO DO: Add other hadron mult over threshold? '''


'''
NOTE: While script is labelled as 'signal' characterization, methods can be used for
characterizing certain backgrounds as well.

INCLUDED METHODS:
 - get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, seg, signal_dict)
 - muon_characterization(spill_id, vert_id, ghdr, gstack, traj, seg, muon_dict)
 - hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, seg, hadron_dict)
 '''


''' Purpose: Fill Python dictionary with basic GENIE truth-level vertex information
             and 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), 
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but signal_dict (DICT) is full after
             method runs'''
def get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict):

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]

    mom = truth_level_summ['lep_mom'] # Truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    vtx = truth_level_summ['vertex'] # Truth-level vertex information
    vtx_x = vtx[0][0]
    vtx_y = vtx[0][1]
    vtx_z = vtx[0][2]
                                                                                                                                          
    signal_dict[(spill_id,vert_id)]=dict(
        nu_energy=float(nu_energy),
        q2 = float(q2),
        mom=float(mom), 
        ang=float(ang),
        vtx_x = float(vtx_x), 
        vtx_y = float(vtx_y), 
        vtx_z = float(vtx_z))
    return


''' Purpose: Fill Python dictionary with basic GENIE and edep-sim truth-level vertex information
             related to FS muon
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), 
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but muon_dict (DICT) is full after
             method runs'''
def muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict):

    traj_vert_mask = traj['vertex_id']==vert_id 
    final_states = traj[traj_vert_mask] # Get trajectories associated with vertex

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask] # Get GENIE truth info associated with vertex

    mom = truth_level_summ['lep_mom'] # Truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    nu_int_type = truth.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism
    nu_pdg = truth_level_summ['nu_pdg']

    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0. # Set contained and total track energies and lengths to 0

    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg']) # Get set of PDG IDs for particles associated with vertex

    exclude_track_ids = set() # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    for fs in final_states:

        # Choose nu_mu_bar or nu_mu vertices
        if (abs(fs['pdg_id']) != 13): continue
        
        track_id = fs['traj_id']
        if track_id in exclude_track_ids: continue

        pdg = fs['pdg_id'] # *** pdg ***     

        track_id_set = particle_assoc.same_pdg_connected_trajectories(pdg, track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting

        is_primary = truth.is_primary_particle(track_id_set, final_states, traj, ghdr) 

        if is_primary == False: continue # Only look at final state particles

        track_id_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states,traj, ghdr)
        final_state_vertex_tid_mask = final_states['traj_id'] == track_id_at_vertex
        fs_at_vertex = final_states[final_state_vertex_tid_mask]

        parent_pdg = truth.find_parent_pdg(fs_at_vertex['parent_id'],vert_id, traj, ghdr)# *** parent pdg ***

        if len(track_id_set)>1:
            print("Length of Track ID Set:", len(track_id_set))

        for tid in track_id_set:

            total_edep += kinematics.total_edep_charged_e(tid,traj,seg) # *** total visible energy ***
            contained_edep+= kinematics.fv_edep_charged_e(tid, traj, seg)
            
            contained_length+=kinematics.fv_edep_charged_length(tid, traj, seg) # *** total visible track length ***
            total_length+=kinematics.total_edep_charged_length(tid, traj, seg)
        
        # Characterize Muon Endpoint/Containment
        if len(track_id_set)>1:
            start_tid = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states,traj, ghdr)
            start_tid_mask = final_states['traj_id']==start_tid
            muon_start_traj = final_states[start_tid_mask]
            start_pt = muon_start_traj['xyz_start']
            
            end_tid = particle_assoc.find_forward_primary_particle_end_trajectory(track_id_set, final_states,traj, ghdr)
            end_tid_mask = final_states['traj_id']==end_tid
            muon_end_traj = final_states[end_tid_mask]
            end_pt = muon_end_traj['xyz_end']
        else:
            end_pt = fs['xyz_end']
            start_pt = fs['xyz_start']
        end_pt_loc = geo_methods.particle_end_loc(start_pt, end_pt)

    # Save collected info in input muon_dict                                                                                                                         
    muon_dict[(spill_id,vert_id)]=dict(
        pdg=int(pdg),
        nu_pdg = int(nu_pdg), 
        parent_pdg=int(parent_pdg),
        total_edep=float(total_edep),
        contained_edep=float(contained_edep),
        total_length=float(total_length),
        contained_length=float(contained_length),
        mom=float(mom), 
        ang=float(ang),
        nu_energy=float(nu_energy),
        q2 = float(q2),
        end_pt_loc = str(end_pt_loc),
        muon_start = [float(i) for i in start_pt],
        muon_end = [float(i) for i in end_pt],
        nu_int_type=str(nu_int_type))
    return


''' Purpose: Fill Python dictionary with basic GENIE and edep-sim truth-level vertex information
             related to hadrons 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), threshold length in cm (FLOAT)
             empty Python dictionary (DICT)
    Outputs: Nothing returned, but hadron_dict (DICT) is full after
             method runs'''
def hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, threshold, hadron_dict):
        
    #print("\nHADRONS:")
    traj_vert_mask = traj['vertex_id']==vert_id
    final_states = traj[traj_vert_mask] # Trajectories associated with vertex

    leptons_abs_pdg = [11, 12, 13, 14, 15, 16] # List of lepton PDG IDs (abs value)

    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask] # Get GENIE truth info associated with vertex
    nu_pdg = truth_level_summ['nu_pdg']
    
    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_vert = gstack[gstack_vert_mask] # Particle ID information associated with vertex

    gstack_vert_fs_mask = gstack_vert['part_status']==1 # Excludes initial state particles
    gstack_vert_fs = gstack_vert[gstack_vert_fs_mask]['part_pdg'] # Final state particle PDG IDs

    gstack_vert_fs_hadrons = [fsp for fsp in gstack_vert_fs if abs(fsp) not in leptons_abs_pdg and fsp != 22] # LIST of f.s. hadron PDG IDs 
    gstack_vert_fs_pdg_set = set(gstack_vert_fs_hadrons) # SET of f.s. hadron PDG IDs

    #print("Vertex PDG Stack:", gstack_vert['part_pdg'])
    #print("Final State Hadrons:", gstack_vert_fs_hadrons)
    #print("Final State PDG Stack:", gstack_vert_fs)
    #print("Vertex PDG Stack Set:", gstack_pdg_set)
    nu_int_type = truth.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism (Truth)

    hadron_mult = len(gstack_vert_fs_hadrons) # F.S. hadron multiplicity
    n_mult = 0
    p_mult = 0
    other_had_mult = 0

    # Get f.s. hadron sub-category multiplicities
    for had in range(hadron_mult):

        if gstack_vert_fs_hadrons[had] == 2112:
            n_mult+=1 # Neutron multiplicity
        elif gstack_vert_fs_hadrons[had] == 2212:
            p_mult+=1 # Proton multiplicity
        else:
            other_had_mult+=1

    # Set contained and total track energies and lengths to 0; Set hadron + proton multiplicities over threshold to 0.
    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    total_edep_over_thresh=0.; contained_edep_over_thresh=0.
    max_proton_contained_length=0.; max_proton_total_length=0. 
    lead_proton_ang_wrt_beam=0.; lead_proton_momentum=0.; sub_lead_proton_ang_wrt_beam=0.; sub_lead_proton_momentum=0.
    lead_proton_traj_at_vertex = 0.; sub_lead_proton_traj_at_vertex = 0.
    hadron_mult_over_thresh = 0.; p_mult_over_thresh = 0.
    p_ke = 0.

    exclude_track_ids = set() # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    for fs in final_states:

        if abs(fs['pdg_id']) in leptons_abs_pdg: continue # No leptons
        if fs['pdg_id'] > 1000000000: continue # No nuclei
        if fs['pdg_id'] == 22: continue # No photons
        
        track_id = fs['traj_id']
        if track_id in exclude_track_ids: 
            #print("Excluding a track because it has already been studied.")
            continue 

        track_id_set = particle_assoc.same_pdg_connected_trajectories(fs['pdg_id'], track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting
        #if fs['pdg_id'] == 2112: print("\nTrack ID Set:", track_id_set)

        proton_contained_length = 0.; proton_total_length=0. # Reset proton track lengths
        fs_total_edep =0.; fs_contained_edep=0. # Reset Edep for individual final states
        for tid in track_id_set:

            total_edep_temp = kinematics.total_edep_charged_e(tid,traj,seg)
            contained_edep_temp = kinematics.fv_edep_charged_e(tid, traj, seg)
            fs_total_edep+=total_edep_temp
            total_edep += total_edep_temp# *** total visible energy ***
            fs_contained_edep+= contained_edep_temp
            contained_edep+= contained_edep_temp

            contained_length+=kinematics.fv_edep_charged_length(tid, traj, seg) # *** total contained length for protons ***
            total_length+=kinematics.total_edep_charged_length(tid, traj, seg)
            
            if fs['pdg_id'] == 2212:
                proton_contained_length+=kinematics.fv_edep_charged_length(tid, traj, seg) # *** total contained length for protons ***
                proton_total_length+=kinematics.total_edep_charged_length(tid, traj, seg)

        if truth.is_primary_particle(track_id_set, final_states, traj, ghdr) and contained_length > threshold \
            and fs['pdg_id'] not in pdg_defs.neutral_hadron_pdg_dict.keys():
            hadron_mult_over_thresh +=1
            if fs['pdg_id'] == 2212: 
                p_mult_over_thresh += 1
                p_traj_id_at_vertex = particle_assoc.find_trajectory_at_vertex(track_id_set, final_states, traj, ghdr)
                p_mom = kinematics.truth_primary_particle_momentum(track_id_set, final_states, traj, ghdr)
                p_ke += kinematics.truth_primary_particle_kinetic_energy(fs['pdg_id'],track_id_set, final_states, traj, ghdr)
                total_edep_over_thresh += fs_total_edep # *** total visible energy ***
                contained_edep_over_thresh+= fs_contained_edep
                if p_mom > lead_proton_momentum:
                    sub_lead_proton_traj_at_vertex = lead_proton_traj_at_vertex
                    sub_lead_proton_momentum = lead_proton_momentum
                    sub_lead_proton_ang_wrt_beam = lead_proton_ang_wrt_beam

                    lead_proton_traj_at_vertex = p_traj_id_at_vertex
                    lead_proton_momentum = p_mom
                    lead_proton_ang_wrt_beam = kinematics.angle_wrt_beam_direction(track_id_set, final_states, traj, ghdr)
                elif p_mom <= lead_proton_momentum and p_mom > sub_lead_proton_momentum:
                    sub_lead_proton_traj_at_vertex = p_traj_id_at_vertex
                    sub_lead_proton_momentum = p_mom
                    sub_lead_proton_ang_wrt_beam = kinematics.angle_wrt_beam_direction(track_id_set, final_states, traj, ghdr)
                
                if proton_contained_length > max_proton_contained_length:
                    max_proton_contained_length = proton_contained_length # Update max contained proton length in vertex
                    max_proton_total_length = proton_total_length # Update max total proton length in vertex

    angle_between_lead_and_sublead_protons = 0. # Initialize angle between leading and subleading protons
    if p_mult_over_thresh >= 2:
        angle_between_lead_and_sublead_protons = kinematics.angle_between_two_trajectories(lead_proton_traj_at_vertex, \
                                                                                          sub_lead_proton_traj_at_vertex, final_states)

    # Save collected info in input hadron_dict
    hadron_dict[(spill_id,vert_id)]=dict(
        nu_pdg = int(nu_pdg),
        hadron_mult = int(hadron_mult),
        neutron_mult = int(n_mult),
        proton_mult = int(p_mult),
        other_had_mult = int(other_had_mult),
        hadron_mult_over_thresh = int(hadron_mult_over_thresh),
        proton_mult_over_thresh = int(p_mult_over_thresh),
        hadron_pdg = [int(i) for i in gstack_vert_fs_hadrons],
        hadron_pdg_set = [int(i) for i in list(gstack_vert_fs_pdg_set)],
        total_edep=float(total_edep),
        contained_edep=float(contained_edep),
        max_p_total_length=float(max_proton_total_length),
        max_p_contained_length=float(max_proton_contained_length),
        lead_proton_momentum = float(lead_proton_momentum), 
        sub_lead_proton_momentum = float(sub_lead_proton_momentum), 
        lead_proton_ang_wrt_beam = float(lead_proton_ang_wrt_beam), 
        sub_lead_proton_ang_wrt_beam = float(sub_lead_proton_ang_wrt_beam),
        sub_lead_proton_angle_with_lead_proton = float(angle_between_lead_and_sublead_protons),
        primary_protons_total_ke = float(p_ke),
        total_edep_over_thresh = float(total_edep_over_thresh),
        contained_edep_over_thresh = float(contained_edep_over_thresh),
        nu_int_type=str(nu_int_type))
    return
                                                 
