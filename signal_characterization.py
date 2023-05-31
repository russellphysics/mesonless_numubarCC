import matplotlib
import matplotlib.pyplot as plt
import twoBytwo_defs
import numpy as np
import auxiliary

'''
NOTE: While script is labelled as 'signal' characterization, methods can be used for
characterizing certain backgrounds as well.

INCLUDED METHODS:
 - get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, seg, signal_dict)
 - muon_characterization(spill_id, vert_id, ghdr, gstack, traj, seg, muon_dict, wrong_sign)
 - hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, seg, hadron_dict, wrong_sign)
 '''


''' Purpose: Fill Python dictionary with basic GENIE truth-level vertex information
             and 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), empty Python dictionary (DICT), 
             signifier of looking for right sign nu_mu_bar vertices or wrong sign nu_mu vertices (BOOL)
    Outputs: Nothing returned, but signal_dict (DICT) is full after
             method runs'''
def get_truth_dict(spill_id, vert_id, ghdr, gstack, traj, vert, seg, signal_dict):

    ghdr_vert_mask = ghdr['vertexID']==vert_id
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
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), empty Python dictionary (DICT), 
             signifier of looking for right sign nu_mu_bar vertices or wrong sign nu_mu vertices (BOOL)
    Outputs: Nothing returned, but muon_dict (DICT) is full after
             method runs'''
def muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, muon_dict, wrong_sign):

    traj_vert_mask = traj['vertexID']==vert_id 
    final_states = traj[traj_vert_mask] # Get trajectories associated with vertex

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask] # Get GENIE truth info associated with vertex

    mom = truth_level_summ['lep_mom'] # Truth-level outgoing muon momentum
    ang = truth_level_summ['lep_ang'] *np.pi / 180. # Truth-level muon angle with beam
    nu_energy = truth_level_summ['Enu'] # Truth-level neutrino energy
    q2 = truth_level_summ['Q2'] # Truth-level interaction 4-momentum squared
    nu_int_type = auxiliary.nu_int_type(ghdr, vert_id) # Neutrino interaction mechanism

    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0. # Set contained and total track energies and lengths to 0

    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg']) # Get set of PDG IDs for particles associated with vertex

    exclude_track_ids = set() # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    for fs in final_states:

        # Choose nu_mu_bar or nu_mu vertices
        if wrong_sign==False and (fs['pdgId'] != -13): continue
        elif wrong_sign==True and (fs['pdgId'] != 13): continue
        
        track_id = fs['trackID']
        if track_id in exclude_track_ids: continue

        pdg = fs['pdgId'] # *** pdg ***     
        parent_pdg = auxiliary.find_parent_pdg(fs['parentID'],vert_id, traj, ghdr)# *** parent pdg ***

        track_id_set = auxiliary.same_pdg_connected_trajectories(pdg, track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting

        is_primary = auxiliary.is_primary_particle(track_id_set, final_states, traj, ghdr, wrong_sign) 

        if is_primary == False: continue # Only look at final state particles
        for tid in track_id_set:

            total_edep += auxiliary.total_edep_charged_e(tid,traj,seg) # *** total visible energy ***
            contained_edep+= auxiliary.fv_edep_charged_e(tid, traj, seg)
            
            contained_length+=auxiliary.fv_edep_charged_length(tid, traj, seg) # *** total visible track length ***
            total_length+=auxiliary.total_edep_charged_length(tid, traj, seg)


    # Characterize Muon Endpoint/Containment
        end_pt = fs['xyz_end']
        start_pt = fs['xyz_start']
        end_pt_loc = twoBytwo_defs.particle_end_loc(start_pt, end_pt)

    # Save collected info in input muon_dict                                                                                                                         
    muon_dict[(spill_id,vert_id)]=dict(
        pdg=int(pdg),
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
        nu_int_type=str(nu_int_type))
    return


''' Purpose: Fill Python dictionary with basic GENIE and edep-sim truth-level vertex information
             related to hadrons 
    Inputs : Spill ID (INT), Vertex ID (INT), genie_hdr dataset (HDF5 DATASET), 
             genie_stack dataset (HDF5 DATASET), edep-sim trajectories dataset (HDF5 DATASET), 
             vertex dataset (HDF5 DATASET), edep-sim segements dataset (HDF5 DATASET), empty Python dictionary (DICT), 
             signifier of looking for right sign nu_mu_bar vertices or wrong sign nu_mu vertices (BOOL)
    Outputs: Nothing returned, but hadron_dict (DICT) is full after
             method runs'''
def hadron_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, hadron_dict, wrong_sign):
        
    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask] # Trajectories associated with vertex

    leptons_abs_pdg = [11, 12, 13, 14, 15, 16] # List of lepton PDG IDs (abs value)
    
    gstack_vert_mask = gstack['vertexID']==vert_id
    gstack_vert = gstack[gstack_vert_mask] # Particle ID information associated with vertex

    gstack_vert_fs_mask = gstack_vert['part_status']==1 # Excludes initial state particles
    gstack_vert_fs = gstack_vert[gstack_vert_fs_mask]['part_pdg'] # Final state particle PDG IDs

    gstack_vert_fs_hadrons = [fsp for fsp in gstack_vert_fs if abs(fsp) not in leptons_abs_pdg and fsp != 22] # LIST of f.s. hadron PDG IDs 
    gstack_vert_fs_pdg_set = set(gstack_vert_fs_hadrons) # SET of f.s. hadron PDG IDs

    #print("Vertex PDG Stack:", gstack_vert['part_pdg'])
    #print("Final State Hadrons:", gstack_vert_fs_hadrons)
    #print("Final State PDG Stack:", gstack_vert_fs)
    #print("Vertex PDG Stack Set:", gstack_pdg_set)

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

    # Set contained and total track energies and lengths to 0
    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.; max_proton_contained_length=0.; max_proton_total_length=0. 
    
    exclude_track_ids = set() # Create set of track IDs to exclude to eliminate redundancies
                              # (i.e. if they've been identified as being from the same particle as an earlier track)
    for fs in final_states:

        if abs(fs['pdgId']) in leptons_abs_pdg: continue # No leptons
        if fs['pdgId'] > 1000000000: continue # No nuclei
        if fs['pdgId'] == 22: continue # No photons
        
        track_id = fs['trackID']
        if track_id in exclude_track_ids: continue 

        track_id_set = auxiliary.same_pdg_connected_trajectories(fs['pdgId'], track_id, final_states, traj, ghdr)
        exclude_track_ids.update(track_id_set) # Exclude track IDs associated with same particle from future counting

        proton_contained_length = 0.; proton_total_length=0. # Reset proton track lengths
        for tid in track_id_set:

            total_edep += auxiliary.total_edep_charged_e(tid,traj,seg) # *** total visible energy ***
            contained_edep+= auxiliary.fv_edep_charged_e(tid, traj, seg)
            
            if fs['pdgId'] == 2212:
                proton_contained_length+=auxiliary.fv_edep_charged_length(tid, traj, seg) # *** total contained length for protons ***
                proton_total_length+=auxiliary.total_edep_charged_length(tid, traj, seg)

        if fs['pdgId'] == 2212 and proton_contained_length > max_proton_contained_length and \
            auxiliary.is_primary_particle(track_id_set, final_states, traj, ghdr, wrong_sign):
            max_proton_contained_length = proton_contained_length # Update max contained proton length in vertex
            max_proton_total_length = proton_total_length # Update max total proton length in vertex

    # Save collected info in input hadron_dict
    hadron_dict[(spill_id,vert_id)]=dict(
        hadron_mult = int(hadron_mult),
        neutron_mult = int(n_mult),
        proton_mult = int(p_mult),
        other_had_mult = int(other_had_mult),
        hadron_pdg = [int(i) for i in gstack_vert_fs_hadrons],
        hadron_pdg_set = [int(i) for i in list(gstack_vert_fs_pdg_set)],
        total_edep=float(total_edep),
        contained_edep=float(contained_edep),
        max_p_total_length=float(max_proton_total_length),
        max_p_contained_length=float(max_proton_contained_length))
    return
                                                 
