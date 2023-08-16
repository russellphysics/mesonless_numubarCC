################################################################################
##                                                                            ##
##    CONTAINS: Methods related to parsing and extracting particle and event  ##
##              truth information.                                            ##
##    NOTE: "ghdr" refers to a genie_hdr dataset extracted from an hdf5       ##
##          file as seen in common/file_parsing.py                            ##
##                                                                            ##
################################################################################

import particlePDG_defs


####---------------------- SIGNAL BASIC CHARACTERISTICS --------------------####

def signal_nu_pdg(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    ghdr_nu_interaction = ghdr[ghdr_vert_mask]['nu_pdg']
    if abs(ghdr_nu_interaction[0])==particlePDG_defs.nu_mu_pdg: return True
    else: return False


def signal_cc(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertex_id']==vert_id
    return ghdr[ghdr_vert_mask]['isCC'][0]


def signal_meson_status(gstack, vert_id):
    gstack_vert_mask = gstack['vertex_id']==vert_id
    gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    if len(particlePDG_defs.meson_pdg.intersection(gstack_pdg_set))==0: return True
    else: return False


def nu_int_type(ghdr, vert_id):
    ghdr_vert_mask = ghdr['vertex_id']==vert_id
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

####--------------------- PARTICLE ORIGIN INFORMATION ----------------------####

def is_primary_particle(traj_id_set, vertex_assoc_traj,traj, ghdr):
    is_prim = False

    for tid in traj_id_set:

        particle_mask = vertex_assoc_traj['traj_id'] == tid
        parent_pdg = find_parent_pdg(vertex_assoc_traj[particle_mask]['parent_id'],
                                     vertex_assoc_traj[particle_mask]['vertex_id'],
                                     traj, ghdr)
        #print("Parent PDG:", parent_pdg)
        if abs(parent_pdg)==14:
            is_prim = True
            break
        else: continue

    return is_prim


def find_parent_pdg(parent_id, vertex_id, traj, ghdr):
    if parent_id==-1:
        ghdr_mask = ghdr['vertex_id']==vertex_id
        parent_pdg=ghdr[ghdr_mask]['nu_pdg']
    else:
        parent_mask = traj['traj_id']==parent_id
        parent_pdg = traj[parent_mask]['pdg_id']
    if parent_pdg==[]: parent_pdg=[0]
    return parent_pdg

