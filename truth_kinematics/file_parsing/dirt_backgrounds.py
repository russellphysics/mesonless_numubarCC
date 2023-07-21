import numpy as np
sys.path.append('../../common')
import geometry_methods as geo_methods

def dirt_muon_characterization(spill_id, vert_id, ghdr, gstack, traj, vert, seg, dirt_dict):

    traj_vert_mask = traj['vertexID']==vert_id
    final_states = traj[traj_vert_mask]

    ghdr_vert_mask = ghdr['vertexID']==vert_id
    truth_level_summ = ghdr[ghdr_vert_mask]


    total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.
    

    #print("PDG IDs of F.S. Particles:", final_states['pdgId'])
    #gstack_vert_mask = gstack['vertexID']==vert_id
    #gstack_pdg_set = set(gstack[gstack_vert_mask]['part_pdg'])
    #print("Event PDG Stack:", gstack_pdg_set)
    
    for fs in final_states:
        if fs['pdgId'] not in [13, -13]: continue # [111,211,-211]: continue

        #cut on mu length
        if np.linalg.norm(np.subtract(fs['xyz_end'],fs['xyz_start']))< 10: continue

        ##place a cut for backgrounds
        #print('muon start point', fs['xyz_start'])
        #print('muon end point', fs['xyz_end'])
        if( geo_methods.fiducialized_vertex(fs['xyz_start'])):
            mu_length = np.linalg.norm(np.subtract(fs['xyz_end'],fs['xyz_start']))
            #print('This is a dirt background')
            #print('Neutrino vertex: ', truth_level_summ['vertex'])
            #print('muon start point: ', fs['xyz_start'])
            #print('muon end point: ', fs['xyz_end'])
        else:
            #print('started out of FV')
            continue

        pdg = fs['pdgId'] # *** pdg ***
        
        parent_id = fs['parentID']
        parent = final_states['trackID']==parent_id

        #parent_pdg = final_states[parent]['pdgId'] # *** parent pdg ***

        if parent_id ==-1:
            print('mu parent was a neutrino')
            ghdr_mask = ghdr['vertexID']=fs['vertexID']
            parent_pdg = ghdr[ghdr_mask]['nu_pdg']
            parent_pdg = parent_pdg.tolist()[0]

        else: ##probably a meson
            parent_pdg = final_states[parent]['pdgId'] # *** parent pdg ***
            parent_pdg = parent_pdg.tolist()[0]
            parent_true_length= np.linalg.norm(np.subtract(final_states[parent]['xyz_end'],final_states[parent]['xyz_start']))
            if parent_true_length>3: continue
            
            print('mu start E: ', fs['E_start'])
            print('mu length: ', np.linalg.norm(np.subtract(fs['xyz_end'],fs['xyz_start'])))
            print('check parent pdg: ', parent_pdg)
            print('parent start E: ',final_states[parent]['E_start'])
            print('parent length [cm]: ', parent_true_length)
            #check grandparent                                                     
            grandparent_id = final_states[parent]['parentID']

            if grandparent_id ==-1:
                print('mu grandparent was a nutrino')
                continue
            else:
               grandparent_mask = final_states['trackID']==grandparent_id
               grandparent_pdg = final_states[grandparent_mask]['pdgId']
               grandparent_pdg = grandparent_pdg.tolist()[0]
               print('mu grandparent pdg:', grandparent_pdg)
               grandparent_length = np.linalg.norm(np.subtract(final_states[grandparent_mask]['xyz_end'],final_states[grandparent_mask]['xyz_start']))
               if grandparent_pdg == 2112:
                   print('This is very likely a dirt background')
                   print('neutron length: ', grandparent_length)
               else:
                   print('mu grandparent length: ', grandparent_length)
                   if grandparent_length>3: continue
                   print('need to check higher up in the family')
                   if final_states[grandparent_mask]['parentID'] == -1: continue
                   grandparent_mask = final_states['trackID']= final_states[grandparent_mask]['parentID']
                   grandparent_pdg = final_states[grandparent_mask]['pdgId']
                   grandparent_pdg = grandparent_pdg.tolist()[0]
                   print('mu grandparent2 pdg:', grandparent_pdg)
                   grandparent_length = np.linalg.norm(np.subtract(final_states[grandparent_mask]['xyz_end'],final_states[grandparent_mask]['xyz_start']))
                   print('mu grandparent2 length: ', grandparent_length)

                   if grandparent_length>3: continue
                   else:
                       print('still check higher up')


        track_id = fs['trackID']
        total_edep=0.; contained_edep=0.; total_length=0.; contained_length=0.

        #### checking dirt muons
        if abs(pdg)==13:
            seg_id_mask = seg['trackID']==track_id
            total_edep = sum(seg[seg_id_mask]['dE']) # *** total visible energy ***
            contained_edep = 0
            contained_length = 0
            
            for sg in seg[seg_id_mask]:
                total_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                      (sg['y_start']-sg['y_end'])**2+
                                      (sg['z_start']-sg['z_end'])**2) # *** total length *** This is not correct

            for sg in seg[seg_id_mask]:
                if geo_methods.fiducialized_vertex( [(sg['x_start']+sg['x_end'])/2.,
                                         (sg['y_start']+sg['y_end'])/2.,
                                         (sg['z_start']+sg['z_end'])/2.] ):
                    contained_edep+=sg['dE'] # *** contained visible energy ***

                    contained_length+=np.sqrt((sg['x_start']-sg['x_end'])**2+
                                              (sg['y_start']-sg['y_end'])**2+
                                              (sg['z_start']-sg['z_end'])**2) # *** contained length ***

        if contained_edep>5:
            print(pdg,'\t',parent_pdg,'\t',total_edep,' MeV\t',contained_edep,' MeV\t', total_length,' cm\t',contained_length,' cm')


            dirt_dict[(spill_id,vert_id, track_id)]=dict(
                pdg=int(pdg),
                parent_pdg=parent_pdg,
                total_edep=total_edep,
                contained_edep=contained_edep,
                total_length=total_length,
                contained_length=contained_length,
                true_mom=truth_level_summ['lep_mom'],
                true_angle= truth_level_summ['lep_ang'],
                true_energy=truth_level_summ['Elep'],
                nu_energy=truth_level_summ['Enu'],
                q_sq = truth_level_summ['Q2'])
            
            ## end if
    return


