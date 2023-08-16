################################################################################
##                                                                            ##
##    CONTAINS: Methods related to ProtoDUNE-ND geometry.                     ##
##              Geometry definitions are copied straight from                 ##
## https://github.com/DUNE/2x2_sim/blob/main/validation/edepsim_validation.py ##    
##                                                                            ##
################################################################################

import numpy as np
import geometry_defs


####------------------ POSITION LOCATION CLASSIFICATION --------------------####

def fiducialized_particle_origin(traj, vert_id):
    traj_vert_mask = traj['vertex_id']==vert_id
    final_states = traj[traj_vert_mask]
    for fs in final_states:
        if fiducialized_vertex(fs['xyz_start'])==True:
            return True
    return False


def fiducialized_vertex(vert_pos):
    flag=False; x_drift_flag=False; y_vertical_flag=False; z_beam_flag=False
    for i in range(3):
        for i_bounds, bounds in enumerate(geometry_defs.tpc_bounds(i)):
            if vert_pos[i]>bounds[0] and vert_pos[i]<bounds[1]:
                if i==0: x_drift_flag=True; break
                if i==1: y_vertical_flag=True
                if i==2: z_beam_flag=True
    if x_drift_flag==True and y_vertical_flag==True and z_beam_flag==True: flag\
=True
    return flag


def tpc_vertex(vert_pos):
    temp=[]
    for i in range(3): temp.append(geometry_defs.tpc_bounds(i).tolist())
    tpc_fv={}
    for i in range(8): tpc_fv[i]=False
    tpc=0
    enclosed=False
    for x in range(4):
        for y in range(1):
            for z in range(2):
                if vert_pos[0]>temp[0][x][0] and vert_pos[0]<temp[0][x][1] and\
                   vert_pos[1]>temp[1][y][0] and vert_pos[1]<temp[1][y][1] and\
                   vert_pos[2]>temp[2][z][0] and vert_pos[2]<temp[2][z][1]:
                    tpc_fv[tpc]=True
                    return tpc_fv
                tpc+=1
    return tpc_fv


def minerva_vertex(vert_pos):
    upstream=False
    flag=False; x_flag=False; y_flag=False; z_upstream_flag=False; z_downstream_flag=False
    for i in range(3):
        ctr=0
        for i_bounds, bounds in enumerate(geometry_defs.MINERvA_bounds(i)):
            if i==0 and vert_pos[i]>bounds[0] and vert_pos[i]<bounds[1]:
                x_flag=True
            if i==1 and vert_pos[i]>bounds[0] and vert_pos[i]<bounds[1]:
                y_flag=True
            if i==2 and ctr==0 and \
               vert_pos[i]>bounds[0] and vert_pos[i]<bounds[1]:
                z_upstream_flag=True
            if i==2 and ctr==1 and \
               vert_pos[i]>bounds[0] and vert_pos[i]<bounds[1]:
                z_downstream_flag=True
            ctr+=1
    if x_flag==True and y_flag==True and z_upstream_flag==True:
        flag=True; upstream=True
    elif x_flag==True and y_flag==True and z_downstream_flag==True:
        flag=True; upstream=False
    return (flag, upstream)
                                          

####------------------ PARTICLE CONTAINMENT / ENDPOINTS --------------------####

def particle_containment(traj, traj_id):
    mask = traj['traj_id']==traj_id
    start=fiducialized_vertex(traj[mask]['xyz_start'][0].tolist())
    end=fiducialized_vertex(traj[mask]['xyz_end'][0].tolist())
    if start==True and end==True: return 'fc' # fully contained
    elif (start==True and end==False) or (start==False and end==True): return 'pc' # partially contained
    else: return 'tg' # through going


''' Inputs: ([x,y,z] vector) particle trajectory start point
            ([x,y,z] vector) particle trajectory end point
    Output: (string) tells where trajectory ends and/or exits detectors (key for dictionary)'''
def particle_end_loc(particle_start, particle_end):

    ## TO DO: add possibility of particle leaving from side or front of minerva upstream
    end_pt_loc = ''

    if fiducialized_vertex(particle_end):
        end_pt_loc = 'f'
    elif minerva_vertex(particle_end)[0]==True and minerva_vertex(particle_end)[1]==True:
        end_pt_loc = 'u'
    elif minerva_vertex(particle_end)[0]==True and minerva_vertex(particle_end)[1]==False:
        end_pt_loc = 'd'
    else:
        x_MINERvA = geometry_defs.MINERvA_bounds(0)[0]
        y_MINERvA = geometry_defs.MINERvA_bounds(1)[0]
        z_MINERvA_down = geometry_defs.MINERvA_bounds(2)[1]
        z_tpc_down = geometry_defs.tpc_bounds(2)[1]

        # Check whether endpoint Z is between 2x2 and MINERvA Downstream
        if particle_end[2]<z_MINERvA_down[0] and particle_end[2]>z_tpc_down[1]:
                end_pt_loc = 'p'

        # Check leaving back or side of MINERvA
        else:
            traj_vector = particle_end - particle_start

            if particle_end[2]>z_MINERvA_down[0] and particle_end[2]<z_MINERvA_down[1]:

                traj_param_front = (particle_end[2] - z_MINERvA_down[0])/traj_vector[2]
                minerva_down_front_intersect = particle_end + traj_param_front*traj_vector

                if minerva_down_front_intersect[0] > x_MINERvA[0] and minerva_down_front_intersect[0] < x_MINERvA[1] and\
                    minerva_down_front_intersect[1] > y_MINERvA[0] and minerva_down_front_intersect[1] < y_MINERvA[1]:

                    end_pt_loc = 's'

                else:
                    end_pt_loc = 'p'

            elif particle_end[2]>z_MINERvA_down[1]:

                traj_param_back = (z_MINERvA_down[1] - particle_end[2])/traj_vector[2]
                minerva_down_back_intersect = particle_end + traj_param_back*traj_vector

                traj_param_front = (z_MINERvA_down[0] - particle_end[2])/traj_vector[2]
                minerva_down_front_intersect = particle_end + traj_param_front*traj_vector

                if (minerva_down_back_intersect[2]-z_MINERvA_down[1]) > 0.01: 
                    print("STOP: MATH ERROR IN INTERSECTION CALCULATION!")
                    print('MINERvA back Z intersect:', round(minerva_down_back_intersect[2], 2))
                    print('MINERvA back Z:', z_MINERvA_down[1])
                if minerva_down_back_intersect[0] > x_MINERvA[0] and minerva_down_back_intersect[0] < x_MINERvA[1] and\
                    minerva_down_back_intersect[1] > y_MINERvA[0] and minerva_down_back_intersect[1] < y_MINERvA[1]:
                    
                    end_pt_loc = 'b'

                elif minerva_down_front_intersect[0] > x_MINERvA[0] and minerva_down_front_intersect[0] < x_MINERvA[1] and\
                    minerva_down_front_intersect[1] > y_MINERvA[0] and minerva_down_front_intersect[1] < y_MINERvA[1]:
                    
                    end_pt_loc = 's'
                
                else: 
                    end_pt_loc = 'p'

    return end_pt_loc