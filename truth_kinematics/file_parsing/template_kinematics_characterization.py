import h5py
import argparse
import numpy as np
import sys
sys.path.append('../../common')
import file_parsing
import geometry_defs as geo_defs
import geometry_methods as geo_methods
import particlePDG_defs as pdg_defs
import truth_methods as truth
import singleParticleAssociation_methods as particle_assoc
import kinematicVariable_methods as kinematics



def main(sim_file, input_type):

    sim_h5 = h5py.File(sim_file,'r')
#    print_keys_attributes(sim_h5)
#    return

    ### partition file by spill
    unique_spill = np.unique(sim_h5['trajectories']['spillID'])
    for spill_id in unique_spill:

        ghdr, gstack, traj, vert, seg = file_parsing.get_spill_data(sim_h5, spill_id)

        ### partition by vertex ID within beam spill
        for v_i in range(len(vert['vertexID'])):
            vert_pos= [vert['x_vert'][v_i], vert['y_vert'][v_i], vert['z_vert'][v_i]]
            vert_in_active_LAr = geo_methods.fiducialized_vertex( vert_pos )

            nu_mu_bar = truth.signal_nu_pdg(ghdr, vert['vertexID'][v_i])
            is_cc = truth.signal_cc(ghdr, vert['vertexID'][v_i])
            pionless = truth.signal_meson_status(gstack, vert['vertexID'][v_i])
            fv_particle_origin=geo_methods.fiducialized_particle_origin(traj, vert['vertexID'][v_i])


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--sim_file', default=None, required=True, type=str, help='''string corresponding to the path of the edep-sim ouput simulation file to be considered''')
    parser.add_argument('-t', '--input_type', default='edep', choices=['edep', 'larnd'], type=str, help='''string corresponding to the output file type: edep or larnd''')
    args = parser.parse_args()
    main(**vars(args))
