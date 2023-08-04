################################################################################
##                                                                            ##
##    CONTAINS: Methods to create a JSON dictionary describing NC PID         ##
##              background events and certain associated plots.               ##
##                                                                            ##
################################################################################

import matplotlib.pyplot as plt
import argparse
import numpy as np
import json
from mpl_toolkits.axes_grid1.inset_locator import (inset_axes, InsetPosition, mark_inset)
import sys
sys.path.append('../../common')
import file_parsing


def files_processed(processed_files, total_files=1023, \
                    production_pot=1e19, target_pot=2.5e19):
    return target_pot/((processed_files*production_pot)/total_files)



def charged_pion_threshold(d, threshold, scale_factor):
    background_dict=dict()
    cp_count={} # charged pion count above tracking threshold per neutrino interaction
    for key in d.keys():
        if d[key]['pdg']==abs(211):
            contained_length=d[key]['contained_length']
            if contained_length>threshold:
                spill_id = key.split('-')[0]
                vertex_id = key.split('-')[1]
                temp=(spill_id, vertex_id)
                if temp not in cp_count: cp_count[temp]=0.
                cp_count[temp]+=1

    cp_length={} # single track MIP length per neutrino interaction
    for key in d.keys():
        spill_id = key.split('-')[0]
        vertex_id = key.split('-')[1]
        temp=(spill_id, vertex_id)
        if temp in cp_count.keys():
            if cp_count[temp]!=1: continue
            contained_length=d[key]['contained_length']
            cp_length[temp]=contained_length            

    count_bgd=0
    for key in cp_length.keys():
        if cp_length[key]>threshold: count_bgd+=1
    print(count_bgd,' irreducible NC backgrounds with ',threshold,' cm tracking threshold\n',count_bgd*scale_factor,' irreducible events in 2.5x10^19 POT')

    candidate={}
    for key in d.keys():
        spill_id = key.split('-')[0]
        vertex_id = key.split('-')[1]
        temp=(spill_id, vertex_id)
        if temp in cp_length.keys():
            if cp_length[temp]>threshold:
                print("Candidate keys:", candidate.keys())
                print(d[key])
                ##### end_pt_loc not yet implemented here, so the following lines break the code:
                #if d[key]['end_pt_loc'] not in candidate.keys():
                #    candidate[d[key]['end_pt_loc']]=[]
                #candidate[d[key]['end_pt_loc']].append(d[key]['mom'])
                    
                background_dict[temp]=dict(
                    nu_energy=d[key]['nu_energy'],
                    q2=d[key]['q2'],
                    mom=d[key]['mom'],
                    ang=d[key]['ang'],
                    vtx_x=d[key]['vtx_x'],
                    vtx_y=d[key]['vtx_y'],
                    vtx_z=d[key]['vtx_z']
                    )

    bins=np.linspace(0,200,41)
    fig, ax = plt.subplots(figsize=(6,6))

    meson_length=[cp_length[key] for key in cp_length.keys()]
    weight=[scale_factor]*len(meson_length)
    ax.hist(meson_length, bins=bins, weights=weight, histtype='step', color='b')
    ax.axvline(x=threshold, linestyle='dashed', color='orange', label='3 cm tracking threshold')
    ax.set_xlim(0,200)
    ax.legend(loc='lower right')
    ax.set_xlabel(r'Maximum 2x2-contained $\pi^{\pm}$ Track Length [cm]')
    ax.set_ylabel(r'$\nu$ Interactions / 5 cm')
    
    axtwin = ax.twinx()
    axtwin.hist(meson_length, bins=bins, weights=weight, histtype='step',\
                cumulative=True, density=True, linestyle='solid', color='r')
    axtwin.set_ylabel('Cumulative Probability')

    ax.text(115,1670,r'NuMI ME RHC 2.5$\times$10$^{19}$ POT')
    ax.text(0,1670,r'$\nu$NC')

    plt.savefig('irreducible_nc_beam_bgd.png')

    return background_dict

    

def main(nc_json_file, tracking_threshold, n_files_processed):
    f = open(nc_json_file)
    nc_pion_dict=json.load(f)
    scale_factor = files_processed(n_files_processed)
    bgd_dict = charged_pion_threshold(nc_pion_dict, tracking_threshold, \
                                      scale_factor)
    file_parsing.save_dict_to_json(bgd_dict, 'nc_pid_bkg_dict', True)
    


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-nc', '--nc_json_file', default=None, required=False, type=str, \
                        help='''string corresponding to the path of the NC pion backgrounds JSON file''')
    parser.add_argument('-t', '--tracking_threshold', default=3., required=False, type=float, \
                        help='''Tracking threshold in track length [cm]''')
    parser.add_argument('-n', '--n_files_processed', default=1, required=True, type=int, \
                        help='''File count of number of files processed in production sample''')
    args = parser.parse_args()
    main(**vars(args))
