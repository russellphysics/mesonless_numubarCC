################################################################################
##                                                                            ##
##    CONTAINS: Methods to create stacked histograms describing kinematics    ##
##              of signal and all background events using previously created  ##
##              JSON dictionaries.                                            ##
##                                                                            ##
################################################################################

import numpy as np
import matplotlib.pyplot as plt
import json
import argparse



def files_processed(processed_files, total_files=1023, \
                    production_pot=1e19, target_pot=2.5e19):
    return target_pot/((processed_files*production_pot)/total_files)



def plot_stacked_histo(signal, signal_factor, \
                       cc_threshold, cc_threshold_factor, \
                       nc_pid, nc_pid_factor, dirt, dirt_factor, \
                       metric, bins, xlabel, ylabel, figname, leg_location, \
                       xlim,yscale):
    fig, ax = plt.subplots(figsize=(6,6))

    s = [signal[key][metric]/1e3 for key in signal.keys()]
    d = [dirt[key][metric]/1e3 for key in dirt.keys()]
    t = [cc_threshold[key][metric]/1e3 for key in cc_threshold.keys()]
    p = [nc_pid[key][metric]/1e3 for key in nc_pid.keys()]
    if metric=='q2':
        s = [signal[key][metric]/1e6 for key in signal.keys()]
        d = [dirt[key][metric]/1e6 for key in dirt.keys()]
        t = [cc_threshold[key][metric]/1e6 for key in cc_threshold.keys()]
        p = [nc_pid[key][metric]/1e6 for key in nc_pid.keys()]        
        
    s_weight = [signal_factor]*len(s)
    d_weight = [dirt_factor]*len(d)
    t_weight = [cc_threshold_factor]*len(t)
    p_weight = [nc_pid_factor]*len(p)

    print('signal: ',len(s)*signal_factor,'\n',
          'dirt background: ',len(d)*dirt_factor,'\n',
          'threshold background: ',len(t)*cc_threshold_factor,'\n',
          'pid background: ',len(p)*nc_pid_factor)

    ax.hist(s, bins=bins, weights=s_weight, stacked=True, histtype='bar',\
            label=r'mesonless $\bar{\nu}_\mu$ CC')
    ax.hist(p, bins=bins, weights=p_weight,  stacked=True, histtype='bar',\
            label=r'$\nu$NC pid backgrounds')
    ax.hist(t, bins=bins, weights=t_weight,  stacked=True, histtype='bar',\
            label=r'$\nu$CC threshold backgrounds')
    ax.hist(d, bins=bins, weights=d_weight,  stacked=True, histtype='bar',\
            label=r'dirt backgrounds')
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(xlim)
    ax.set_yscale(yscale)
    ax.legend(loc=leg_location)
    ax.grid(True)
    plt.savefig(figname+'.png')
    


def main(signal, n_signal, dirt, n_dirt, \
         cc_threshold, n_cc_threshold, nc_pid, n_nc_pid):
    f = open(signal)
    signal_dict=json.load(f)
    signal_sf = files_processed(n_signal)

    f = open(dirt)
    dirt_dict=json.load(f)
    dirt_sf = files_processed(n_dirt)

    f = open(cc_threshold)
    cc_threshold_dict=json.load(f)
    cc_threshold_sf = files_processed(n_cc_threshold)

    f = open(nc_pid)
    nc_pid_dict=json.load(f)
    nc_pid_sf = files_processed(n_nc_pid)                     

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'nu_energy', np.linspace(0,10,21), \
                       r'$\nu$ Energy [GeV]', r'$\nu$ Interactions / 500 MeV',\
                       'stacked_nu_energy', 'center right',(0,10),'linear')

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'q2', np.linspace(0,5,26), \
                       r'$Q^2$ [GeV$^2$]', r'$\nu$ Interactions / 200 MeV$^2$',\
                       'stacked_q2', 'upper right',(0,5),'linear')

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'mom', np.linspace(0,10,21), \
                       r'Muon Candidate Momentum [GeV/c]', \
                       r'$\nu$ Interactions / 500 MeV/c',\
                       'stacked_mu_momentum', 'center right',(0,10),'linear')

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'ang', np.linspace(-0.45,0.75,61), \
                       r'$\theta_\mu$ [radians]', r'$\nu$ Interactions / 0.02 radians',\
                       'stacked_mu_angle', 'upper right',(-0.45,0.75),'linear')

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'vtx_x', np.linspace(-600,600,51), \
                       r'$\nu$ Vertex X Position [cm]', r'$\nu$ Interactions',\
                       'stacked_vertex_x', 'upper right',(-600,600),'linear')

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'vtx_y', np.linspace(-500,500,51), \
                       r'$\nu$ Vertex Y Position [cm]', r'$\nu$ Interactions',\
                       'stacked_vertex_y', 'upper right',(-500,500),'linear')

    plot_stacked_histo(signal_dict, signal_sf, \
                       cc_threshold_dict, cc_threshold_sf, \
                       nc_pid_dict, nc_pid_sf, \
                       dirt_dict, dirt_sf, \
                       'vtx_z', np.linspace(-2000,1000,51), \
                       r'$\nu$ Vertex Z Position [cm]', r'$\nu$ Interactions',\
                       'stacked_vertex_z', 'upper right',(-2000,1000),'linear')


    
if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--signal', default='signal_dict.json', \
                        type=str, help='''signal JSON''')
    parser.add_argument('-ns','--n_signal', default=50, type=int, \
                        help='''number of files processed for signal JSON''')
    parser.add_argument('-cc','--cc_threshold', default='cc_threshold_bkg_dict.json', \
                        type=str, help='''nuCC threshold background JSON''')
    parser.add_argument('-ncc','--n_cc_threshold', default=100, type=int, \
                        help='''number of files processed for nuCC threshold background JSON''')
    parser.add_argument('-nc','--nc_pid', default='nc_pid_bkg_dict.json', \
                        type=str, help='''nuNC PID background JSON''')
    parser.add_argument('-nnc','--n_nc_pid', default=100, type=int, \
                        help='''number of files processed for nuNC pid background JSON''')
    parser.add_argument('-d','--dirt', default='dirt_bkg_dict.json', \
                        type=str, help='''dirt background JSON''')
    parser.add_argument('-nd','--n_dirt', default=1000, type=int, \
                        help='''number of files processed for dirt background JSON''')
    args = parser.parse_args()
    main(**vars(args))
