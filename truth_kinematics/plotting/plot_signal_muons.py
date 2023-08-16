################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots describing muons in signal             ##
##              events using a muon dictionary created using methods in       ##
##              /truth_kinematics/file_parsing/signal_characterization.py     ##
##              and a scale factor for scaling event counts to those expected ##
##              with 2.5e19 POT.                                              ##
##                                                                            ##
################################################################################

import matplotlib.pyplot as plt
import numpy as np
import sys
sys.path.append('../file_parsing')
sys.path.append('../../common')
import geometry_defs as geo_defs

# PLOT: Muon kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_muons(d, scale_factor, sig_bkg = 0):
    
    # DEFINE: Plotting muon kinematics for signal or background events
    sample_type = ''
    sample_title = ''
    if sig_bkg == 0: 
        sample_type = 'signal'
        sample_title = 'Signal'
    elif sig_bkg == 1:
        sample_type = 'dirt_bkg'
        sample_title = 'Dirt Background'
    elif sig_bkg == 2:
        sample_type = 'beam_bkg'
        sample_title = 'Beam Background'
    else: 
        return "Error: plot_muons function given undefined signal/background definition"
    
        
                                                  
    # PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(figsize=(8,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    #data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,400,20))
    #counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,400,20))
    ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='Total', histtype='step')
    #ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='Contained',histtype='step', linestyle='--')
    ax0.set_xlabel('Total Visible Muon Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    #ax0.set_title(r'Muon Energy')
    #ax0.legend()
    ax0.set_yscale('log')
    ax0.grid(True)
    plt.savefig(sample_type+"_events_muon_visible_energy.png")
    plt.close(fig0)

    # PLOT: muon energy containment fraction ## UNNECESSARY PLOT AT THE MOMENT
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    counts1, bins1 = np.histogram(data1, bins=np.linspace(0,1,20))
    ax1.hist(bins1[:-1], bins=bins1, weights = counts1*scale_factor, histtype='step')
    ax1.set_xlabel('Visible Muon Energy Containment Fraction')
    ax1.set_ylabel('Count / 0.05')
    ax1.grid(True)     
    ax0.set_yscale('log')  
    plt.savefig(sample_type+"_events_muon_containment_fraction.png")
    plt.close(fig1)    
    
    # PLOT: truth-level outgoing muon (lepton) angle 
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.array([d[key]['ang'] for key in d.keys()])
    counts2, bins2 =np.histogram(data2, bins=np.linspace(0.00,0.80,41))
    ax2.hist(bins2[:-1], bins=bins2, weights = counts2*scale_factor, histtype='step')
    ax2.set_xlabel(r"Outgoing Muon Angle with Beam Direction [Rad]")
    ax2.set_ylabel("Count / 0.02 Rad")
    plt.savefig(sample_type+"_events_outgoing_muon_angle_truth.png")
    plt.close(fig2)   

    # PLOT: truth-level outgoing muon (lepton) momentum 
    fig3, ax3 = plt.subplots(figsize=(6,4))
    data3 = np.array([d[key]['mom'] for key in d.keys()])/1000.
    counts3, bins3 =np.histogram(data3, bins=np.linspace(0,15,31))
    ax3.hist(bins3[:-1], bins=bins3, weights = counts3*scale_factor, histtype='step')
    ax3.set_xlabel(r"Outgoing Muon Momentum [GeV/c]")
    ax3.set_ylabel("Count / 0.5 GeV/c")
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth.png")
    plt.close(fig3)  

    # PLOT: truth-level 4-momentum squared of interaction
    fig4, ax4 = plt.subplots(figsize=(6,4))
    data4 = np.array([d[key]['q2'] for key in d.keys()]) / 1000000.
    counts4, bins4 =np.histogram(data4, bins=np.linspace(0,5,51))
    ax4.hist(bins4[:-1], bins=bins4, weights = counts4*scale_factor, histtype='step')
    ax4.set_xlabel(r"Q$^2$ [GeV$^2$/c$^2$]")
    ax4.set_ylabel(r"Count / 0.1 GeV$^2$/c$^2$") 
    plt.savefig(sample_type+"_events_qsq_truth.png")
    plt.close(fig4)  

    # PLOT: truth-level neutrino energy of interaction
    fig5, ax5 = plt.subplots(figsize=(6,4))
    data5 = np.array([d[key]['nu_energy'] for key in d.keys()]) / 1000.
    counts5, bins5 =np.histogram(data5, bins=np.linspace(0,15,31))
    ax5.hist(bins5[:-1], bins=bins5, weights = counts5*scale_factor, histtype='step')
    ax5.set_xlabel(r"Incident Neutrino/Antineutrino Energy [GeV]")
    ax5.set_ylabel("Count / 0.5 GeV") 
    plt.savefig(sample_type+"_events_nu_energy_truth.png")
    plt.close(fig5)      

    # PLOT: truth-level neutrino energy of interaction STACKED HIST BY END_PT_LOC
    loc_labels = [geo_defs.particle_end_loc_dict[k] for k in geo_defs.particle_end_loc_dict.keys()]
    data6f = []; data6u = []; data6d = []; data6b = []; data6s = []; data6p = []
    data7f = []; data7u = []; data7d = []; data7b = []; data7s = []; data7p = []
    data8f = []; data8u = []; data8d = []; data8b = []; data8s = []; data8p = []
    data9f = []; data9u = []; data9d = []; data9b = []; data9s = []; data9p = []
    for key in d.keys():
        if d[key]['end_pt_loc'] == 'f':
            data6f.append(d[key]['nu_energy'] / 1000.)
            data7f.append(d[key]['q2'] / 1000000.)
            data8f.append(d[key]['ang'])
            data9f.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 'd':
            data6d.append(d[key]['nu_energy'] / 1000.)
            data7d.append(d[key]['q2'] / 1000000.)
            data8d.append(d[key]['ang'])
            data9d.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 'b':
            data6b.append(d[key]['nu_energy'] / 1000.)
            data7b.append(d[key]['q2'] / 1000000.)
            data8b.append(d[key]['ang'])
            data9b.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 's':
            data6s.append(d[key]['nu_energy'] / 1000.)
            data7s.append(d[key]['q2'] / 1000000.)   
            data8s.append(d[key]['ang'])
            data9s.append(d[key]['mom'] / 1000.)          
        elif d[key]['end_pt_loc'] == 'p':
            data6p.append(d[key]['nu_energy'] / 1000.)
            data7p.append(d[key]['q2'] / 1000000.)
            data8p.append(d[key]['ang'])
            data9p.append(d[key]['mom'] / 1000.)
        elif d[key]['end_pt_loc'] == 'u':
            data6u.append(d[key]['nu_energy'] / 1000.)
            data7u.append(d[key]['q2'] / 1000000.)
            data8u.append(d[key]['ang'])
            data9u.append(d[key]['mom'] / 1000.)

    print("Minimum momentum of muons punching through MINERvA [GeV/c]:", np.min(data9b))

    fig6, ax6 = plt.subplots(figsize=(9,6))
    bins6 = np.linspace(0,15,31)
    counts6f, bins6f =np.histogram(np.array(data6f), bins=bins6)
    counts6d, bins6d =np.histogram(np.array(data6d), bins=bins6)
    counts6b, bins6b =np.histogram(np.array(data6b), bins=bins6)
    counts6s, bins6s =np.histogram(np.array(data6s), bins=bins6)
    counts6p, bins6p =np.histogram(np.array(data6p), bins=bins6)
    counts6u, bins6u =np.histogram(np.array(data6u), bins=bins6)
    ax6.hist((bins6f[:-1],bins6d[:-1],bins6b[:-1],bins6s[:-1],bins6p[:-1],bins6u[:-1]), bins=bins6, \
        weights = (counts6f*scale_factor,counts6d*scale_factor,counts6b*scale_factor,counts6s*scale_factor,counts6p*scale_factor,counts6u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax6.set_xlabel(r"Incident Neutrino/Antineutrino Energy [GeV]")
    ax6.set_title(sample_title+' Event Neutrino/Antineutrino Energy Spectrum by Muon Track End Behavior')
    ax6.set_ylabel("Count / 0.5 GeV") 
    ax6.legend(loc='upper right')
    plt.savefig(sample_type+"_events_nu_energy_truth_stacked_by_muon_end_loc.png")
    plt.close(fig6)   

    fig7, ax7 = plt.subplots(figsize=(9,6))
    bins7 = np.linspace(0,5,51)
    counts7f, bins7f =np.histogram(np.array(data7f), bins=bins7)
    counts7d, bins7d =np.histogram(np.array(data7d), bins=bins7)
    counts7b, bins7b =np.histogram(np.array(data7b), bins=bins7)
    counts7s, bins7s =np.histogram(np.array(data7s), bins=bins7)
    counts7p, bins7p =np.histogram(np.array(data7p), bins=bins7)
    counts7u, bins7u =np.histogram(np.array(data7u), bins=bins7)
    ax7.hist((bins7f[:-1],bins7d[:-1],bins7b[:-1],bins7s[:-1],bins7p[:-1],bins7u[:-1]), bins=bins7, \
        weights = (counts7f*scale_factor,counts7d*scale_factor,counts7b*scale_factor,counts7s*scale_factor,counts7p*scale_factor,counts7u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax7.legend(loc='upper right')
    ax7.set_ylabel(r"Count / 0.1 GeV$^2$/c$^2$") 
    ax7.set_xlabel(r"Q$^2$ [GeV$^2$/c$^2$]")
    ax7.set_title(sample_title+r' Event Q$^2$ by Muon Track End Behavior')
    plt.savefig(sample_type+"_events_qsq_truth_stacked_by_muon_end_loc.png") 
    plt.close(fig7) 

    fig8, ax8 = plt.subplots(figsize=(9,6))
    bins8 = bins=np.linspace(0.00,0.80,41)
    counts8f, bins8f =np.histogram(np.array(data8f), bins=bins8)
    counts8d, bins8d =np.histogram(np.array(data8d), bins=bins8)
    counts8b, bins8b =np.histogram(np.array(data8b), bins=bins8)
    counts8s, bins8s =np.histogram(np.array(data8s), bins=bins8)
    counts8p, bins8p =np.histogram(np.array(data8p), bins=bins8)
    counts8u, bins8u =np.histogram(np.array(data8u), bins=bins8)
    ax8.hist((bins8f[:-1],bins8d[:-1],bins8b[:-1],bins8s[:-1],bins8p[:-1],bins8u[:-1]), bins=bins8, \
        weights = (counts8f*scale_factor,counts8d*scale_factor,counts8b*scale_factor,counts8s*scale_factor,counts8p*scale_factor,counts8u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax8.legend(loc='upper right')
    ax8.set_title(sample_title+' Event Outgoing Muon Angle by Muon Track End Behavior')
    ax8.set_xlabel(r"Outgoing Muon Angle with Beam Direction")
    ax8.set_ylabel("Count / 0.02 Rad")
    plt.savefig(sample_type+"_events_muon_angle_truth_stacked_by_muon_end_loc.png") 
    plt.close(fig8) 

    fig9, ax9 = plt.subplots(figsize=(9,6))
    bins9 = np.linspace(0,15,31)
    counts9f, bins9f =np.histogram(np.array(data9f), bins=bins9)
    counts9d, bins9d =np.histogram(np.array(data9d), bins=bins9)
    counts9b, bins9b =np.histogram(np.array(data9b), bins=bins9)
    counts9s, bins9s =np.histogram(np.array(data9s), bins=bins9)
    counts9p, bins9p =np.histogram(np.array(data9p), bins=bins9)
    counts9u, bins9u =np.histogram(np.array(data9u), bins=bins9)
    ax9.hist((bins9f[:-1],bins9d[:-1],bins9b[:-1],bins9s[:-1],bins9p[:-1],bins9u[:-1]), bins=bins9, \
        weights = (counts9f*scale_factor,counts9d*scale_factor,counts9b*scale_factor,counts9s*scale_factor,counts9p*scale_factor,counts9u*scale_factor), histtype='bar', label=loc_labels, stacked='True')
    ax9.legend(loc='upper right')
    ax9.set_title(sample_title+' Event Outgoing Muon Momentum by Muon Track End Behavior')
    ax9.set_xlabel(r"Outgoing Muon Momentum [GeV/c]")
    ax9.set_ylabel("Count / 0.5 GeV/c")
    plt.savefig(sample_type+"_events_outgoing_muon_momentum_truth_stacked_by_muon_end_loc.png")
    plt.close(fig9) 

     # PLOT: truth-level neutrino energy of interaction STACKED HIST BY END_PT_LOC
    loc_labels = ['Other', 'COH', 'DIS', 'RES', 'MEC', 'QES']
    data10qes = []; data10mec = []; data10res = []; data10dis = []; data10coh = []; data10und = []
    for key in d.keys():
        if d[key]['nu_int_type'] == 'QES':
            data10qes.append(d[key]['nu_energy'] / 1000.)
        elif d[key]['nu_int_type'] == 'MEC':
            data10mec.append(d[key]['nu_energy'] / 1000.)
        elif d[key]['nu_int_type'] == 'RES':
            data10res.append(d[key]['nu_energy'] / 1000.)
        elif d[key]['nu_int_type'] == 'DIS':
            data10dis.append(d[key]['nu_energy'] / 1000.)        
        elif d[key]['nu_int_type'] == 'COH':
            data10coh.append(d[key]['nu_energy'] / 1000.)
        elif d[key]['nu_int_type'] == 'UND':
            data10und.append(d[key]['nu_energy'] / 1000.)

    # PLOT: truth-level neutrino energy of interaction STACKED HIST BY Neutrino Interaction Mechanism
    fig10, ax10 = plt.subplots(figsize=(9,6))
    bins10 = np.linspace(0,15,31)
    counts10und, bins10und =np.histogram(np.array(data10und), bins=bins10)
    counts10coh, bins10coh =np.histogram(np.array(data10coh), bins=bins10)
    counts10dis, bins10dis =np.histogram(np.array(data10dis), bins=bins10)
    counts10res, bins10res =np.histogram(np.array(data10res), bins=bins10)
    counts10mec, bins10mec =np.histogram(np.array(data10mec), bins=bins10)
    counts10qes, bins10qes =np.histogram(np.array(data10qes), bins=bins10)
    ax10.hist((bins10und[:-1],bins10coh[:-1],bins10dis[:-1],bins10res[:-1],bins10mec[:-1],bins10qes[:-1]), bins=bins10, \
        weights = (counts10und*scale_factor,counts10coh*scale_factor,counts10dis*scale_factor,counts10res*scale_factor,counts10mec*scale_factor,counts10qes*scale_factor), \
        histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
    ax10.set_xlabel(r"Incident Neutrino Energy [GeV]")
    ax10.set_title(sample_title+' Event Neutrino Energy Spectrum by Neutrino Interaction Mechanism')
    ax10.set_ylabel("Count / 0.5 GeV") 
    ax10.set_xlim(0,15)
    ax10.legend(loc='upper right')
    plt.savefig(sample_type+"_events_nu_energy_truth_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig10)

    # Get Detector Boundaries for Plotting
    # TPCs
    tpc_bounds_x = geo_defs.tpc_bounds(0)
    tpc_bounds_y = geo_defs.tpc_bounds(1)[0] # only one set of dims
    tpc_bounds_z = geo_defs.tpc_bounds(2)
    #print("TPC Bounds X:", tpc_bounds_x)
    #print("TPC Bounds Y:", tpc_bounds_y)
    #print("TPC Bounds Z:", tpc_bounds_z)
    # MINERvA
    MINERvA_bounds_x = geo_defs.MINERvA_bounds(0)[0] # only one set of dims
    MINERvA_bounds_y = geo_defs.MINERvA_bounds(1)[0] # only one set of dims
    MINERvA_bounds_z = geo_defs.MINERvA_bounds(2)
    #print("MINERvA Bounds X:", MINERvA_bounds_x)
    #print("MINERvA Bounds Y:", MINERvA_bounds_y)
    #print("MINERvA Bounds Z:", MINERvA_bounds_z)


    # PLOT: truth-level muon start location
    fig11, ax11 = plt.subplots(figsize=(8,6))
    data11x = np.array([d[key]['muon_start'][0] for key in d.keys()])
    data11y = np.array([d[key]['muon_start'][1] for key in d.keys()])
    bins11xa = np.linspace(-500,500,101)
    bins11ya = np.linspace(-500-268,500-268,101)
    counts11, bins11x, bins11y = np.histogram2d(np.array(data11x), np.array(data11y), bins=(bins11xa, bins11ya))
    b11xmesh, b11ymesh = np.meshgrid(bins11x, bins11y)
    counts11 = counts11.T # NOTE: Hist2D doesn't follow Cartesian coords, so need to transpose counts for plotting
    plot_start = ax11.pcolormesh(b11xmesh, b11ymesh,counts11*scale_factor)
    ax11.set_xlabel(r"X [cm]")
    ax11.set_title(sample_title+' Event Muon Truth XY Start Position')
    ax11.set_ylabel(r"Y [cm]") 
    cbar = fig11.colorbar(mappable=plot_start, ax=ax11)
    cbar.set_label(r"Events / 100 cm$^2$")
    # TPC Dimensions in XY
    ax11.text(tpc_bounds_x[0][0], tpc_bounds_y[1]+10, r'2x2', color='red',weight='bold')
    ax11.plot(np.array(tpc_bounds_x[0]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[1]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[2]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[3]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[0]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[1]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[2]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax11.plot(np.array(tpc_bounds_x[3]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[0][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[0][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[1][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[1][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[2][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[2][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[3][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax11.vlines(tpc_bounds_x[3][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    # MINERvA modules
    ax11.text(MINERvA_bounds_x[0], MINERvA_bounds_y[1]+10, r'MINER$\nu$A', color='magenta',weight='bold')
    ax11.plot(np.array(MINERvA_bounds_x), np.full(2, MINERvA_bounds_y[0]), color='magenta', linestyle="-", linewidth=1.)
    ax11.plot(np.array(MINERvA_bounds_x), np.full(2, MINERvA_bounds_y[1]), color='magenta', linestyle="-", linewidth=1.)
    ax11.vlines(MINERvA_bounds_x[0],MINERvA_bounds_y[0],MINERvA_bounds_y[1], color='magenta', linestyle="-", linewidth=1.)
    ax11.vlines(MINERvA_bounds_x[1],MINERvA_bounds_y[0],MINERvA_bounds_y[1], color='magenta', linestyle="-", linewidth=1.)
    plt.savefig(sample_type+"_events_muon_start_xy_truth_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig11)

    # PLOT: truth-level muon start location
    fig12, ax12 = plt.subplots(figsize=(8,6))
    data12x = np.array([d[key]['muon_end'][0] for key in d.keys()])
    data12y = np.array([d[key]['muon_end'][1] for key in d.keys()])
    bins12xa = np.linspace(-500,500,101)
    bins12ya = np.linspace(-500-268,500-268,101)
    counts12, bins12x, bins12y = np.histogram2d(np.array(data12x), np.array(data12y), bins=(bins12xa, bins12ya))
    b12xmesh, b12ymesh = np.meshgrid(bins12x, bins12y)
    counts12 = counts12.T # NOTE: Hist2D doesn't follow Cartesian coords, so need to transpose counts for plotting
    plot_end = ax12.pcolormesh(b12xmesh, b12ymesh,counts12*scale_factor)
    ax12.set_xlabel(r"X [cm]")
    ax12.set_title(sample_title+' Event Muon Truth XY End Position')
    ax12.set_ylabel(r"Y [cm]") 
    cbar = fig12.colorbar(mappable=plot_end, ax=ax12)
    cbar.set_label(r"Events / 100 cm$^2$")
    # TPC Dimensions in XY
    ax12.text(tpc_bounds_x[0][0], tpc_bounds_y[1]+10, r'2x2', color='red',weight='bold')
    ax12.plot(np.array(tpc_bounds_x[0]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[1]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[2]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[3]), np.full(2, tpc_bounds_y[0]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[0]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[1]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[2]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax12.plot(np.array(tpc_bounds_x[3]), np.full(2, tpc_bounds_y[1]), color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[0][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[0][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[1][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[1][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[2][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[2][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[3][0],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    ax12.vlines(tpc_bounds_x[3][1],tpc_bounds_y[0],tpc_bounds_y[1], color='red', linestyle="-", linewidth=1)
    # MINERvA modules
    ax12.text(MINERvA_bounds_x[0], MINERvA_bounds_y[1]+10, r'MINER$\nu$A', color='magenta',weight='bold')
    ax12.plot(np.array(MINERvA_bounds_x), np.full(2, MINERvA_bounds_y[0]), color='magenta', linestyle="-", linewidth=1.)
    ax12.plot(np.array(MINERvA_bounds_x), np.full(2, MINERvA_bounds_y[1]), color='magenta', linestyle="-", linewidth=1.)
    ax12.vlines(MINERvA_bounds_x[0],MINERvA_bounds_y[0],MINERvA_bounds_y[1], color='magenta', linestyle="-", linewidth=1.)
    ax12.vlines(MINERvA_bounds_x[1],MINERvA_bounds_y[0],MINERvA_bounds_y[1], color='magenta', linestyle="-", linewidth=1.)
    plt.savefig(sample_type+"_events_muon_end_xy_truth_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig12)

    # Neutrino vs. Antineutrino Events
    fig13, ax13 = plt.subplots(figsize=(6,4))
    nu_pdg_list=[d[key]['nu_pdg'] for key in d.keys()]
    parent_pdg_list=[d[key]['parent_pdg'] for key in d.keys()]
    print("Parent PDG values:", set(parent_pdg_list))
    print("Nu PDG values:", set(nu_pdg_list))
    nu_pdg_set=set(pdg for pdg in nu_pdg_list)
    nu_pdg_count=[(pdg, nu_pdg_list.count(pdg)) for pdg in nu_pdg_set]
    nu_pdg_fraction=[100*(i[1]/len(np.array(nu_pdg_list))) for i in nu_pdg_count]
    if nu_pdg_count[0][0] == 14:
        nu_pdg_labels=['Neutrino']
        if len(nu_pdg_set) > 1:
            nu_pdg_labels.append('Antineutrino')
    else:
        nu_pdg_labels=['Antineutrino']
        if len(nu_pdg_set) > 1:
            nu_pdg_labels.append('Neutrino')
    ax13.pie(nu_pdg_fraction, labels=nu_pdg_labels, autopct='%1.1f%%')
    ax13.set_title(sample_title+r" Event Neutrino vs. Antineutrino Breakdown")
    plt.savefig(sample_type+"_events_neutrino_vs_antineutrino_truth.png")
    plt.close(fig13)    