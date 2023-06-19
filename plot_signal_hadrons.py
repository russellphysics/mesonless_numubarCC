import matplotlib
import matplotlib.pyplot as plt
import h5py
import glob
import json
import argparse
import numpy as np
import twoBytwo_defs
import auxiliary
import signal_characterization as sig_char

# PLOT: Hadron kinematics
#       sig_bkg is an int such that 0 == signal, 1 == 'dirt' backgrounds, 2 == 'beam' backgrounds
def plot_hadrons(d, scale_factor, sig_bkg = 0):
    
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
        return "Error: plot_hadrons function given undefined signal/background definition"
    
    # PLOT: total visible energy + contained visible energy
    fig0, ax0 = plt.subplots(figsize=(8,4))
    data0tot = np.array([d[key]['total_edep'] for key in d.keys()])
    #data0cont = np.array([d[key]['contained_edep'] for key in d.keys()])
    counts0tot, bins0tot = np.histogram(data0tot, bins=np.linspace(0,800,40))
    #counts0cont, bins0cont = np.histogram(data0cont, bins=np.linspace(0,800,40))
    ax0.hist(bins0tot[:-1], bins=bins0tot, weights = counts0tot*scale_factor, label='Total', histtype='step')
    #ax0.hist(bins0cont[:-1], bins=bins0cont, weights = counts0cont*scale_factor, label='Contained',histtype='step', linestyle='--')
    ax0.set_xlabel('Total Visible Hadron Energy [MeV]')
    ax0.set_ylabel('Count / 20 MeV')
    ax0.set_yscale('log')
    #ax0.legend()
    ax0.grid(True)
    plt.savefig(sample_type+"_events_hadron_visible_energy.png")
    plt.close(fig0)

    # PLOT: hadron energy containment fraction ## UNNECESSARY PLOT AT THE MOMENT
    fig1, ax1 = plt.subplots(figsize=(6,4))
    data1 = np.array([d[key]['contained_edep']/d[key]['total_edep'] for key in d.keys() if d[key]['total_edep']!=0])
    counts1, bins1 = np.histogram(data1, bins=np.linspace(0,1,20))
    ax1.hist(bins1[:-1], bins=bins1, weights = counts1*scale_factor, histtype='step')
    ax1.set_xlabel('Visible Hadron Energy Containment Fraction')
    ax1.set_ylabel('Count / 0.05')
    ax1.grid(True)       
    ax0.set_yscale('log')
    plt.savefig(sample_type+"_events_hadron_energy_containment_fraction.png")
    plt.close(fig1)

    # PLOT: hadron multiplicity
    fig2, ax2 = plt.subplots(figsize=(6,4))
    data2 = np.array([d[key]['hadron_mult'] for key in d.keys()])
    counts2, bins2 =np.histogram(data2, bins=np.linspace(0,25,26))
    ax2.hist(bins2[:-1], bins=bins2, weights = counts2*scale_factor, histtype='step')
    ax2.set_xlabel(r"Primary Hadron Multiplicity")
    ax2.set_ylabel("Count / Hadron") 
    plt.savefig(sample_type+"_events_hadron_multiplicity_truth.png")
    plt.close(fig2)    

    # PLOT: truth-level 4-momentum squared of interaction
    #       ** no scale factor applied because we're looking at fractions anyways ** 
    fig3, ax3 = plt.subplots(figsize=(6,4))
    hadron_fs_pdg_list=[sorted(d[key]['hadron_pdg_set']) for key in d.keys()]
    hadron_fs_pdg_set=set(tuple(pdg) for pdg in hadron_fs_pdg_list)
    #print("Hadron PDG List:", hadron_fs_pdg_list)
    #print("Hadron PDG Set:", hadron_fs_pdg_set)
    hadron_fs_pdg_count=[(pdg_set, hadron_fs_pdg_list.count(list(pdg_set))) for pdg_set in hadron_fs_pdg_set]
    hadron_fs_pdg_fraction=[100*(i[1]/len(data2)) for i in hadron_fs_pdg_count]
    hadron_fs_pdg_labels=['+'.join(str(auxiliary.hadron_pdg_dict[j]) for j in i[0]) for i in hadron_fs_pdg_count]
    #print("Number of Events:", len(hadron_fs_pdg_list))
    #print("Hadron FS PDG Count:", hadron_fs_pdg_count)
    #print("Hadron FS PDG Fractions:", hadron_fs_pdg_fraction)
    #print("Hadron FS PDG Labels:", hadron_fs_pdg_labels)
    ax3.pie(hadron_fs_pdg_fraction, labels=hadron_fs_pdg_labels, autopct='%1.1f%%')
    ax3.set_title(r"Final State Hadrons in "+sample_title+" Events")
    plt.savefig(sample_type+"_events_hadron_pdg_ids_truth.png")
    plt.close(fig3)    

    # PLOT: other hadron multiplicity
    fig4, ax4 = plt.subplots(figsize=(6,4))
    data4 = np.array([d[key]['other_had_mult'] for key in d.keys()])
    counts4, bins4 =np.histogram(data4, bins=np.linspace(0,10,11))
    ax4.hist(bins4[:-1], bins=bins4, weights = counts4*scale_factor, histtype='step')
    ax4.set_xlabel(r"Other Primary Hadron Multiplicity")
    ax4.set_ylabel("Count / Other Hadron") 
    plt.savefig(sample_type+"_events_other_hadron_multiplicity_truth.png")
    plt.close(fig4)    

    # PLOT: neutron multiplicity
    fig5, ax5 = plt.subplots(figsize=(6,4))
    data5 = np.array([d[key]['neutron_mult'] for key in d.keys()])
    counts5, bins5 =np.histogram(data5, bins=np.linspace(0,20,21))
    ax5.hist(bins5[:-1], bins=bins5, weights = counts5*scale_factor, histtype='step')
    ax5.set_xlabel(r"Primary Neutron Multiplicity")
    ax5.set_ylabel("Count / Neutron") 
    plt.savefig(sample_type+"_events_neutron_multiplicity_truth.png")
    plt.close(fig5)    

    # PLOT: proton multiplicity
    fig6, ax6 = plt.subplots(figsize=(6,4))
    data6 = np.array([d[key]['proton_mult'] for key in d.keys()])
    counts6, bins6 =np.histogram(data6, bins=np.linspace(0,20,21))
    ax6.hist(bins6[:-1], bins=bins6, weights = counts6*scale_factor, histtype='step')
    ax6.set_xlabel(r"Primary Proton Multiplicity")
    ax6.set_ylabel("Count / Proton") 
    plt.savefig(sample_type+"_events_proton_multiplicity_truth.png")
    plt.close(fig6)   

    
    # PLOT: Fractions of Events with diff numbers of protons
    #       ** no scale factor applied because we're looking at fractions anyways ** 
    fig7, ax7 = plt.subplots(figsize=(6,4))
    p_mult_list = []
    total_p_events = 0
    for key in d.keys():
        #if d[key]['other_had_mult']==0 and d[key]['proton_mult']>0 and d[key]['neutron_mult']>0:
        if d[key]['proton_mult']>0:
            p_mult_list.append(d[key]['proton_mult'])
            total_p_events+=1
    if total_p_events >0:
        p_mult_count=[(mult_count, p_mult_list.count(mult_count)) for mult_count in np.arange(np.max(np.array(p_mult_list)))]
        #print("P mult count:", p_mult_count)
        p_mult_fraction=[100*(i[1]/total_p_events) for i in p_mult_count if i[1]>0 and i[0]<6]
        p_mult_labels=[str(i[0]) for i in p_mult_count if i[1]>0 and i[0]<6]
        many_p = 0
        for num, count in p_mult_count:
            if num < 6: continue
            if count==0: continue
            many_p+=count
        p_mult_fraction.append(100*(many_p/total_p_events))
        p_mult_labels.append('>5')
        #print("P mult labels:", p_mult_labels)
        ax7.pie(p_mult_fraction, labels=p_mult_labels, autopct='%1.1f%%')
        ax7.set_title(r"Primary Proton Multiplicity in "+sample_title+"\nEvents with Protons")
        plt.savefig(sample_type+"_events_proton_mult_in_p_events_truth.png")
    plt.close(fig7)    

    # PLOT: Max proton length for n p events
    fig8, ax8 = plt.subplots(figsize=(8,4))
    #p_tot_lens = []
    p_cont_lens = []
    events_w_protons = 0
    for key in d.keys():
        if d[key]['proton_mult']>0:
            #p_tot_lens.append(d[key]['max_p_total_length'])
            p_cont_lens.append(d[key]['max_p_contained_length'])
            events_w_protons+=1
    if events_w_protons >0:  
        #data8tot = np.array(p_tot_lens)
        data8cont = np.array(p_cont_lens)
        #counts8tot, bins8tot = np.histogram(data8tot, bins=np.linspace(0,60,60))
        counts8cont, bins8cont = np.histogram(data8cont, bins=np.linspace(0,60,60))
        #ax8.hist(bins8tot[:-1], bins=bins8tot, weights = counts8tot*scale_factor, label='Total', histtype='step')
        ax8.hist(bins8cont[:-1], bins=bins8cont, weights = counts8cont*scale_factor, label='Contained',histtype='step', linestyle='--')
        ax8.set_xlabel(r"Length [cm]")
        ax8.set_title("Contained Length of Longest Proton Track in "+sample_title+" Events with Protons")
        ax8.set_ylabel("Events / cm") 
        #ax8.legend()
        plt.savefig(sample_type+"_events_max_proton_length_in_p_events_truth.png")   
    plt.close(fig8)

    # PLOT: hadron multiplicity above threshold
    fig9, ax9 = plt.subplots(figsize=(6,4))
    data9 = np.array([d[key]['hadron_mult_over_thresh'] for key in d.keys()])
    counts9, bins9 =np.histogram(data9, bins=np.linspace(0,12,13))
    ax9.hist(bins9[:-1], bins=bins9, weights = counts9*scale_factor, histtype='step')
    ax9.set_xlabel(r"Primary Hadron Multiplicity Above Threshold")
    ax9.set_ylabel("Events / Hadron") 
    plt.savefig(sample_type+"_events_hadron_multiplicity_above_threshold.png")
    plt.close(fig9)  

    # PLOT: proton multiplicity above threshold
    fig10, ax10 = plt.subplots(figsize=(6,4))
    data10 = np.array([d[key]['proton_mult_over_thresh'] for key in d.keys()])
    counts10, bins10 =np.histogram(data10, bins=np.linspace(0,12,13))
    ax10.hist(bins10[:-1], bins=bins10, weights = counts10*scale_factor, histtype='step')
    ax10.set_xlabel(r"Primary Proton Multiplicity Above Threshold")
    ax10.set_ylabel("Events / Proton") 
    plt.savefig(sample_type+"_events_proton_multiplicity_above_threshold.png")
    plt.close(fig10)   

    # PLOT: Lead proton momentum for events with protons
    fig11, ax11 = plt.subplots(figsize=(8,4))
    p_lead_mom = []
    events_w_protons = 0
    for key in d.keys():
        if d[key]['proton_mult_over_thresh']>0:
            p_lead_mom.append(d[key]['lead_proton_momentum'])
            events_w_protons+=1
    if events_w_protons >0:  
        data11tot = np.array(p_lead_mom)
        counts11tot, bins11tot = np.histogram(data11tot, bins=np.linspace(0,2000,41))
        ax11.hist(bins11tot[:-1], bins=bins11tot, weights = counts11tot*scale_factor, histtype='step')
        ax11.set_xlabel(r"Momentum [MeV/c]")
        ax11.set_title("Leading Proton Momentum in\n"+sample_title+" Events with Protons Above Threshold")
        ax11.set_ylabel("Events / 50 MeV/c") 
        plt.savefig(sample_type+"_lead_proton_momentum_events_with_protons_above_threshold.png")   
    plt.close(fig11)

    # PLOT: Sub-leading proton momentum for events with 2+ protons
    fig12, ax12 = plt.subplots(figsize=(8,4))
    p_sublead_mom = []
    events_w_greq_2_protons = 0
    for key in d.keys():
        if d[key]['proton_mult_over_thresh']>1:
            p_sublead_mom.append(d[key]['sub_lead_proton_momentum'])
            events_w_greq_2_protons+=1
    if events_w_greq_2_protons >0:  
        data12tot = np.array(p_sublead_mom)
        counts12tot, bins12tot = np.histogram(data12tot, bins=np.linspace(0,2000,41))
        ax12.hist(bins12tot[:-1], bins=bins12tot, weights = counts12tot*scale_factor, histtype='step')
        ax12.set_xlabel(r"Momentum [MeV/c]")
        ax12.set_title("Subleading Proton Momentum in \n"+sample_title\
                       +" Events with 2+ Protons Above Threshold")
        ax12.set_ylabel("Events / 50 MeV/c") 
        plt.savefig(sample_type+"_sublead_proton_momentum_events_with_greq_2_protons_above_threshold.png")   
    plt.close(fig12)

    # PLOT: Lead proton angle wrt beam for events with protons
    fig13, ax13 = plt.subplots(figsize=(8,4))
    p_lead_ang_wrt_beam = []
    events_w_protons = 0
    for key in d.keys():
        if d[key]['proton_mult_over_thresh']>0:
            p_lead_ang_wrt_beam.append(d[key]['lead_proton_ang_wrt_beam'])
            events_w_protons+=1
    if events_w_protons >0:  
        data13tot = np.array(p_lead_ang_wrt_beam)
        counts13tot, bins13tot = np.histogram(data13tot, bins=np.linspace(0,3.2,33))
        ax13.hist(bins13tot[:-1], bins=bins13tot, weights = counts13tot*scale_factor, histtype='step')
        ax13.set_xlabel(r"Angle [Rad]")
        ax13.set_title("Leading Proton Angle with respect to Beam Direction in\n"+sample_title\
                       +" Events with Protons Above Threshold")
        ax13.set_ylabel("Events / 0.1 Rad") 
        plt.savefig(sample_type+"_lead_proton_angle_wrt_beam_direction_events_with_protons_above_threshold.png")   
    plt.close(fig13)

    # PLOT: Subleading proton angle wrt beam for events with protons
    fig14, ax14 = plt.subplots(figsize=(8,4))
    p_sublead_ang_wrt_beam = []
    events_w_greq_2_protons = 0
    for key in d.keys():
        if d[key]['proton_mult_over_thresh']>1:
            p_sublead_ang_wrt_beam.append(d[key]['sub_lead_proton_ang_wrt_beam'])
            events_w_greq_2_protons+=1
    if events_w_greq_2_protons >0:  
        data14tot = np.array(p_sublead_ang_wrt_beam)
        counts14tot, bins14tot = np.histogram(data14tot, bins=np.linspace(0,3.2,33))
        ax14.hist(bins14tot[:-1], bins=bins14tot, weights = counts14tot*scale_factor, histtype='step')
        ax14.set_xlabel(r"Angle [Rad]")
        ax14.set_title("Subleading Proton Angle with respect to Beam Direction in\n"+sample_title\
                       +" Events with 2+ Protons Above Threshold")
        ax14.set_ylabel("Events / 0.1 Rad") 
        plt.savefig(sample_type+"_sublead_proton_angle_wrt_beam_direction_events_with_greq_2_protons_above_threshold.png")   
    plt.close(fig14)

    # PLOT: Subleading proton angle wrt leading proton for events with protons
    fig15, ax15 = plt.subplots(figsize=(8,4))
    p_sublead_ang_wrt_lead_proton = []
    events_w_greq_2_protons = 0
    for key in d.keys():
        if d[key]['proton_mult_over_thresh']>1:
            p_sublead_ang_wrt_lead_proton.append(d[key]['sub_lead_proton_angle_with_lead_proton'])
            events_w_greq_2_protons+=1
    if events_w_greq_2_protons >0:  
        data15tot = np.array(p_sublead_ang_wrt_lead_proton)
        counts15tot, bins15tot = np.histogram(data15tot, bins=np.linspace(0,1.6,17))
        ax15.hist(bins15tot[:-1], bins=bins15tot, weights = counts15tot*scale_factor, histtype='step')
        ax15.set_xlabel(r"Angle [Rad]")
        ax15.set_title("Subleading Proton Angle with respect to Leading Proton Direction\n in "+sample_title \
                       +" Events with 2+ Protons Above Threshold")
        ax15.set_ylabel("Events / 0.1 Rad") 
        plt.savefig(sample_type+"_sublead_proton_angle_wrt_lead_proton_events_with_greq_2_protons_above_threshold.png")   
    plt.close(fig15)

    # PLOT: proton multiplicity truth vs. over threshold
    fig16, ax16 = plt.subplots(figsize=(6,4))
    data16_tr = np.array([d[key]['proton_mult'] for key in d.keys()])
    data16_thresh= np.array([d[key]['proton_mult_over_thresh'] for key in d.keys()])
    counts16_tr, bins16_tr =np.histogram(data16_tr, bins=np.linspace(0,20,21))
    counts16_thresh, bins16_thresh =np.histogram(data16_thresh, bins=np.linspace(0,20,21))
    ax16.hist(bins16_tr[:-1], bins=bins16_tr, weights = counts16_tr*scale_factor, label="Truth", histtype='step')
    ax16.hist(bins16_thresh[:-1], bins=bins16_thresh, weights = counts16_thresh*scale_factor, label="Over Threshold", histtype='step', linestyle='--')
    ax16.set_xlabel(r"Primary Proton Multiplicity")
    ax16.set_ylabel("Events / Proton") 
    plt.legend()
    plt.savefig(sample_type+"_events_proton_multiplicity_truth_vs_over_threshold.png")
    plt.close(fig16)   

    # PLOT: Primary Proton K.E. for events with Protons Over Threshold
    fig17, ax17 = plt.subplots(figsize=(8,4))
    p_ke = []
    events_w_protons = 0
    for key in d.keys():
        if d[key]['proton_mult_over_thresh']>0:
            p_ke.append(d[key]['primary_protons_total_ke'])
            events_w_protons+=1
    if events_w_protons >0:  
        data17tot = np.array(p_ke)
        counts17tot, bins17tot = np.histogram(data17tot, bins=np.linspace(0,2000,41))
        ax17.hist(bins17tot[:-1], bins=bins17tot, weights = counts17tot*scale_factor, histtype='step')
        ax17.set_xlabel(r"Kinetic Energy [MeV]")
        ax17.set_title("Total Primary Proton Kinetic Energy in\n"+sample_title+" Events with Protons Above Threshold")
        ax17.set_ylabel("Events / 50 MeV") 
        plt.savefig(sample_type+"_primary_proton_ke_events_with_protons_above_threshold.png")   
    plt.close(fig17)

    # PLOT: Truth KE vs. Contained KE for primary protons
    fig18, ax18 = plt.subplots(figsize=(8,4))
    data18tot = np.array([d[key]['primary_protons_total_ke'] for key in d.keys()])
    data18cont = np.array([d[key]['contained_edep_over_thresh'] for key in d.keys()])
    counts18tot, bins18tot = np.histogram(data18tot, bins=np.linspace(0,2000,101))
    counts18cont, bins18cont = np.histogram(data18cont, bins=np.linspace(0,2000,101))
    ax18.hist(bins18tot[:-1], bins=bins18tot, weights = counts18tot*scale_factor, label='Total', histtype='step')
    ax18.hist(bins18cont[:-1], bins=bins18cont, weights = counts18cont*scale_factor, label='Contained',histtype='step', linestyle='--')
    ax18.set_xlabel('Primary Proton Energy [MeV]')
    ax18.set_ylabel('Events / 20 MeV')
    ax18.set_yscale('log')
    ax18.legend()
    plt.savefig(sample_type+"_events_primary_proton_truth_ke_vs_contained_ke.png")
    plt.close(fig18)

     # PLOT: truth-level hadron information STACKED HISTs BY neutrino interaction mechanism
    loc_labels = ['Other', 'COH', 'DIS', 'RES', 'MEC', 'QES']
    data19qes = []; data19mec = []; data19res = []; data19dis = []; data19coh = []; data19und = []
    data20qes = []; data20mec = []; data20res = []; data20dis = []; data20coh = []; data20und = []
    data21qes = []; data21mec = []; data21res = []; data21dis = []; data21coh = []; data21und = []
    data22qes = []; data22mec = []; data22res = []; data22dis = []; data22coh = []; data22und = []
    data23qes = []; data23mec = []; data23res = []; data23dis = []; data23coh = []; data23und = []
    data24qes = []; data24mec = []; data24res = []; data24dis = []; data24coh = []; data24und = []
    data25qes = []; data25mec = []; data25res = []; data25dis = []; data25coh = []; data25und = []
    data26qes = []; data26mec = []; data26res = []; data26dis = []; data26coh = []; data26und = []
    data27qes = []; data27mec = []; data27res = []; data27dis = []; data27coh = []; data27und = []
    data28qes = []; data28mec = []; data28res = []; data28dis = []; data28coh = []; data28und = []
    data29qes = []; data29mec = []; data29res = []; data29dis = []; data29coh = []; data29und = []
    events_with_protons = 0; events_with_protons_over_thresh = 0; events_with_geq_2_protons_over_thresh =0
    for key in d.keys():
        if d[key]['nu_int_type'] == 'QES':
            data19qes.append(d[key]['proton_mult'])
            data20qes.append(d[key]['proton_mult_over_thresh'])
            data21qes.append(d[key]['hadron_mult'])
            data22qes.append(d[key]['hadron_mult_over_thresh'])
            if d[key]['proton_mult']>0:
                data23qes.append(d[key]['max_p_contained_length'])
                events_with_protons +=1
            if d[key]['proton_mult_over_thresh'] > 0:
                data24qes.append(d[key]['primary_protons_total_ke'])
                data25qes.append(d[key]['lead_proton_momentum'])
                data26qes.append(d[key]['lead_proton_ang_wrt_beam'])
                events_with_protons_over_thresh +=1
                if d[key]['proton_mult_over_thresh'] > 1:
                    data27qes.append(d[key]['sub_lead_proton_momentum'])
                    data28qes.append(d[key]['sub_lead_proton_ang_wrt_beam'])
                    data29qes.append(d[key]['sub_lead_proton_angle_with_lead_proton'])
                    events_with_geq_2_protons_over_thresh +=1
        elif d[key]['nu_int_type'] == 'MEC':
            data19mec.append(d[key]['proton_mult'])
            data20mec.append(d[key]['proton_mult_over_thresh'])
            data21mec.append(d[key]['hadron_mult'])
            data22mec.append(d[key]['hadron_mult_over_thresh'])
            if d[key]['proton_mult']>0:
                data23mec.append(d[key]['max_p_contained_length'])
                events_with_protons +=1
            if d[key]['proton_mult_over_thresh'] > 0:
                data24mec.append(d[key]['primary_protons_total_ke'])
                data25mec.append(d[key]['lead_proton_momentum'])
                data26mec.append(d[key]['lead_proton_ang_wrt_beam'])
                events_with_protons_over_thresh +=1
                if d[key]['proton_mult_over_thresh'] > 1:
                    data27mec.append(d[key]['sub_lead_proton_momentum'])
                    data28mec.append(d[key]['sub_lead_proton_ang_wrt_beam'])
                    data29mec.append(d[key]['sub_lead_proton_angle_with_lead_proton'])
                    events_with_geq_2_protons_over_thresh +=1
        elif d[key]['nu_int_type'] == 'RES':
            data19res.append(d[key]['proton_mult'])
            data20res.append(d[key]['proton_mult_over_thresh'])
            data21res.append(d[key]['hadron_mult'])
            data22res.append(d[key]['hadron_mult_over_thresh'])
            if d[key]['proton_mult']>0:
                data23res.append(d[key]['max_p_contained_length'])
                events_with_protons +=1
            if d[key]['proton_mult_over_thresh'] > 0:
                data24res.append(d[key]['primary_protons_total_ke'])
                data25res.append(d[key]['lead_proton_momentum'])
                data26res.append(d[key]['lead_proton_ang_wrt_beam'])
                events_with_protons_over_thresh +=1
                if d[key]['proton_mult_over_thresh'] > 1:
                    data27res.append(d[key]['sub_lead_proton_momentum'])
                    data28res.append(d[key]['sub_lead_proton_ang_wrt_beam'])
                    data29res.append(d[key]['sub_lead_proton_angle_with_lead_proton'])
                    events_with_geq_2_protons_over_thresh +=1
        elif d[key]['nu_int_type'] == 'DIS':
            data19dis.append(d[key]['proton_mult'])
            data20dis.append(d[key]['proton_mult_over_thresh'])
            data21dis.append(d[key]['hadron_mult'])
            data22dis.append(d[key]['hadron_mult_over_thresh'])
            if d[key]['proton_mult']>0:
                data23dis.append(d[key]['max_p_contained_length'])
                events_with_protons +=1
            if d[key]['proton_mult_over_thresh'] > 0:
                data24dis.append(d[key]['primary_protons_total_ke'])
                data25dis.append(d[key]['lead_proton_momentum'])
                data26dis.append(d[key]['lead_proton_ang_wrt_beam'])
                events_with_protons_over_thresh +=1
                if d[key]['proton_mult_over_thresh'] > 1:
                    data27dis.append(d[key]['sub_lead_proton_momentum'])
                    data28dis.append(d[key]['sub_lead_proton_ang_wrt_beam'])
                    data29dis.append(d[key]['sub_lead_proton_angle_with_lead_proton'])   
                    events_with_geq_2_protons_over_thresh +=1
        elif d[key]['nu_int_type'] == 'COH':
            data19coh.append(d[key]['proton_mult'])
            data20coh.append(d[key]['proton_mult_over_thresh'])
            data21coh.append(d[key]['hadron_mult'])
            data22coh.append(d[key]['hadron_mult_over_thresh'])
            if d[key]['proton_mult']>0:
                data23coh.append(d[key]['max_p_contained_length'])
                events_with_protons +=1
            if d[key]['proton_mult_over_thresh'] > 0:
                data24coh.append(d[key]['primary_protons_total_ke'])
                data25coh.append(d[key]['lead_proton_momentum'])
                data26coh.append(d[key]['lead_proton_ang_wrt_beam'])
                events_with_protons_over_thresh +=1
                if d[key]['proton_mult_over_thresh'] > 1:
                    data27coh.append(d[key]['sub_lead_proton_momentum'])
                    data28coh.append(d[key]['sub_lead_proton_ang_wrt_beam'])
                    data29coh.append(d[key]['sub_lead_proton_angle_with_lead_proton'])
                    events_with_geq_2_protons_over_thresh +=1
        elif d[key]['nu_int_type'] == 'UND':
            data19und.append(d[key]['proton_mult'])
            data20und.append(d[key]['proton_mult_over_thresh'])
            data21und.append(d[key]['hadron_mult'])
            data22und.append(d[key]['hadron_mult_over_thresh'])
            if d[key]['proton_mult']>0:
                data23und.append(d[key]['max_p_contained_length'])
                events_with_protons +=1
            if d[key]['proton_mult_over_thresh'] > 0:
                data24und.append(d[key]['primary_protons_total_ke'])
                data25und.append(d[key]['lead_proton_momentum'])
                data26und.append(d[key]['lead_proton_ang_wrt_beam'])
                events_with_protons_over_thresh +=1
                if d[key]['proton_mult_over_thresh'] > 1:
                    data27und.append(d[key]['sub_lead_proton_momentum'])
                    data28und.append(d[key]['sub_lead_proton_ang_wrt_beam'])
                    data29und.append(d[key]['sub_lead_proton_angle_with_lead_proton'])
                    events_with_geq_2_protons_over_thresh +=1

    # Proton Multiplicity by Neutrino Interaction Mechanism
    fig19, ax19 = plt.subplots(figsize=(6,4))
    bins19 = np.linspace(0,17,18)
    counts19und, bins19und =np.histogram(np.array(data19und), bins=bins19)
    counts19coh, bins19coh =np.histogram(np.array(data19coh), bins=bins19)
    counts19dis, bins19dis =np.histogram(np.array(data19dis), bins=bins19)
    counts19res, bins19res =np.histogram(np.array(data19res), bins=bins19)
    counts19mec, bins19mec =np.histogram(np.array(data19mec), bins=bins19)
    counts19qes, bins19qes =np.histogram(np.array(data19qes), bins=bins19)
    ax19.hist((bins19und[:-1],bins19coh[:-1],bins19dis[:-1],bins19res[:-1],bins19mec[:-1],bins19qes[:-1]), bins=bins19, \
        weights = (counts19und*scale_factor,counts19coh*scale_factor,counts19dis*scale_factor,counts19res*scale_factor,counts19mec*scale_factor,counts19qes*scale_factor), \
        histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
    ax19.set_xlabel(r"Primary Proton Multiplicity")
    ax19.set_title(sample_title+' Event Primary Proton Multiplicity by Neutrino Interaction Mechanism')
    ax19.set_ylabel("Events / Proton") 
    ax19.legend(loc='upper right')
    plt.savefig(sample_type+"_events_proton_multiplicity_truth_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig19)

    # Proton Multiplicity Over Threshold by Neutrino Interaction Mechanism
    fig20, ax20 = plt.subplots(figsize=(6,4))
    bins20 = np.linspace(0,17,18)
    counts20und, bins20und =np.histogram(np.array(data20und), bins=bins20)
    counts20coh, bins20coh =np.histogram(np.array(data20coh), bins=bins20)
    counts20dis, bins20dis =np.histogram(np.array(data20dis), bins=bins20)
    counts20res, bins20res =np.histogram(np.array(data20res), bins=bins20)
    counts20mec, bins20mec =np.histogram(np.array(data20mec), bins=bins20)
    counts20qes, bins20qes =np.histogram(np.array(data20qes), bins=bins20)
    ax20.hist((bins20und[:-1],bins20coh[:-1],bins20dis[:-1],bins20res[:-1],bins20mec[:-1],bins20qes[:-1]), bins=bins20, \
        weights = (counts20und*scale_factor,counts20coh*scale_factor,counts20dis*scale_factor,counts20res*scale_factor,counts20mec*scale_factor,counts20qes*scale_factor), \
        histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
    ax20.set_xlabel(r"Primary Proton Multiplicity Above Threshold")
    ax20.set_title(sample_title+' Event Primary Proton Multiplicity \nAbove Threshold by Neutrino Interaction Mechanism')
    ax20.set_ylabel("Events / Proton") 
    ax20.legend(loc='upper right')
    plt.savefig(sample_type+"_events_proton_multiplicity_over_threshold_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig20)

    # Hadron Multiplicity by Neutrino Interaction Mechanism
    fig21, ax21 = plt.subplots(figsize=(6,4))
    bins21 = np.linspace(0,25,26)
    counts21und, bins21und =np.histogram(np.array(data21und), bins=bins21)
    counts21coh, bins21coh =np.histogram(np.array(data21coh), bins=bins21)
    counts21dis, bins21dis =np.histogram(np.array(data21dis), bins=bins21)
    counts21res, bins21res =np.histogram(np.array(data21res), bins=bins21)
    counts21mec, bins21mec =np.histogram(np.array(data21mec), bins=bins21)
    counts21qes, bins21qes =np.histogram(np.array(data21qes), bins=bins21)
    ax21.hist((bins21und[:-1],bins21coh[:-1],bins21dis[:-1],bins21res[:-1],bins21mec[:-1],bins21qes[:-1]), bins=bins21, \
        weights = (counts21und*scale_factor,counts21coh*scale_factor,counts21dis*scale_factor,counts21res*scale_factor,counts21mec*scale_factor,counts21qes*scale_factor), \
        histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
    ax21.set_xlabel(r"Primary Hadron Multiplicity")
    ax21.set_title(sample_title+' Event Primary Hadron Multiplicity by Neutrino Interaction Mechanism')
    ax21.set_ylabel("Events / Hadron") 
    ax21.legend(loc='upper right')
    plt.savefig(sample_type+"_events_hadron_multiplicity_truth_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig21)

    # Hadron Multiplicity Over Threshold by Neutrino Interaction Mechanism
    fig22, ax22 = plt.subplots(figsize=(6,4))
    bins22 = np.linspace(0,25,26)
    counts22und, bins22und =np.histogram(np.array(data22und), bins=bins22)
    counts22coh, bins22coh =np.histogram(np.array(data22coh), bins=bins22)
    counts22dis, bins22dis =np.histogram(np.array(data22dis), bins=bins22)
    counts22res, bins22res =np.histogram(np.array(data22res), bins=bins22)
    counts22mec, bins22mec =np.histogram(np.array(data22mec), bins=bins22)
    counts22qes, bins22qes =np.histogram(np.array(data22qes), bins=bins22)
    ax22.hist((bins22und[:-1],bins22coh[:-1],bins22dis[:-1],bins22res[:-1],bins22mec[:-1],bins22qes[:-1]), bins=bins22, \
        weights = (counts22und*scale_factor,counts22coh*scale_factor,counts22dis*scale_factor,counts22res*scale_factor,counts22mec*scale_factor,counts22qes*scale_factor), \
        histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
    ax22.set_xlabel(r"Primary Hadron Multiplicity Above Threshold")
    ax22.set_title(sample_title+' Event Primary Hadron Multiplicity \nAbove Threshold by Neutrino Interaction Mechanism')
    ax22.set_ylabel("Events / Hadron") 
    ax22.legend(loc='upper right')
    plt.savefig(sample_type+"_events_hadron_multiplicity_above_threshold_stacked_by_neutrino_interaction_mechanism.png")
    plt.close(fig22)

    if events_with_protons > 0:
        # Max Contained Proton Length Events with Protons
        fig23, ax23 = plt.subplots(figsize=(8,4))
        bins23 = np.linspace(0,60,61)
        counts23und, bins23und =np.histogram(np.array(data23und), bins=bins23)
        counts23coh, bins23coh =np.histogram(np.array(data23coh), bins=bins23)
        counts23dis, bins23dis =np.histogram(np.array(data23dis), bins=bins23)
        counts23res, bins23res =np.histogram(np.array(data23res), bins=bins23)
        counts23mec, bins23mec =np.histogram(np.array(data23mec), bins=bins23)
        counts23qes, bins23qes =np.histogram(np.array(data23qes), bins=bins23)
        ax23.hist((bins23und[:-1],bins23coh[:-1],bins23dis[:-1],bins23res[:-1],bins23mec[:-1],bins23qes[:-1]), bins=bins23, \
            weights = (counts23und*scale_factor,counts23coh*scale_factor,counts23dis*scale_factor,counts23res*scale_factor,counts23mec*scale_factor,counts23qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax23.set_xlabel(r"Length [cm]")
        ax23.set_title("Contained Length of Longest Proton Track in\n"+sample_title+' Events with Protons by Neutrino Interaction Mechanism')
        ax23.set_ylabel("Events / cm") 
        ax23.legend(loc='upper right')
        plt.savefig(sample_type+"_events_max_contained_proton_length_events_with_protons_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig23)

    if events_with_protons_over_thresh>0:
        # KE of Primary Protons Over Threshold
        fig24, ax24 = plt.subplots(figsize=(8,4))
        bins24 = np.linspace(0,2000,41)
        counts24und, bins24und =np.histogram(np.array(data24und), bins=bins24)
        counts24coh, bins24coh =np.histogram(np.array(data24coh), bins=bins24)
        counts24dis, bins24dis =np.histogram(np.array(data24dis), bins=bins24)
        counts24res, bins24res =np.histogram(np.array(data24res), bins=bins24)
        counts24mec, bins24mec =np.histogram(np.array(data24mec), bins=bins24)
        counts24qes, bins24qes =np.histogram(np.array(data24qes), bins=bins24)
        ax24.hist((bins24und[:-1],bins24coh[:-1],bins24dis[:-1],bins24res[:-1],bins24mec[:-1],bins24qes[:-1]), bins=bins24, \
            weights = (counts24und*scale_factor,counts24coh*scale_factor,counts24dis*scale_factor,counts24res*scale_factor,counts24mec*scale_factor,counts24qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax24.set_xlabel(r"Kinetic Energy [MeV]")
        ax24.set_title('Total Primary Proton Kinetic Energy Above Threshold\nby Neutrino Interaction Mechanism in'+sample_title+'Events with Protons')
        ax24.set_ylabel("Events / 50 MeV") 
        ax24.legend(loc='upper right')
        plt.savefig(sample_type+"_events_primary_proton_ke_truth_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig24)

        # Leading Proton Momentum 
        fig25, ax25 = plt.subplots(figsize=(8,4))
        bins25 = np.linspace(0,2000,41)
        counts25und, bins25und =np.histogram(np.array(data25und), bins=bins25)
        counts25coh, bins25coh =np.histogram(np.array(data25coh), bins=bins25)
        counts25dis, bins25dis =np.histogram(np.array(data25dis), bins=bins25)
        counts25res, bins25res =np.histogram(np.array(data25res), bins=bins25)
        counts25mec, bins25mec =np.histogram(np.array(data25mec), bins=bins25)
        counts25qes, bins25qes =np.histogram(np.array(data25qes), bins=bins25)
        ax25.hist((bins25und[:-1],bins25coh[:-1],bins25dis[:-1],bins25res[:-1],bins25mec[:-1],bins25qes[:-1]), bins=bins25, \
            weights = (counts25und*scale_factor,counts25coh*scale_factor,counts25dis*scale_factor,counts25res*scale_factor,counts25mec*scale_factor,counts25qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax25.set_xlabel(r"Momentum [MeV/c]")
        ax25.set_title("Leading Proton Momentum in "+sample_title+' Events\n with Protons Above Threshold by Neutrino Interaction Mechanism')
        ax25.set_ylabel("Events / 50 MeV/c") 
        ax25.legend(loc='upper right')
        plt.savefig(sample_type+"_events_leading_proton_momentum_truth_above_threshold_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig25)

        # Leading Proton Opening Angle
        fig26, ax26 = plt.subplots(figsize=(8,4))
        bins26 = np.linspace(0,3.2,33)
        counts26und, bins26und =np.histogram(np.array(data26und), bins=bins26)
        counts26coh, bins26coh =np.histogram(np.array(data26coh), bins=bins26)
        counts26dis, bins26dis =np.histogram(np.array(data26dis), bins=bins26)
        counts26res, bins26res =np.histogram(np.array(data26res), bins=bins26)
        counts26mec, bins26mec =np.histogram(np.array(data26mec), bins=bins26)
        counts26qes, bins26qes =np.histogram(np.array(data26qes), bins=bins26)
        ax26.hist((bins26und[:-1],bins26coh[:-1],bins26dis[:-1],bins26res[:-1],bins26mec[:-1],bins26qes[:-1]), bins=bins26, \
            weights = (counts26und*scale_factor,counts26coh*scale_factor,counts26dis*scale_factor,counts26res*scale_factor,counts26mec*scale_factor,counts26qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax26.set_xlabel(r"Angle [Rad]")
        ax26.set_title("Leading Proton Angle with respect to Beam Direction in\n"+sample_title\
                       +" Events with Protons Above Threshold Stacked by Neutrino Interaction Mechanism")
        ax26.set_ylabel("Events / 0.1 Rad") 
        ax26.legend(loc='upper right')
        plt.savefig(sample_type+"_events_lead_proton_over_thresh_angle_wrt_beam_truth_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig26)

    if events_with_geq_2_protons_over_thresh>0:
        # Sub-Leading Proton Momentum 
        fig27, ax27 = plt.subplots(figsize=(8,4))
        bins27 = np.linspace(0,2000,41)
        counts27und, bins27und =np.histogram(np.array(data27und), bins=bins27)
        counts27coh, bins27coh =np.histogram(np.array(data27coh), bins=bins27)
        counts27dis, bins27dis =np.histogram(np.array(data27dis), bins=bins27)
        counts27res, bins27res =np.histogram(np.array(data27res), bins=bins27)
        counts27mec, bins27mec =np.histogram(np.array(data27mec), bins=bins27)
        counts27qes, bins27qes =np.histogram(np.array(data27qes), bins=bins27)
        ax27.hist((bins27und[:-1],bins27coh[:-1],bins27dis[:-1],bins27res[:-1],bins27mec[:-1],bins27qes[:-1]), bins=bins27, \
            weights = (counts27und*scale_factor,counts27coh*scale_factor,counts27dis*scale_factor,counts27res*scale_factor,counts27mec*scale_factor,counts27qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax27.set_xlabel(r"Momentum [MeV/c]")
        ax27.set_title("Subleading Proton Momentum in "+sample_title+' Events\n with 2+ Protons Above Threshold by Neutrino Interaction Mechanism')
        ax27.set_ylabel("Events / 50 MeV/c") 
        ax27.legend(loc='upper right')
        plt.savefig(sample_type+"_events_sublead_proton_momentum_geq_2protons_over_thresh_truth_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig27)

        # Subleading Proton Opening Angle
        fig28, ax28 = plt.subplots(figsize=(8,4))
        bins28 = np.linspace(0,3.2,33)
        counts28und, bins28und =np.histogram(np.array(data28und), bins=bins28)
        counts28coh, bins28coh =np.histogram(np.array(data28coh), bins=bins28)
        counts28dis, bins28dis =np.histogram(np.array(data28dis), bins=bins28)
        counts28res, bins28res =np.histogram(np.array(data28res), bins=bins28)
        counts28mec, bins28mec =np.histogram(np.array(data28mec), bins=bins28)
        counts28qes, bins28qes =np.histogram(np.array(data28qes), bins=bins28)
        ax28.hist((bins28und[:-1],bins28coh[:-1],bins28dis[:-1],bins28res[:-1],bins28mec[:-1],bins28qes[:-1]), bins=bins28, \
            weights = (counts28und*scale_factor,counts28coh*scale_factor,counts28dis*scale_factor,counts28res*scale_factor,counts28mec*scale_factor,counts28qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax28.set_xlabel(r"Angle [Rad]")
        ax28.set_title("Subleading Proton Angle with respect to Beam Direction in\n"+sample_title\
                       +" Events with 2+ Protons Above Threshold Stacked by Neutrino Interaction Mechanism")
        ax28.set_ylabel("Events / 0.1 Rad") 
        ax28.legend(loc='upper right')
        plt.savefig(sample_type+"_events_sublead_proton_angle_events_w_geq_2_protons_over_thresh_wrt_beam_direction_truth_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig28)

        # Subleading proton angle wrt leading proton for events with protons over threshold
        fig29, ax29 = plt.subplots(figsize=(8,4))
        bins29 = np.linspace(0,1.6,17)
        counts29und, bins29und =np.histogram(np.array(data29und), bins=bins29)
        counts29coh, bins29coh =np.histogram(np.array(data29coh), bins=bins29)
        counts29dis, bins29dis =np.histogram(np.array(data29dis), bins=bins29)
        counts29res, bins29res =np.histogram(np.array(data29res), bins=bins29)
        counts29mec, bins29mec =np.histogram(np.array(data29mec), bins=bins29)
        counts29qes, bins29qes =np.histogram(np.array(data29qes), bins=bins29)
        ax29.hist((bins29und[:-1],bins29coh[:-1],bins29dis[:-1],bins29res[:-1],bins29mec[:-1],bins29qes[:-1]), bins=bins29, \
            weights = (counts29und*scale_factor,counts29coh*scale_factor,counts29dis*scale_factor,counts29res*scale_factor,counts29mec*scale_factor,counts29qes*scale_factor), \
            histtype='bar', label=loc_labels, stacked='True', color=['brown', 'orange', 'red', 'purple', 'blue', 'green', ])
        ax29.set_xlabel(r"Angle [Rad]")
        ax29.set_title("Subleading Proton Angle with respect to Leading Proton Direction in "+sample_title \
                       +" Events \nwith 2+ Protons Above Threshold Stacked by Neutrino Interaction Mechanism")
        ax29.set_ylabel("Events / 0.1 Rad") 
        ax29.legend(loc='upper right')
        plt.savefig(sample_type+"_events_sublead_proton_angle_wrt_lead_proton_events_w_geq_2_protons_over_thresh_stacked_by_neutrino_interaction_mechanism.png")
        plt.close(fig29)