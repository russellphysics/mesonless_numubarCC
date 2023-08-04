################################################################################
##                                                                            ##
##    CONTAINS: Script to create plots describing muons and hadrons in signal ##
##              events, using muon and hadron dictionaries created using the  ##
##              methods in /truth_kinematics/file_parsing/                    ##
##              signal_characterization.py and a scale factor for scaling     ##
##              event counts to those expected with 2.5e19 POT.               ##
##                                                                            ##
################################################################################

import json
import argparse
from plot_signal_muons import plot_muons
from plot_signal_hadrons import plot_hadrons

def main(scale_factor, muon_json_file, hadron_json_file):

    muon_file = open(muon_json_file)
    muon_dict=json.load(muon_file)

    hadron_file = open(hadron_json_file)
    hadron_dict=json.load(hadron_file)

    plot_muons(muon_dict, scale_factor, sig_bkg = 0)
    plot_hadrons(hadron_dict, scale_factor, sig_bkg = 0)


if __name__=='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sf', '--scale_factor', default=1, required=True, type=float, \
                        help='''Scale factor related to input dictionary''')
    parser.add_argument('-mu', '--muon_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the muon info JSON file''')
    parser.add_argument('-had', '--hadron_json_file', default=None, required=True, type=str, \
                        help='''string corresponding to the path of the hadron info JSON file''')
    args = parser.parse_args()
    main(**vars(args))
