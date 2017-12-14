#!/usr/bin/env python

""" MultiQC module to parse output from Salmon """

from __future__ import print_function
from collections import OrderedDict
from GCModel import GCModel
from SeqModel import SeqModel
import json
import logging
import os
import sys
import numpy

from multiqc import config
from multiqc.plots import linegraph
from multiqc.plots import heatmap
from multiqc.modules.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(name='Salmon', anchor='salmon',
                href='http://combine-lab.github.io/salmon/',
                info="is a tool for quantifying the expression of transcripts using RNA-seq data.")

        #Intializing Variables
        self.salmon_meta = dict()
        self.salmon_gc_lower = {}
        self.salmon_gc_middle = {}
        self.salmon_gc_upper = {}
        self.salmon_gc_average = {}
        self.salmon_seq_3_a = {}
        self.salmon_seq_3_c = {}
        self.salmon_seq_3_g = {}
        self.salmon_seq_3_t = {}
        self.salmon_seq_5_a = {}
        self.salmon_seq_5_c = {}
        self.salmon_seq_5_g = {}
        self.salmon_seq_5_t = {}
        self.salmon_seq_5_avg = {}
        self.salmon_seq_3_avg = {}
        self.heat_map_gc_bias_lower_values = []
        self.heat_map_gc_bias_middle_values = []
        self.heat_map_gc_bias_upper_values = []
        self.heat_map_gc_bias_average_values = []
        self.heat_map_seq_3_values = []
        self.heat_map_seq_5_values = []
        self.sample_names = []

        # Parse meta information. JSON win!
        for f in self.find_log_files('salmon/meta'):
            # Get the s_name from the parent directory
            s_name = os.path.basename( os.path.dirname(f['root']) )
            s_name = self.clean_s_name(s_name, f['root'])
            sample_name = os.path.basename(os.path.dirname(os.path.dirname(f['root'])))
            if s_name=='bias':
                meta_info_json=json.loads(f['f'])
                if meta_info_json['gc_bias_correct'] == True:

                    #######################################################
                    # Calculating GC Bias
                    self.sample_names.append(sample_name)

                    gcModel = GCModel()
                    gcModel.from_file(os.path.dirname(f['root']))
                    gc_bias = numpy.divide(gcModel.obs_,gcModel.exp_).tolist()
                    weights = numpy.divide(gcModel.obs_weights_,gcModel.exp_weights_).tolist()
                    average_values = OrderedDict()

                    #Lower section of reads
                    lower_values = OrderedDict()
                    for j in range(len(gc_bias[0])):
                        lower_values[j*4.45] = gc_bias[0][j] * weights[0]
                        average_values[j*4.45] = average_values.get(j,0) + lower_values[j*4.45]
                    self.salmon_gc_lower[sample_name] = lower_values

                    #Middle section of reads
                    middle_values = OrderedDict()
                    for j in range(len(gc_bias[1])):
                        middle_values[j*4.45] = gc_bias[1][j] * weights[1]
                        average_values[j*4.45] = average_values.get(j, 0) + middle_values[j*4.45]
                    self.salmon_gc_middle[sample_name] = middle_values

                    #Upper section of reads
                    upper_values = OrderedDict()
                    for j in range(len(gc_bias[2])):
                        upper_values[j*4.45] = gc_bias[2][j] * weights[2]
                        average_values[j*4.45] = average_values.get(j, 0) + upper_values[j*4.45]
                    self.salmon_gc_upper[sample_name] = upper_values

                    #Averaging the values over all three section of reads
                    for j in range(len(gc_bias[0])):
                        average_values[j*4.45] /= 3


                    self.heat_map_gc_bias_lower_values.append(lower_values.values())
                    self.heat_map_gc_bias_middle_values.append(middle_values.values())
                    self.heat_map_gc_bias_upper_values.append(upper_values.values())
                    self.heat_map_gc_bias_average_values.append(average_values.values()) 
                    self.salmon_gc_average[sample_name] = average_values

                    #############################################################
                    # Calculating Sequence Bias
                    seqModel = SeqModel()
                    seqModel.from_file(os.path.dirname(f['root']))
                    seq_bias_3 = numpy.divide(seqModel.obs3_, seqModel.exp3_).tolist()
                    seq_bias_5 = numpy.divide(seqModel.obs5_, seqModel.exp5_).tolist()
                    average_3_values = OrderedDict()
                    average_5_values = OrderedDict()

                    seq_bias_3_a = OrderedDict()
                    seq_bias_3_c = OrderedDict()
                    seq_bias_3_g = OrderedDict()
                    seq_bias_3_t = OrderedDict()
                    for j in range(len(seq_bias_3[0])):
                        seq_bias_3_a[j] = seq_bias_3[0][j]
                        seq_bias_3_c[j] = seq_bias_3[1][j]
                        seq_bias_3_g[j] = seq_bias_3[2][j]
                        seq_bias_3_t[j] = seq_bias_3[3][j]

                    self.salmon_seq_3_a[sample_name] = seq_bias_3_a
                    self.salmon_seq_3_c[sample_name] = seq_bias_3_c
                    self.salmon_seq_3_g[sample_name] = seq_bias_3_g
                    self.salmon_seq_3_t[sample_name] = seq_bias_3_t

                    seq_bias_5_a = OrderedDict()
                    seq_bias_5_c = OrderedDict()
                    seq_bias_5_g = OrderedDict()
                    seq_bias_5_t = OrderedDict()
                    for j in range(len(seq_bias_5[0])):
                        seq_bias_5_a[j] = seq_bias_5[0][j]
                        seq_bias_5_c[j] = seq_bias_5[1][j]
                        seq_bias_5_g[j] = seq_bias_5[2][j]
                        seq_bias_5_t[j] = seq_bias_5[3][j]

                    self.salmon_seq_5_a[sample_name] = seq_bias_5_a
                    self.salmon_seq_5_c[sample_name] = seq_bias_5_c
                    self.salmon_seq_5_g[sample_name] = seq_bias_5_g
                    self.salmon_seq_5_t[sample_name] = seq_bias_5_t

                    for i in range(len(seq_bias_3[0])):
                        for j in range(0, 4):
                            average_3_values[i] = average_3_values.get(i, 0) + seq_bias_3[j][i]
                        average_3_values[i] /= 4

                    for i in range(len(seq_bias_5[0])):
                        for j in range(0, 4):
                            average_5_values[i] = average_5_values.get(i, 0) + seq_bias_5[j][i]
                        average_5_values[i] /= 4

                    self.salmon_seq_3_avg[sample_name] = average_3_values
                    self.salmon_seq_5_avg[sample_name] = average_5_values
                    self.heat_map_seq_3_values.append(average_3_values.values()) 
                    self.heat_map_seq_5_values.append(average_5_values.values()) 


            self.salmon_meta[sample_name] = json.loads(f['f'])

        ############### Configs for distribution plots #################
        pconfig_GCBias_Begin = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: GC Bias in Beginning of Reads',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'xCeiling':100,
            'xFloor':0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_GCBias_Middle = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: GC Bias in Middle of Reads',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_GCBias_Last = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: GC Bias in Last of Reads',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_GCBias_Average = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: GC Bias Average Distribution for All Samples',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_3_A = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 3 A Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_3_C = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 3 C Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_3_G = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 3 G Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_3_T = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 3 T Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_5_A = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 5 A Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_5_C = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 5 C Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_5_G = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 5 G Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_5_T = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 5 T Base',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_3_Average = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 3 Average',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }
        pconfig_SeqBias_5_Average = {
            'smooth_points': 500,
            'id': 'salmon_plot',
            'title': 'Salmon: Sequence Bias 5 Average',
            'ylab': 'Ratio',
            'xlab': 'Bias',
            'ymin': 0,
            'xmin': 0,
            'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
        }

        ########### Distribution Plots ##################################
        self.add_section(name='GC Bias First Row', plot = linegraph.plot(self.salmon_gc_lower, pconfig_GCBias_Begin))
        self.add_section(name='GC Bias Middle Row', plot=linegraph.plot(self.salmon_gc_middle, pconfig_GCBias_Middle))
        self.add_section(name='GC Bias Last Row', plot=linegraph.plot(self.salmon_gc_upper, pconfig_GCBias_Last))
        self.add_section(name='GC Bias Average', plot=linegraph.plot(self.salmon_gc_average, pconfig_GCBias_Average))
        self.add_section(name='Sequence 3 A-Base', plot=linegraph.plot(self.salmon_seq_3_a, pconfig_SeqBias_3_A))
        self.add_section(name='Sequence 3 C-Base', plot=linegraph.plot(self.salmon_seq_3_c, pconfig_SeqBias_3_C))
        self.add_section(name='Sequence 3 G-Base', plot=linegraph.plot(self.salmon_seq_3_g, pconfig_SeqBias_3_G))
        self.add_section(name='Sequence 3 T-Base', plot=linegraph.plot(self.salmon_seq_3_t, pconfig_SeqBias_3_T))
        self.add_section(name='Sequence 3 Average', plot=linegraph.plot(self.salmon_seq_3_avg, pconfig_SeqBias_3_Average))
        self.add_section(name='Sequence 5 A-Base', plot=linegraph.plot(self.salmon_seq_5_a, pconfig_SeqBias_5_A))
        self.add_section(name='Sequence 5 C-Base', plot=linegraph.plot(self.salmon_seq_5_c, pconfig_SeqBias_5_C))
        self.add_section(name='Sequence 5 G-Base', plot=linegraph.plot(self.salmon_seq_5_g, pconfig_SeqBias_5_G))
        self.add_section(name='Sequence 5 T-Base', plot=linegraph.plot(self.salmon_seq_5_t, pconfig_SeqBias_5_T))
        self.add_section(name='Sequence 5 Average', plot=linegraph.plot(self.salmon_seq_5_avg, pconfig_SeqBias_5_Average))

        ################### HeatMaps #####################################
        self.add_section(name='GC Bias First Row Heatmap', plot=heatmap.plot(numpy.corrcoef(self.heat_map_gc_bias_lower_values),self.sample_names, self.sample_names))
        self.add_section(name='GC Bias Middle Row Heatmap', plot=heatmap.plot(numpy.corrcoef(self.heat_map_gc_bias_middle_values),self.sample_names, self.sample_names))
        self.add_section(name='GC Bias Last Row Heatmap', plot=heatmap.plot(numpy.corrcoef(self.heat_map_gc_bias_upper_values),self.sample_names, self.sample_names))
        self.add_section(name='GC Bias Average Heatmap', plot=heatmap.plot(numpy.corrcoef(self.heat_map_gc_bias_average_values),self.sample_names, self.sample_names))
        self.add_section(name='Sequence 3 Heatmap', plot=heatmap.plot(numpy.corrcoef(self.heat_map_seq_3_values),self.sample_names, self.sample_names))
        self.add_section(name='Sequence 5 Heatmap', plot=heatmap.plot(numpy.corrcoef(self.heat_map_seq_5_values),self.sample_names, self.sample_names))

        # Parse Fragment Length Distribution logs
        self.salmon_fld = dict()
        for f in self.find_log_files('salmon/fld'):
            # Get the s_name from the parent directory
            if os.path.basename(f['root']) == 'libParams':
                s_name = os.path.basename( os.path.dirname(f['root']) )
                s_name = self.clean_s_name(s_name, f['root'])
                parsed = OrderedDict()
                for i, v in enumerate(f['f'].split()):
                    parsed[i] = float(v)
                if len(parsed) > 0:
                    if s_name in self.salmon_fld:
                        log.debug("Duplicate sample name found! Overwriting: {}".format(s_name))
                    self.add_data_source(f, s_name)
                    self.salmon_fld[s_name] = parsed

        # Filter to strip out ignored sample names
        self.salmon_meta = self.ignore_samples(self.salmon_meta)
        self.salmon_fld = self.ignore_samples(self.salmon_fld)

        if len(self.salmon_meta) == 0 and len(self.salmon_fld) == 0:
            raise UserWarning

        if len(self.salmon_meta) > 0:
            log.info("Found {} meta reports".format(len(self.salmon_meta)))
            self.write_data_file(self.salmon_meta, 'multiqc_salmon')
        if len(self.salmon_fld) > 0:
            log.info("Found {} fragment length distributions".format(len(self.salmon_fld)))

        # Add alignment rate to the general stats table
        headers = OrderedDict()
        headers['percent_mapped'] = {
                'title': '% Aligned',
                'description': '% Mapped reads',
                'max': 100,
                'min': 0,
                'suffix': '%',
                'scale': 'YlGn'
                }
        headers['num_mapped'] = {
                'title': 'M Aligned',
                'description': 'Mapped reads (millions)',
                'min': 0,
                'scale': 'PuRd',
                'modify': lambda x: float(x) / 1000000,
                'shared_key': 'read_count'
                }
        self.general_stats_addcols(self.salmon_meta, headers)

        # Fragment length distribution plot
        pconfig = {
                'smooth_points': 500,
                'id': 'salmon_plot',
                'title': 'Salmon: Fragment Length Distribution',
                'ylab': 'Fraction',
                'xlab': 'Fragment Length (bp)',
                'ymin': 0,
                'xmin': 0,
                'tt_label': '<b>{point.x:,.0f} bp</b>: {point.y:,.0f}',
                }
        self.add_section( plot = linegraph.plot(self.salmon_fld, pconfig) )


        
        


