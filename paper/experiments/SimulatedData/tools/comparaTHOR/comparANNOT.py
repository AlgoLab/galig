#!/usr/bin/python3

import argparse
import sys
import os
import csv
import glob

ANNOTATION_DIR = ""
RESULTS_DIR = ""

EVENT_TYPES = ('A5', 'A3', 'IR', 'ES')

spladder_ev_map = {'alt_5prime'       : 'A5',
                   'alt_3prime'       : 'A3',
                   'exon_skip'        : 'ES',
                   'intron_retention' : 'IR',
                   'mult_exon_skip'   : 'NA'}

rmats_ev_map = { 'A5SS' : 'A5',
                 'A3SS' : 'A3',
                 'SE'   : 'ES',
                 'RI'   : 'IR'}

suppa_ev_map = {'A5' : 'A5',
                'A3' : 'A3',
                'SE' : 'ES',
                'RI' : 'IR'}

def evaluate_asgal_annot(args, gene_events):
    asgal_introns = {ev_type : () for ev_type in EVENT_TYPES}
    with open(os.path.join(RESULTS_DIR, args.sample, args.chromosome, 'asgal_annot',
                           args.gene,   args.gene + '.events.csv'),
              newline='') as asgal_csv:
        results = csv.reader(asgal_csv)
        next(results, None) # Drop header
        for row in results:
            asgal_introns[row[0]] += ((int(row[1]) - 1, int(row[2]) + 1),)

    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(gene_events[ev_type].values()))
        tp=len(set(asgal_introns[ev_type]) & set(gene_events[ev_type].values()))
        fp=len(set(asgal_introns[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(gene_events[ev_type].values()) - set(asgal_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def evaluate_spladder_annot(args, gene_events):
    spladder_introns = {ev_type : () for ev_type in EVENT_TYPES}
    spladder_files_list = glob.glob(os.path.join(RESULTS_DIR,      args.sample,
                                                 args.chromosome,  'spladder_annot',
                                                 args.gene,        '*.txt'))
    for spladder_file in spladder_files_list:
        with open(spladder_file) as spladder_current_event_file:
            next(spladder_current_event_file, None) # Drop header
            for row in spladder_current_event_file:
                _, strand, event_id, _, *positions = row.strip('\n').split('\t')
                if event_id.startswith('mult'):
                    # here spladder outputs multiple positions in columns 2 and 3
                    # We can set them to 0 and avoid problems later on
                    positions[2], positions[3] = 0, 0
                positions = tuple(int(p) for p in positions[:6])
                spladder_introns = spladder_parse_line(spladder_introns,
                                                       event_id,
                                                       strand,
                                                       positions)
    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(gene_events[ev_type].values()))
        tp=len(set(spladder_introns[ev_type]) & set(gene_events[ev_type].values()))
        fp=len(set(spladder_introns[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(gene_events[ev_type].values()) - set(spladder_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def evaluate_rmats_annot(args, gene_events):
    # We do not consider MXE
    rmats_introns = {ev_type : () for ev_type in EVENT_TYPES}
    rmats_files_list = glob.glob(os.path.join(RESULTS_DIR,   args.sample, args.chromosome,
                                              'rMATS_annot', args.gene,
                                              'fromGTF.[!(n,M)]*.txt'))
    for rmats_file in rmats_files_list:
        with open(rmats_file) as rmats_current_event_file:
            next(rmats_current_event_file, None) # Drop header
            rmats_ev_type = rmats_file[rmats_file.rfind('.', 0, rmats_file.rfind('.')) + 1:
                                       rmats_file.rfind('.')]
            ev_type = rmats_ev_map[rmats_ev_type]
            for row in rmats_current_event_file:
                _, _, _, _, strand, *positions = row.strip('\n').split('\t')
                positions = tuple(int(p) for p in positions)
                rmats_introns = rmats_parse_line(rmats_introns,
                                                 ev_type,
                                                 strand,
                                                 positions)
    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(gene_events[ev_type].values()))
        tp=len(set(rmats_introns[ev_type]) & set(gene_events[ev_type].values()))
        fp=len(set(rmats_introns[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(gene_events[ev_type].values()) - set(rmats_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def suppa_parse_line(suppa_introns, ev_type, strand, exon_list):
    if ev_type is 'IR':
        p1, p2 = tuple(int(p) for p in exon_list[1].split('-'))
        suppa_introns[ev_type] += ((p1, p2),)
    elif ev_type is 'ES':
        p = tuple(int(p) for exon in exon_list for p in exon.split('-'))
        suppa_introns[ev_type] += ((p[0], p[3]),)
    elif ev_type is 'IR':
        p = tuple(int(p) for exon in exon_list for p in exon.split('-'))
        suppa_introns[ev_type] += ((p[1], p[2]),)
    elif ev_type in ['A3', 'A5']:
        p = tuple(int(p) for exon in exon_list for p in exon.split('-'))
        suppa_introns[ev_type] += ((p[0], p[1]),
                                   (p[2], p[3]))
    return suppa_introns

def evaluate_suppa_annot(args, gene_events):
    # We do not consider: MX (mutually exclusive), AF (alternative first exon), AL (alternative last exon)
    suppa_introns = {ev_type : () for ev_type in EVENT_TYPES}
    with open(os.path.join(RESULTS_DIR, args.sample, args.chromosome, 'suppa_annot',
                           args.gene,   'iso_tpm.psi')) as suppa_file:
        next(suppa_file, None) # Drop header
        for line in suppa_file:
            transcript, rest = line.split('\t')[0].split(';')
            event_id, chr, *exon_list, strand = rest.split(':')
            if event_id in ['MX', 'AF', 'AL']:
                continue
            ev_type = suppa_ev_map[event_id]
            suppa_introns = suppa_parse_line(suppa_introns,
                                             ev_type,
                                             strand,
                                             exon_list)
    results = {ev_type : {} for ev_type in EVENT_TYPES}
    for ev_type in EVENT_TYPES:
        nelems = len(set(gene_events[ev_type].values()))
        tp=len(set(suppa_introns[ev_type]) & set(gene_events[ev_type].values()))
        fp=len(set(suppa_introns[ev_type]) - set(gene_events[ev_type].values()))
        fn=len(set(gene_events[ev_type].values()) - set(suppa_introns[ev_type]))
        results[ev_type] = {'nelems' : nelems, 'tp' : tp, 'fp' : fp, 'fn' : fn}
    return results

def spladder_parse_line(spladder_introns, event_id, strand, positions):
    ev_type = spladder_ev_map[event_id[0:event_id.rfind('_')]]
    if ev_type in ['ES', 'IR']:
        # Exon skipping, mult exon skipping, intron retention
        spladder_introns[ev_type] += ((positions[1], positions[4]),)
    elif ev_type is 'A5':
        if strand is '+':
            spladder_introns[ev_type] += ((positions[3], positions[0]),
                                          (positions[5], positions[0]))
        elif strand is '-':
            spladder_introns[ev_type] += ((positions[1], positions[2]),
                                          (positions[1], positions[4]))
    elif ev_type is 'A3':
        if strand is '+':
            spladder_introns[ev_type] += ((positions[1], positions[2]),
                                          (positions[1], positions[4]))
        elif strand is '-':
            spladder_introns[ev_type] += ((positions[3], positions[0]),
                                          (positions[5], positions[0]))
    return spladder_introns

def rmats_parse_line(rmats_introns, ev_type, strand, positions):
    if ev_type in ['ES', 'IR']:
        # Exon skipping, mult exon skipping, intron retention
        rmats_introns[ev_type] += ((positions[3], positions[4] + 1),)
    elif ev_type is 'A5':
        if strand is '+':
            rmats_introns[ev_type] += ((positions[1], positions[4] + 1),
                                       (positions[3], positions[4] + 1))
        elif strand is '-':
            rmats_introns[ev_type] += ((positions[5], positions[0] + 1),
                                       (positions[5], positions[2] + 1))
    elif ev_type is 'A3':
        if strand is '+':
            rmats_introns[ev_type] += ((positions[5], positions[0] + 1),
                                       (positions[5], positions[2] + 1))
        elif strand is '-':
            rmats_introns[ev_type] += ((positions[1], positions[4] + 1),
                                       (positions[3], positions[4] + 1))
    return rmats_introns

def add_alternative_donor_acceptor(gene_events, ev_type, ev_num, positions):
    if positions[0] < positions[1] < positions[2]:
        gene_events[ev_type][(ev_type, ev_num, 'i')] = tuple(sorted((positions[0], positions[2])))
        gene_events[ev_type][(ev_type, ev_num, 'e')] = tuple(sorted((positions[1], positions[2])))
    else:
        gene_events[ev_type][(ev_type, ev_num, 'e')] = tuple(sorted((positions[0], positions[2])))
        gene_events[ev_type][(ev_type, ev_num, 'i')] = tuple(sorted((positions[1], positions[2])))
    return gene_events

def parse_gene_events(args):
    gene_events = {ev_type : {} for ev_type in EVENT_TYPES}
    with open(os.path.join(ANNOTATION_DIR,
                           args.chromosome,
                           args.gene + '.events')) as gene_events_file:
        for line in gene_events_file:
            _, ev_type, *positions, ev_id = line.strip('\n').split(' ')
            positions = tuple(int(p) for p in positions)
            _, ev_num = ev_id[:2], ev_id[2:]
            if not ev_type.startswith('A'):
                gene_events[ev_type][(ev_type, ev_num)] = (positions[0], positions[1])
            else:
                gene_events = add_alternative_donor_acceptor(gene_events, ev_type,
                                                             ev_num, positions)
    return gene_events

def main():
    global ANNOTATION_DIR, RESULTS_DIR

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-f', '--folder',      dest='base_folder', help='Data folder',     required=True)
    parser.add_argument('-c', '--chr',         dest='chromosome',  help='Chromosome name', required=True)
    parser.add_argument('-g', '--gene',        dest='gene',        help='Gene name',       required=True)
    parser.add_argument('-s', '--sample-size', dest='sample',      help='Sample Size',     required=True)

    args = parser.parse_args()

    ANNOTATION_DIR = os.path.join(args.base_folder, 'ReducedAnnotations')
    RESULTS_DIR = os.path.join(args.base_folder, 'Results')
    
    gene_events = parse_gene_events(args)

    print("SIZE,GENOME,GENE,TOOL,EV_TYPE,NELEMS,TP,FP,FN")

    results = evaluate_asgal_annot(args, gene_events)
    for ev_type in EVENT_TYPES:
        print("{},{},{},asgal-annot,{},{},{},{},{}".format(args.sample, args.chromosome,
                                                           args.gene,
                                                           ev_type,
                                                           results[ev_type]['nelems'],
                                                           results[ev_type]['tp'],
                                                           results[ev_type]['fp'],
                                                           results[ev_type]['fn']))

    results = evaluate_spladder_annot(args, gene_events)
    for ev_type in EVENT_TYPES:
        print("{},{},{},spladder-annot,{},{},{},{},{}".format(args.sample, args.chromosome,
                                                              args.gene,
                                                              ev_type,
                                                              results[ev_type]['nelems'],
                                                              results[ev_type]['tp'],
                                                              results[ev_type]['fp'],
                                                              results[ev_type]['fn']))

    results =  evaluate_rmats_annot(args, gene_events)
    for ev_type in EVENT_TYPES:
        print("{},{},{},rmats-annot,{},{},{},{},{}".format(args.sample, args.chromosome,
                                                           args.gene,
                                                           ev_type,
                                                           results[ev_type]['nelems'],
                                                           results[ev_type]['tp'],
                                                           results[ev_type]['fp'],
                                                           results[ev_type]['fn']))

    results = evaluate_suppa_annot(args, gene_events)
    for ev_type in EVENT_TYPES:
        print("{},{},{},suppa-annot,{},{},{},{},{}".format(args.sample, args.chromosome,
                                                           args.gene,
                                                           ev_type,
                                                           results[ev_type]['nelems'],
                                                           results[ev_type]['tp'],
                                                           results[ev_type]['fp'],
                                                           results[ev_type]['fn']))


if __name__ == "__main__":
    main()

