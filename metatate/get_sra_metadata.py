#! /usr/bin/env python

import sys
import os
import re
import xml.etree.ElementTree as ET
from xml.dom import minidom

from Bio import Entrez

def prettify(elem):
    """Return a pretty-printed XML string for the Element.
    """
    rough_string = ET.tostring(elem, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    return reparsed.toprettyxml(indent="\t").encode('utf-8')

def get_rundata(xmlfile):
    rundata = {}
    tree = ET.parse(xmlfile)
    root = tree.getroot()
    
    samp = root.findall('./EXPERIMENT_PACKAGE/SAMPLE')
    assert len(samp)==1
    samp = samp[0]
    rundata['sample_accession'] = samp.get('accession')
    rundata['sample_alias'] = samp.get('alias')
    if samp.find('TITLE') is not None:
        rundata['sample_title'] = samp.find('TITLE').text
    for sattr in samp.findall("SAMPLE_ATTRIBUTES/SAMPLE_ATTRIBUTE"):
        rundata['sample_%s' % sattr.find("TAG").text] = sattr.find("VALUE").text
    
    run = root.findall('./EXPERIMENT_PACKAGE/RUN_SET/RUN')
    assert len(run)==1
    run = run[0]
    rundata['run_accession'] = run.get('accession')
    rundata['run_alias'] = run.get('alias')
    rundata['run_bases'] = run.find('Bases').get('count')
    rundata['run_readperspot'] = run.find('Statistics').get('nreads')
    rundata['run_nspots'] = run.find('Statistics').get('nspots')
    
    experiment = root.findall('./EXPERIMENT_PACKAGE/EXPERIMENT')
    assert len(experiment) == 1
    experiment = experiment[0]
    rundata['experiment_accession'] = experiment.get('accession')
    rundata['experiment_alias'] = experiment.get('alias')
    libdesc = experiment.find('./DESIGN/LIBRARY_DESCRIPTOR')
    if libdesc is not None:
        for child in libdesc:
            if not list(child):
                rundata[child.tag.lower()] = child.text
            elif len(list(child))==1:
                rundata[child.tag.lower()] = list(child)[0].tag
            else:
                print >>sys.stderr, "WARNING: Ambiguous tag %s" % child.tag
                rundata[child.tag.lower()] = list(child)[0].tag
    
    return rundata

def main(args):
    if args.email is None: sys.exit("Email is required for using Entrez")
    Entrez.email = args.email
    sraID = args.sra_id
    if not os.path.exists(args.dest):
        os.mkdir(args.dest)
    
    # Search entrez for the SRA study ID
    ### This stopped working for some reason:
    # record = Entrez.read(Entrez.esearch(db="sra",term=sraID,retmax=5000))
    # id_list = record['IdList']
    ### Here is a workaround:
    response = Entrez.esearch(db="sra",term=sraID,retmax=5000)
    searchroot = ET.fromstring(response.read())
    id_list = [elem.text for elem in searchroot.findall('IdList/Id')]
    
    # Download the XML for each sample
    xml_files = []
    for sampid in id_list:
        print >>sys.stderr, 'Fetching %s' % sampid
        pr = Entrez.efetch(db="sra",id=sampid).read()
        root = ET.fromstring(pr)
        run_ids = root.findall('./EXPERIMENT_PACKAGE/RUN_SET/RUN/IDENTIFIERS/PRIMARY_ID')
        print >>sys.stderr, 'Runs: %s' % ', '.join(_.text for _ in run_ids)
        xml_files.append('%s/%s.xml' % (args.dest, run_ids[0].text))
        with open(xml_files[-1], 'w') as outh:
            print >>outh, prettify(root)
    
    # Parse XML for summary table
    run_summary = []
    for f in xml_files:
        run_summary.append(get_rundata(f))
    
    # Output the summary table
    columns = sorted(run_summary[0].keys())
    with open('%s/sample_matrix.txt' % args.dest,'w') as outh:
        print >>outh, '\t'.join(_.replace(' ','_') for _ in columns)
        for d in run_summary:
            print >>outh, '\t'.join(d[c] for c in columns)

    # Output URLs to the SRA files
    with open('%s/sample_urls.txt' % args.dest, 'w') as outh:
        for d in run_summary:
            run_accn = d['run_accession']        
            print >>outh, "ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByRun/sra/%s/%s/%s/%s.sra" %  (run_accn[:3],run_accn[:6],run_accn,run_accn)

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description='Get information for SRA project.')
    parser.add_argument('--email', help='Email address, needed for using Entrez')    
    parser.add_argument('--dest', default='metadata', help='Destination directory')   
    parser.add_argument('sra_id')
    main(parser.parse_args())
