import argparse
import os
import sys
import subprocess
import re

def netmhc_run(netmhcpath, peppath, alleles, outpath):
    """
    :param netmhcpath: absolute path of the netmhc execution file
    :param peppath: input peptides file
    :param alleles: HLA alleles separated by ,
    :param outpath: outfile
    :return: None
    """

    cmd = netmhcpath \
        + ' -a ' + alleles \
        + ' -f ' + peppath \
        + ' -inptype 1' \
        + ' > ' + outpath
    p = subprocess.Popen(cmd, shell=True)
    p.wait()


def netmhc_reload(resultpaths):
    """
    :param resultpath: the result of netmhc
    :return: a dictionary {neoepitope: {allele: [affinity, rank, FILTER]}}
    """
    pep_dic = {}
    for resultpath in resultpaths:
        with open(resultpath, 'r') as f:
            for line in f:
                if line.startswith("    0  HLA"):
                    tmpline = re.split(r'\s+', line)
                    allele = tmpline[2]
                    pep = tmpline[3]
                    affinity = float(tmpline[6])
                    rank = float(tmpline[7])
                    if pep not in pep_dic:
                        pep_dic[pep] = {}
                    pep_dic[pep][allele] = [affinity, rank, 'FILTER']
    return pep_dic


def create_arg_parser():
    parser = argparse.ArgumentParser(prog="netMHCstabpan")
    parser.add_argument('-a', '--allele-file', dest='afile', metavar='ALLELE_FILE', required=True,
                        help='File for HLA allele.')
    parser.add_argument('-p', '--pep-file', dest='pfile', metavar='PEP_FILE', default=None,
                        help='File for a list of peptides')
    parser.add_argument('-o', '--out-dir', dest='odir', metavar='OUTPUT_DIR', default=None,
                        help='Output directory')
    args = parser.parse_args()
    return args


def hla_load(filepath):
    """
    :param filepath: the absolute path of a HLA typing file
    :return: a list of HLA alleles joined by ,
    """
    hla_alleles = []
    filename = os.path.basename(filepath)
    with open(filepath, 'r') as f:
        for line in f:
            hla_allele = line.rstrip()
            if not hla_format_check(hla_allele):
                raise IOError("{0} in file {1} is not a supported HLA format.".format(hla_allele, filename))
            hla_alleles.append(hla_allele.replace('*', ''))
    return ','.join(hla_alleles)


def hla_format_check(hla_allele):
    """
    :param hla_allele: a string indicating the hla allele
    :return: whether it is in right format
    """
    legal_pattern = re.compile(r'HLA-[ABC]\*\d{2}:\d{2}$')
    return legal_pattern.match(hla_allele)


def pep_split(peppath, outdir):
    peps = []
    with open(peppath, 'r') as f:
        for line in f:
            peps.append(line.rstrip())
    len2pep = {}
    for pep in peps:
        length = len(pep)
        if length in len2pep:
            len2pep[length].append(pep)
        else:
            len2pep[length] = [pep]
    
    pepfiles = []
    for length in len2pep:
        pepfile = 'pep.' + str(length) + '.txt'
        pepfiles.append(pepfile)
        with open(os.path.join(outdir, pepfile), 'w') as f:
            for pep in len2pep[length]:
                f.write(pep + '\n')

    return pepfiles


def main():
    args = create_arg_parser()
    if not os.path.exists(args.odir):
        os.mkdir(args.odir)

    netmhcpath = '/lustre1/shiyang/software/netMHCstabpan-1.0/netMHCstabpan'
    alleles = hla_load(args.afile)
    pepfiles = pep_split(args.pfile, args.odir)
    opaths = []
    for pepfile in pepfiles:
        peppath = os.path.join(args.odir, pepfile)
        opath = os.path.join(args.odir, pepfile.replace('pep', 'neo'))
        opaths.append(opath)
        netmhc_run(netmhcpath, peppath, alleles, opath)

    dic_neo = netmhc_reload(opaths)
    print(dic_neo)
    with open(os.path.join(args.odir, 'neo.all.txt'), 'w') as f:
        f.write('Neoantigen' + '\t' + 'Allele' + '\t' + 'Stability' + '\t' + 'Rank_stab' + '\n')
        for pep in dic_neo:
            for allele in dic_neo[pep]:
                stab = str(dic_neo[pep][allele][0])
                rank = str(dic_neo[pep][allele][1])
                f.write('\t'.join([pep, allele, stab, rank]) + '\n')


if __name__ == '__main__':
    main()
