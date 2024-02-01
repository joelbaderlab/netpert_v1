#!/usr/bin/env python

import sys
import os.path
import argparse
import logging
import pandas as pd
import numpy as np
from math import log, e, isclose
from scipy import linalg,integrate
from collections import defaultdict

logging.basicConfig(format='%(name)s %(levelname)s: %(message)s')
logger = logging.getLogger('NetPert')
logger.setLevel(logging.INFO)

class CMDParser():
    #Command line parser

    def __init__(self):
        parser = argparse.ArgumentParser()
        subparser = parser.add_subparsers()

        #NetPert analysis
        analysis_parser = subparser.add_parser('analysis', help='Run NetPert.')
        analysis_parser.add_argument('project_dir',help='Project directory path.')
        analysis_parser.add_argument('driver',help='Network driver gene symbol.')
        analysis_parser.add_argument('-s',choices=['human','mouse'],default='human',help='Species. Default human')
        analysis_parser.set_defaults(func = runAnalysis)

        args = parser.parse_args()
        self.__dict__.update(args.__dict__)
        return

def runAnalysis(args,output='write'):
    """
    Runs NetPert analysis
    """
    print('*************************************************')
    logger.info('Running NetPert analysis...')
    print('*************************************************')
    
    #Compile interactions into one file and store it in databases dir, if needed.
    getNetworkInteractions(args)

    #Get driver, intermediates, and network respones.
    driver,intermediates,responses = getNetworkGenes(args,verbose=True)

    #Prioritize genes
    prioritization_parts = getPrioritization(args,driver,intermediates,responses,r=0.5,nt=2)

    #Output
    if output == 'write':
        #write to file
        writeWeights(args,driver,intermediates,responses,prioritization_parts,drug_target_max=5)
        return None
    elif output == 'return':
        #return results
        return(prioritization_parts)
    else:
        logger.error('runAnalysis returns NetPert prioritization results or writes to file')
        sys.exit()

def writeWeights(args,driver,intermediates,responses,prioritization_parts,drug_target_max=5):
    """
    Writes NetPert weights to project directory.
    Adds Repurposing Hub compounds.
    """

    v2weight = prioritization_parts['NetPert']['weights']
    v2rank = prioritization_parts['NetPert']['ranks']
    vertices_ranked = prioritization_parts['NetPert']['vertices_ranked']

    #map human gene to repurposing hub drug
    hgene2repurp,repurp2hgene = getRepurpDrugs()

    #map mouse gene to repurposing hub drug
    v2repurp = getNetworkGenes2RepurpDrugs(args,driver,intermediates,responses,hgene2repurp)

    #get gene pathway type
    v2pathType = getPathType(args,driver,intermediates,responses)
    v2DIIR = getPathTypeDIIR(args,driver,intermediates,responses,v2pathType)

    #write
    filepath = args.project_dir + '/geneRankings.tsv'
    fp = open(filepath,'w')
    fields = ['gene','type','DIIR','NetPert_rank','NetPert_weight','repurposinghub_cpds']
    fp.write('\t'.join(fields) + '\n')
    for v in vertices_ranked:
        repurps = dict()
        repurplist = list()
        if v in v2repurp:
            repurps = v2repurp[v]
        for drug in repurps.keys():
            l = len(repurp2hgene[drug])
            if l <= drug_target_max:
                repurplist.append(drug)
        repurpstr = ','.join(repurplist)

        if v2DIIR[v]:
            diir = 'DIIR'
        else:
            diir = ''

        toks = [v,v2pathType[v],diir,str(v2rank[v]),"{:.2e}".format(v2weight[v]),repurpstr]
        fp.write('\t'.join(toks) + '\n')
    fp.close()
    logger.info('Wrote gene rankings to %s.',filepath)
    return None

def getPathTypeDIIR(args,driver,intermediates,responses,v2pathType):
    """
    Returns dict mapping of gene to boolean:
    True = gene is a DIIR gene
    False = gene is not a DIIR gene
    """

    ret = dict()

    #get network edges
    edge_dict = readEdges(args)
    edge_ab,edge_ba = edgeTuplesToDicts(edge_dict)

    #DIIR genes are one of the following:
    #DI -> DIR
    #DI -> IR (both are DIIR)
    #DIR -> IR

    for v in v2pathType:
        ret[v] = False
        typ = v2pathType[v]

        if typ == 'DI' and v in edge_ab:
            for des in edge_ab[v]:
                des_typ = v2pathType[des]
                if des_typ == 'DIR' or des_typ == 'IR':
                    ret[v] = True
                    continue
        elif typ == 'IR' and v in edge_ba:
            for src in edge_ba[v]:
                src_typ = v2pathType[src]
                if src_typ == 'DIR' or src_typ == 'DI':
                    ret[v] = True
                    continue

    return(ret)


def getPathType(args,driver,intermediates,responses):
    """
    Returns dict mapping gene to pathway type:
    D = driver
    R = response
    DIR = intermediate with edge from driver and edge to response
    DI = intermediate with edge from driver
    IR = intermediate with edge to response
    I = intermediate no direct edges to/from driver/response
    """

    ret = dict()

    #driver
    ret[driver] = 'D'

    #responses
    for r in responses:
        ret[r] = 'R'

    #intermediates
    #get network edges
    edge_dict = readEdges(args)
    edge_ab,edge_ba = edgeTuplesToDicts(edge_dict)

    for v in intermediates:
        ret[v] = 'I'
        #find DI genes
        if (v in edge_ba) and (driver in edge_ba[v].keys()):
            ret[v] = 'DI'
            #find DIR genes
            if v in edge_ab:
                for des in edge_ab[v].keys():
                    if des in responses:
                        ret[v] = 'DIR'
            continue
        #find IR genes
        elif v in edge_ab:
            for des in edge_ab[v].keys():
                if des in responses:
                    ret[v] = 'IR'

    return(ret)

def getNetworkGenes2RepurpDrugs(args,driver,intermediates,responses,hgene2repurp):
    """
    Returns a dict of dict's mapping any gene to drug repurposing compounds targeting the protein of the human gene.
    """

    g2repurp = defaultdict(dict)

    #get mouse to human conversion, if needed
    if args.s == 'mouse':
        mouse2human,_ = getMouse2Human()
    else:
        mouse2human = dict()
    
    allgenes = list()
    allgenes.append(driver)
    allgenes.extend(intermediates)
    allgenes.extend(responses)

    for gene in allgenes:
        if gene in mouse2human.keys():
            gene_h_list = list(mouse2human[gene].keys())
        else:
            gene_h_list = [gene]
        for gene_h in gene_h_list:
            if gene_h in hgene2repurp.keys():
                for drug in hgene2repurp[gene_h]:
                    g2repurp[gene][drug] = True

    return(g2repurp)

def getRepurpDrugs():
    """
    Returns 2 dict mappings: 
    1. human gene symbol to Broad Repurposing Hub drug.
    2. Broad Repurposing Hub drug to human gene symbol.
    """
    hgene2repurp = defaultdict(dict)
    repurp2hgene = defaultdict(dict)
    drugs = dict()
    filepath = './databases/repurposing_drugs_20200324.txt'
    df = pd.read_csv(filepath, sep='\t', engine='python',comment='!')

    for index,row in df.iterrows():
        drug = str(row['pert_iname'])
        genes = str(row['target']).split('|')
        for gene in genes:
            if gene != 'nan':
                hgene2repurp[gene][drug] = True
                repurp2hgene[drug][gene] = True
                drugs[drug] = True

    logger.info('Read %d drugs and %d targets from the Broad Drug Repurposing Hub.', len(drugs),len(hgene2repurp))
    return(hgene2repurp,repurp2hgene)

def getNetworkInteractions(args):
    """
    Compiles interactions into one file and stores it in databases dir, if needed.
    Removes unwanted PPI between a protein and itself.
    """
    logger.info('Compiling %s network interactions from data files...',args.s)
    
    #get filepath
    if args.s == 'human':
        interaction_filepath = './databases/human_interactions.txt'
    elif args.s == 'mouse':
        interaction_filepath = './databases/mouse_interactions.txt'
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    #check if interaction file exists
    if os.path.isfile(interaction_filepath):
        logger.info('Found compiled %s interactions file %s',args.s,interaction_filepath)
    else:
        #read PPI and TF
        if args.s == 'human':
            ppi = readIIDHuman()
            tf = readTF('human')
        elif args.s == 'mouse':
            ppi = readIIDMouse()
            tf = readTF('mouse')
        else:
            logger.error('%s is not a supported species.', args.s)
            sys.exit()
        #combine data and write
        ppi_frame = pd.DataFrame(ppi,columns=['source','target'])
        ppi_frame['type'] = 'PPI'
        tf_frame = pd.DataFrame(tf,columns=['source','target'])
        tf_frame['type'] = 'TF'
        frames = [ppi_frame,tf_frame]
        wholeinteract = pd.concat(frames)
        wholeinteract.to_csv(interaction_filepath,sep='\t',index=False,header=False)
        logger.info('Wrote compiled %s interactions file to %s',args.s,interaction_filepath)

    return None

def getNetworkGenes(args,verbose=False):
    """
    Reads user-defined driver and respones. Checks they are in the network.
    Identifies intermediates from interaction file.
    Returns driver, list of intermediates, and list of network respones.
    """
    if verbose:
        logger.info('Compiling network genes...')

    #Read network interactions
    if args.s == 'human':
        gene_filepath = './databases/human_genes.txt'
        interaction_filepath = './databases/human_interactions.txt'
    elif args.s == 'mouse':
        gene_filepath = './databases/mouse_genes.txt'
        interaction_filepath = './databases/mouse_interactions.txt'
    else:
        logger.error('%s is not a supported species.' % args.s)
        sys.exit()

    if not os.path.isfile(gene_filepath):

        #read interaction file
        wholeinteract = readInteractions(interaction_filepath)

        #get whole gene list
        wholegenelist = wholeinteract[0].tolist() + wholeinteract[1].tolist()
        wholegenelist = list(set(wholegenelist))
        wholegenelist.sort()

        #write to file
        wholegenedf = pd.DataFrame(wholegenelist)
        wholegenedf.to_csv(gene_filepath,sep='\t',index=False,header=False)
        logger.info('Wrote compiled %s gene list file to %s', args.s, gene_filepath)

    #check if gene list is imported
    try:
        wholegenelist
    except NameError:
        wholegenelist = readGenes(gene_filepath)
        if verbose:
            logger.info('Found compiled %s gene list file %s',args.s,gene_filepath)

    if verbose:
        logger.info('Network directory: %s',args.project_dir)
    
    #read driver and check it is in network
    driver = checkDriver(args,wholegenelist,verbose)

    #get responses in network
    responses = getNetworkResponses(args,driver,wholegenelist,verbose)

    #get intermediates in network
    intermediatesFilepath = args.project_dir + '/networkIntermediates.txt' 
    if os.path.isfile(intermediatesFilepath):
        if verbose:
            logger.info('Found network intermediates file %s. Delete this file if you want the intermediates to be compiled again',intermediatesFilepath)
    else:
        #get intermediates
        intermediates = wholegenelist

        #remove driver from intermediates
        intermediates.remove(driver)

        #remove responses from intermediates
        for r in responses:
            intermediates.remove(r)
    
        #output intermediates
        fp = open(intermediatesFilepath, 'w')
        for intermediate in intermediates:
            fp.write(intermediate + '\n')
        fp.close()
        if verbose:
            logger.info('Found %s intermediate genes in the network. Wrote to %s.' % (str(len(intermediates)),intermediatesFilepath))

    #check if intermediates is imported
    try:
        intermediates
    except NameError:
        intermediates = readGenes(intermediatesFilepath)

    #output network size
    if verbose:
        network_size = 1 + len(intermediates) + len(responses)
        logger.info('Network consists of %d genes: the driver %s, %d intermediates, and %d responses.',network_size,driver,len(intermediates),len(responses))

    return(driver,intermediates,responses)

def getPrioritization(args,driver,intermediates,responses,r=0.5,nt=16):
    """
    Calculates NetPert weights and ranks genes.
    """

    #get whole gene list
    wholegenelist = list()
    wholegenelist.append(driver)
    wholegenelist.extend(intermediates)
    wholegenelist.extend(responses)

    #get network edges 
    edge_dict = readEdges(args)
    edge_ab, edge_ba = edgeTuplesToDicts(edge_dict)

    #get genename (vertex) to index dict and vice versa
    v2i, i2v = nameToIndex(wholegenelist)

    #get genename (vertex) to indegree dict and outdegree dict
    v2indeg, v2outdeg = nameToDegree(wholegenelist,edge_ba,edge_ab)

    #get genename (vertex) to type (driver, intermediate, response)
    v2type = nameToType(wholegenelist,driver,intermediates,responses)

    # result[method] = {'weights':weights, 'ranks':ranks, 'vertices_ranked':vertices_ranked}
    results = {}

    #get laplacian matrix
    (adj_mat, lap_mat) = getLaplacian(edge_dict, v2i)

    #get time step size dt
    dt = getTimeStep(driver,v2outdeg,r=r,nt=nt)

    #get G(t) list
    g_list = getGlist(lap_mat,dt,nt)

    #get weights for NetPert
    v2weight = getWeights(args,g_list,v2i,driver,responses,nt,dt,endpoints=False)

    #get rankings for NetPert
    (vertices_ranked,v2rank) = getRankings(v2weight,ties='competition',reverse_order=True)

    results['NetPert'] = {'weights':v2weight,'ranks':v2rank,'vertices_ranked':vertices_ranked}

    return results

def getRankings(v2weight,ties='competition',reverse_order=True):
    """
    Ranks a dict of gene names to value.
    Strategies for handling ties: competition "1,2,2,4", fractional "1,2.5,2.5,4"
    Ties are ordered alphabetically by gene name.
    Returns the following:
    vertices_ranked: list of genes ranked by v2weight. Default order is descending (reverse), and
    v2rank: a dict of gene names to rank. Default order is descending (reverse).
    """
    #sort genes by weight
    pairs = sorted([ (v2weight[v], v) for v in v2weight ], reverse = reverse_order)

    vertices_ranked = list()
    v2rank = dict()

    #get ordinal rank
    (d0,v0) = pairs[0]
    previous_d = d0
    gene_list = list()
    for (d,v) in pairs:
        if isclose(d,previous_d):
            gene_list.append(v)
        else:
            gene_list.sort()
            for g in gene_list:
                vertices_ranked.append(g)
            gene_list = list()
            gene_list.append(v)
            previous_d = d
    gene_list.sort()
    for g in gene_list:
        vertices_ranked.append(g)

    ordinal_rank = 1
    for v in vertices_ranked:
        v2rank[v] = ordinal_rank
        ordinal_rank += 1

    if ties == 'competition':
        (d0,v0) = pairs[0]
        previous_d = d0
        gene_list = list()
        for (d,v) in pairs:
            if isclose(d,previous_d):
                gene_list.append(v)
            else:
                rank_list = list()
                for g in gene_list:
                    rank_list.append( v2rank[g] )
                rank = min(rank_list)
                for g in gene_list:
                    v2rank[g] = rank
                gene_list = list()
                gene_list.append(v)
                previous_d = d
        rank_list = list()
        for g in gene_list:
            rank_list.append( v2rank[g] )
        rank = min(rank_list)
        for g in gene_list:
            v2rank[g] = rank
    elif ties == 'fractional':
        (d0,v0) = pairs[0]
        previous_d = d0
        gene_list = list()
        for (d,v) in pairs:
            if isclose(d,previous_d):
                gene_list.append(v)
            else:
                ranksum = 0
                numgenes = len(gene_list)
                for g in gene_list:
                    ranksum += v2rank[g]
                for g in gene_list:
                    v2rank[g] = ranksum / numgenes
                gene_list = list()
                gene_list.append(v)
                previous_d = d
        ranksum = 0
        numgenes = len(gene_list)
        for g in gene_list:
            ranksum += v2rank[g]
        for g in gene_list:
            v2rank[g] = ranksum / numgenes

    return(vertices_ranked,v2rank)

def getWeights(args,g_list,v2i,driver,responses,nt,dt,endpoints=False):
    """
    Returns the weight for each gene in the network in a dict.
    """
    logger.info('Calculating NetPert weights for network genes...')
    v2w = dict()
    vertices = sorted(v2i.keys())

    for v in vertices:
        v2w[v] = 0.0
    d = driver
    a = v2i[d]
    for r in responses:
        b = v2i[r]
        for v in vertices:
            if not endpoints:
                if v == d or v == r:
                    continue
            x = v2i[v]
            y = [(g_list[nt-i][b,x])*(g_list[i][x,a]) for i in range(nt+1)]
            z = integrate.romb(y,dx=dt)
            v2w[v] = v2w[v] + z
    norm = 0
    for v in vertices:
        norm += v2w[v]
    if norm == 0.:
        norm = 1.
    for v in vertices:
        v2w[v] = v2w[v] / norm

    return(v2w)

def getMouse2Human():
    """
    Returns: 
    1. a dict mapping of mouse gene symbol to a dict mapping of all homolog human gene symbol. 
    2. a dict mapping of human gene symbol to a dict mapping of all homolog mouse gene symbol.
    """

    #Check if mapping file exists
    filepath_m2h = './databases/mouse_to_human_genes.tsv'
    if not os.path.isfile(filepath_m2h):
        #read data
        rawfilepath = './databases/HOM_MouseHumanSequence.rpt'
        data = pd.read_csv(rawfilepath, sep='\t')

        num2species = dict()
        m2h = dict()
        h2m = dict()

        for index,row in data.iterrows():
            num = row['DB Class Key']
            if num not in num2species.keys():
                num2species[num] = dict()
                num2species[num]['mouse'] = dict()
                num2species[num]['human'] = dict()
            species = row['Common Organism Name']
            symbol = row['Symbol']
            if species == 'mouse, laboratory':
                num2species[num]['mouse'][symbol] = True
            if species == 'human':
                num2species[num]['human'][symbol] = True

        for num in num2species.keys():
            for m in num2species[num]['mouse'].keys():
                if m not in m2h.keys():
                    if num2species[num]['human']:
                        m2h[m] = dict()
                for h in num2species[num]['human'].keys():
                    m2h[m][h] = True
                    if h not in h2m.keys():
                        h2m[h] = dict()
                    h2m[h][m] = True

        #write mouse to human mapping
        hgenes = dict()
        f = open(filepath_m2h,'w')
        f.write("mouse\thuman\n")
        for m in m2h.keys():
            for h in m2h[m]:
                f.write("%s\t%s\n" % (m,h))
                hgenes[h] = True
        f.close()
        logger.info('Mapped %d mouse genes to %d human genes. Wrote to %s.',len(m2h),len(hgenes),filepath_m2h)

    else:
        #read the existing mouse to human mapping file
        data_m2h = pd.read_csv(filepath_m2h, sep='\t')

        m2h = dict()
        h2m = dict()

        for index,row in data_m2h.iterrows():
            m = row['mouse']
            h = row['human']
            if m not in m2h.keys():
                m2h[m] = dict()
            m2h[m][h] = True
            if h not in h2m.keys():
                h2m[h] = dict()
            h2m[h][m] = True

        logger.info('Read mapping of %d mouse genes to %d human genes from %s.',len(m2h),len(h2m),filepath_m2h)

    return(m2h,h2m)

def getGlist(lap_mat,dt,nt):
    """
    Returns the G(t) as a list of matrices.
    """
    logger.info(f'Calculating the Green\'s function...')
    ret = [None] * (nt + 1)
    (nr,nc) = lap_mat.shape
    assert(nr == nc), 'Bad laplacian shape: %s' % str(lap_mat.shape)
    my_eye = np.eye(nr)
    epsilon = - dt * lap_mat
    g1 = linalg.expm(epsilon)
    curval = my_eye
    ret[0] = curval
    for i in range(nt):
        curval = np.matmul(curval, g1)
        ret[i + 1] = curval
    return(ret)

def getTimeStep(driver,v2outdeg,r=0.5,nt=16):
    """
    Returns the time step size.
    """
    # r is the relaxation parameter (the reduction in driver activity)
    # nt is the number of time steps to reach r (must be a non-negative power of 2)
    
    #calculate dt (time step size)
    den = 1 - r
    dt = -log(den,e) / (v2outdeg[driver] * nt)

    return(dt)

def getLaplacian(edge_dict,v2i):
    """
    Returns the adjacency matrix and Laplacian.
    """
    nvert = len(v2i)
    adj_mat = np.zeros( (nvert,nvert) )
    indeg_arr = np.zeros(nvert)
    outdeg_arr = np.zeros(nvert)
    lap_mat = np.zeros( (nvert,nvert) )
    for (src,dest) in edge_dict.keys():
        (a,b) = (v2i[src],v2i[dest])
        adj_mat[b,a] += 1.0
        indeg_arr[b] += 1.0
        outdeg_arr[a] += 1.0
    for a in range(nvert):
        lap_mat[a,a] = outdeg_arr[a]
    lap_mat = lap_mat - adj_mat
    return(adj_mat, lap_mat)


def nameToType(genelist,driver,intermediates,responses):
    """
    Returns a dict of gene names to type.
    D: driver.
    I: intermediate.
    R: response.
    """
    v2type = dict()
    for v in genelist:
        if v == driver:
            v2type[v] = 'D'
        elif v in intermediates:
            v2type[v] = 'I'
        elif v in responses:
            v2type[v] = 'R'
        else:
            logger.warning('Gene %s is not in the driver, intermediate, or response lists.',v)
    return(v2type)

def nameToDegree(genelist,edge_ba,edge_ab):
    """
    Returns indegree and outdegree dicts.
    Gene name is the key.
    """
    v2indeg = dict()
    v2outdeg = dict()
    for v in genelist:
        v2indeg[v] = 0
        if v in edge_ba:
            v2indeg[v] = len(edge_ba[v].keys())
        v2outdeg[v] = 0
        if v in edge_ab:
            v2outdeg[v] = len(edge_ab[v].keys())
    return(v2indeg,v2outdeg)


def nameToIndex(genelist):
    """
    Reads a list of gene names.
    Returns dicts of names to indices and vice versa
    """
    n2i = dict()
    i2n = dict()
    names = sorted(genelist)
    for (n,i) in zip( names, list(range(len(names))) ):
        n2i[n] = i
        i2n[i] = n
    return(n2i,i2n)

def edgeTuplesToDicts(edge_dict):
    """
    Reads dict of edges.
    Returns two dict-of-dicts.
    a is source, b is destination.
    ab uses sources as first key.
    ba uses dests as first key.
    """
    edge_ab = dict()
    edge_ba = dict()
    for (a,b) in edge_dict.keys():
        if a not in edge_ab:
            edge_ab[a] = dict()
        edge_ab[a][b] = True
        if b not in edge_ba:
            edge_ba[b] = dict()
        edge_ba[b][a] = True

    return(edge_ab, edge_ba)

def readEdges(args,verbose=False):
    """
    Reads interactions.
    Returns a dictionary of directed edges.
    Undirected edges are recorded as two directed edges in the dictionary.
    """
    if args.s == 'human':
        filepath = './databases/human_interactions.txt'
    elif args.s == 'mouse':
        filepath = './databases/mouse_interactions.txt'

    #check filepath exists
    if not os.path.isfile(filepath):
        logger.error('Missing databases file %s. Run NetPert analysis.',filepath)
        sys.exit()

    ret = dict()
    ndirected = 0
    nundirected = 0
    fp = open(filepath, 'r')
    for line in fp:
        my_list = [ ]
        toks = line.strip().split()
        if (len(toks) < 2) or (len(toks) >3):
            if verbose:
                logger.warning('Bad line, wrong token count: %s', line)
            continue
        (a, b) = (toks[0], toks[1])
        if len(toks)==2:
            my_list.append( (a,b) )
            ndirected += 1
        elif len(toks)==3 and toks[2]=='TF':
            my_list.append( (a,b) )
            ndirected += 1
        elif len(toks)==3 and toks[2]=='PPI':
            my_list.append( (a,b) )
            my_list.append( (b,a) )
            nundirected += 1
        else:
            if verbose:
                logger.warning('Bad line, unknown edge type: %s', line)
            continue
        for key in my_list:
            ret[key] = True

    return(ret)

def checkDriver(args,wholegenelist,verbose=False):
    """
    Reads the driver from the commandline.
    Checks if driver is in the network.
    Outputs driver as a string.
    """
    #read driver
    ret = args.driver

    #check if driver is in the network
    if ret in wholegenelist:
        if verbose:
            logger.info('Found driver %s in the network.' % ret)
    else:
        logger.error('Driver %s was not found in the network. Pick another driver.' % ret)
        sys.exit()

    return(ret)

def readDriver(args):
    """
    Reads the driver from the commandline.
    Outputs driver as a string.
    """
    #read driver
    ret = args.driver

    return(ret)

def readResponses(args,driver):
    """
    Reads a file provided by the user containing the gene symbols of the response genes.
    Removes driver from the response genes.
    Returns a list of response genes.
    """

    responseFilepath = args.project_dir + '/responses.xlsx'
    if os.path.isfile(responseFilepath):
        pass
    else:
        responseFilepath = args.project_dir + '/responses.csv'
        if os.path.isfile(responseFilepath):
            pass
        else:
            responseFilepath = args.project_dir + '/responses.tsv'
            if os.path.isfile(responseFilepath):
                pass
            else:
                logger.erro('No response file found. Must provide a file named responses with .xlsx, .csv, or .tsv extension in directory %s.',args.project_dir)
                sys.exit()

    #check if file is empty
    if not os.stat(responseFilepath).st_size:
        logger.error('The response file %s is empty. It must contain at least one gene symbol.',responseFilepath)
        sys.exit()

    #read responses
    if responseFilepath.split(".")[-1] == 'xlsx':
        xl = pd.ExcelFile(responseFilepath)
        res = len(xl.sheet_names)
        if res != 1:
            logger.error('Response file %s must contain only 1 sheet.',responseFilepath)
            sys.exit()
        responses = pd.read_excel(responseFilepath,header=None)
    elif responseFilepath.split(".")[-1] == 'csv':
        responses = pd.read_csv(responseFilepath,sep=',',header=None)
    else:
        responses = pd.read_csv(responseFilepath,sep='\t',header=None)

    responses = list(set(responses[0].tolist()))

    #remove driver from responses
    if driver in responses:
        responses.remove(driver)

    return(responses)

def getNetworkResponses(args,driver,wholegenelist,verbose=False):
    """
    Reads a file provided by the user containing the gene symbols of the response genes.
    Removes driver from the response genes.
    Checks if each response gene is in the network.
    Outputs list of network response genes called networkResponse.
    Returns a list of network responses.
    """

    #check if network response file exists
    networkResponsesFilepath = args.project_dir + '/networkResponses.txt' 
    if os.path.isfile(networkResponsesFilepath):
        if verbose:
            logger.info('Found network responses file %s. Delete this file if you want the network respones to be compiled again.',networkResponsesFilepath)
        
        #check if file is empty
        if not os.stat(networkResponsesFilepath).st_size:
            logger.error('The network responses file is empty. It must contain at least one gene symbol. Please delete this file and run NetPert again.')
            sys.exit()

        responses = pd.read_csv(networkResponsesFilepath,sep='\t',header=None)
        responses = list(set(responses[0].tolist()))
        responses.sort()

        #remove driver from network responses
        if driver in responses:
            responses.remove(driver)
            if verbose:
                logger.warning('Found driver %s in the network response file. It was removed.' % (driver))
    
        #check if responses are in network and output list
        ret = list()
        for response in responses:
            if response in wholegenelist:
                ret.append(response)
            else:
                logger.error('Not all network response genes in %s are in the network. Please delete this file and run NetPert analysis, again.', networkResponsesFilepath)
                sys.exit()

    else:
        #check if user response file exists
        responseFilepath = args.project_dir + '/responses.xlsx'
        if os.path.isfile(responseFilepath):
            pass
        else:
            responseFilepath = args.project_dir + '/responses.csv'
            if os.path.isfile(responseFilepath):
                pass
            else:
                responseFilepath = args.project_dir + '/responses.tsv'
                if os.path.isfile(responseFilepath):
                    pass
                else:
                    logger.error('No response file found. Must provide a file named responses with .xlsx, .csv, or .tsv extension in directory %s.',args.project_dir)
                    sys.exit()

        #check if file is empty
        if not os.stat(responseFilepath).st_size:
            logger.error('Response file %s is empty. Must contain at least one gene symbol.' % responseFilepath)
            sys.exit()
        
        #read responses and convert to list
        if responseFilepath.split(".")[-1] == 'xlsx':
            xl = pd.ExcelFile(responseFilepath)
            res = len(xl.sheet_names)
            if res != 1:
                logger.error('Response file %s must contain only 1 sheet.',responseFilepath)
                sys.exit()
            responses = pd.read_excel(responseFilepath,header=None)
        elif responseFilepath.split(".")[-1] == 'csv':
            responses = pd.read_csv(responseFilepath,sep=',',header=None)
        else:
            responses = pd.read_csv(responseFilepath,sep='\t',header=None)
        
        responses = list(set(responses[0].tolist()))
        responses.sort()
        if verbose:
            logger.info('Read %s response genes from %s.' % (str(len(responses)),responseFilepath))
    
        #remove driver from responses
        if driver in responses:
            responses.remove(driver)
            if verbose:
                logger.info('Found driver %s in the response list. It was removed. %s response genes remain.' % (driver,str(len(responses))))
        else:
            if verbose:
                logger.info('Did not find driver %s in the response list.' % driver)
    
        #check if responses are in network and output list
        ret = list()
        outfilepath = args.project_dir + '/networkResponses.txt' 
        fp = open(outfilepath, 'w')
        for response in responses:
            if response in wholegenelist:
                ret.append(response)
                fp.write(response + '\n')
        fp.close()
        if verbose:
            logger.info('Found %s of %s response genes in the network. Wrote to %s.' % (str(len(ret)),str(len(responses)),outfilepath))
    return(ret)

def readInteractions(filepath):
    """
    Reads a NetPert Interaction file.
    Returns a pandas dataframe
    """
    ret = pd.read_csv(filepath, sep='\t',header=None)
    return(ret)

def readGenes(filepath):
    """
    Read a NetPert Gene file.
    Returns a list
    """
    ret = pd.read_csv(filepath,sep='\t',header=None)
    ret = list(set(ret[0].tolist()))
    return(ret)

def readIIDMouse():
    """
    Reads mouse ppi file databases/mouse_annotated_PPIs.txt.
    Converts ppi from uniprot id to mgi symbol.
    Returns a list of unique ppi interactions as tuple (a,b).
    a and b are the genes of the two interacting proteins.
    Removes unwanted ppi between a protein and itself.
    """
    filename = './databases/mouse_annotated_PPIs.txt'

    #get mapping from mouse uniprot id to mgi symbol
    uniprot2mgisymbol = get_uniprot2mgisymbol()

    #read the ppi file
    #read_iid_ppi
    logger.info('Reading mouse PPI from %s. Removing PPI between a protein and itself...', filename)
    reftable = pd.read_csv(filename,sep='\t',low_memory=False)
    ret = dict()

    for index,row in reftable.iterrows():
        (a,b) = ( row['uniprot1'], row['uniprot2'] )
        if a == b:
            continue
        if a in uniprot2mgisymbol.keys() and b in uniprot2mgisymbol.keys():
            for genea in uniprot2mgisymbol[a]:
                for geneb in uniprot2mgisymbol[b]:
                    if genea == geneb:
                        continue
                    if (genea,geneb) not in ret and (geneb,genea) not in ret:
                        ret[(genea,geneb)] = True
    logger.info('Mapped %d mouse PPI interactions from uniprot id to mgi symbol.', len(ret))

    return(list(ret.keys()))

def get_uniprot2mgisymbol():
    """
    Returns a dict mapping from mouse uniprot id to mgi symbol.
    """
    #get uniprot to mgi num mapping
    uniprot2mginum = get_uniprot2mginum()
    #get mgi num to mgi symbol mapping
    mginum2symbol = get_mginum2symbol()
    ret = dict()
    counter = 0
    for prot in uniprot2mginum.keys():
        ret[prot] = dict()
        for num in uniprot2mginum[prot]:
            if num in mginum2symbol:
                ret[prot][mginum2symbol[num]] = True
                counter += 1
    logger.info('Mapped %d mouse uniprot ids to %d mgi symbols.', len(ret), counter)
    
    return(ret)

def get_mginum2symbol():
    """
    Reads databases/MRK_List2.rpt'
    Returns a dict mapping from mouse mgi num to symbol.
    """
    mginum2symbol_file = './databases/MRK_List2.rpt'
    logger.info('Reading mgi marker accession number to mgi marker symbol from %s...', mginum2symbol_file)
    reftable = pd.read_csv(mginum2symbol_file,sep='\t')
    ret = dict()
    for index,row in reftable.iterrows():
        num = row['MGI Accession ID']
        symbol = row['Marker Symbol']
        ret[num] = symbol
    logger.info('Mapped %d mgi marker accession numbers to symbols.',len(ret))
    return(ret)

def get_uniprot2mginum():
    """
    Reads databases/MOUSE_10090_idmapping.dat.
    Returns a dict mapping from mouse uniprot id to MGI number.
    """
    uniprot2mginum_file = './databases/MOUSE_10090_idmapping.dat'
    logger.info('Reading uniprot to mgi marker accession number from %s...', uniprot2mginum_file)
    ret = dict()
    ret2 = dict()
    reftable = pd.read_csv(uniprot2mginum_file, sep='\t')
    for index,row in reftable.iterrows():
        name = row.iloc[0]
        if name not in ret:
            ret[name] = dict()
        id_type = row.iloc[1]
        if id_type == 'MGI':
            mgi = row.iloc[2]
            ret[name][mgi] = True
    for prot in ret.keys():
        if ret[prot]:
            ret2[prot] = ret[prot]
    logger.info('Mapped %d of %d uniprot ids to a MGI number.', len(ret2), len(ret))
    return(ret2)

def readIIDHuman():
    """
    Reads ppi from Human IID file in databases dir.
    Returns list of unique ppi interactions as tuple (src,des).
    Removes unwanted PPI between a protein and itself.
    """
    filename = './databases/human_annotated_PPIs.txt'
    logger.info('Reading human PPI from %s. Removing PPI between a protein and itself...', filename)
    ref = pd.read_csv(filename, sep='\t',low_memory=False)

    ret = dict()

    for index,row in ref.iterrows():
        (a,b) = (row['symbol1'],row['symbol2'])
        if len(a.strip().split()) != 1:
            continue
        if len(b.strip().split()) != 1:
            continue
        if a == b:
            continue
        if (a,b) not in ret and (b,a) not in ret:
            ret[(a,b)] = True

    logger.info('Read %d unique human PPI.', len(ret))
    return(list(ret.keys()))

def readTF(species):
    """
    Reads regulator interactions from files with format: GENE1  GENE2   Mode    PMID.
    Files stored in databases dir. 
    File should not have a header.
    """
    if species == 'human': 
        filename = './databases/trrust_rawdata.human.tsv'
    elif species == 'mouse':
        filename = './databases/trrust_rawdata.mouse.tsv'
    else:
        logger.error('%s is not a supported species.',species)
        sys.exit()

    ret = dict()
    fp = open(filename,'r')
    for line in fp:
        toks = line.strip().split('\t')
        if len(toks)!=4:
            logger.warning('Bad line: wrong token count %s',line)
            continue
        (genea,geneb) = (toks[0],toks[1])
        ret[(genea,geneb)] = True
    logger.info('Read %d %s TF regulatory interactions.',len(ret),species)

    return(list(ret.keys()))

def main():
    args = CMDParser()
    args.func(args)
    return None

if __name__=='__main__':
    main()
