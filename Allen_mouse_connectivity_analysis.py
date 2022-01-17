from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
import webbrowser
mcc = MouseConnectivityCache()

# open up a list of all of the experiments
#all_experiments = mcc.get_experiments(dataframe=True) # a dataframe
#structure_tree = mcc.get_structure_tree()
#visp = structure_tree.get_structures_by_name(['Primary visual area'])[0] #[0] if there is more than one structure referenced, but since there is only VISp here the [0] could be removed

def experiment_page (id_experiment) :
    link = "https://connectivity.brain-map.org/projection/experiment/" + str(id_experiment)
    webbrowser.open(link)
experiment_page(638314843)

def experiments_structure (name_list) :
    '''
    Returns Dataframes containing for each structure all its experiments and its info ('gender', 'injection_structures', 'injection_volume', 'injection_x',
       'injection_y', 'injection_z', 'product_id', 'specimen_name', 'strain',
       'structure_abbrev', 'structure_id', 'structure_name', 'transgenic_line',
       'transgenic_line_id', 'id', 'primary_injection_structure')

    Parameters
    ----------
    name_list : list
        List of the structure names we want to take into account (e.g. ['Primary visual area','Hypothalamus'])

    Returns
    -------
    liste: list
        List of Dataframes
    '''

    structure_tree = mcc.get_structure_tree()
    all_structures=structure_tree.get_structures_by_name(name_list)
    liste=[]
    for i in all_structures :
        experiments = mcc.get_experiments(dataframe=True, injection_structure_ids=[i['id']]) #injection_structure_ids is here to take into account only the ids that are in the given structure
        liste.append(experiments)
    return liste
s=experiments_structure(['Primary visual area']) #to read a line we have to put within the loc the id number instead of the number of the line
#print(s[0].loc[156545918])


def cre_experiments (name, cre_type) :
    '''
    Returns Dataframes containing for a given structure and cre line all its experiments and its info ('gender', 'injection_structures', 'injection_volume', 'injection_x',
       'injection_y', 'injection_z', 'product_id', 'specimen_name', 'strain',
       'structure_abbrev', 'structure_id', 'structure_name', 'transgenic_line',
       'transgenic_line_id', 'id', 'primary_injection_structure')

    Parameters
    ----------
    name : str
        name of the structure (e.g. 'Primary visual area','Hypothalamus')

    cre_type : str
        name of the cre-line

    Returns
    -------
    cre_experiments : Dataframe
        all experiments info for a given cre line and structure
    '''

    structure_tree = mcc.get_structure_tree()
    expe = structure_tree.get_structures_by_name([name])[0]
    cre_experiments = mcc.get_experiments(dataframe=True, cre=cre_type, injection_structure_ids=[expe['id']])  # returns only experiments cre line postives
    return cre_experiments
j=cre_experiments('Primary visual area',['Cux2-IRES-Cre']) #(e.g. 'Rbp4-Cre_KL100')
#print(j.loc[263780729])


def stru_tree_name (name) :
    '''
    Returns general info of the given structure and its list of structure_set_ids, i.e. ids of lists giving specific and diverse info in wichi the given structure is related

    Parameters
    ----------
    name : str
        name of the structure (e.g. 'Primary visual area','Hypothalamus')

    Returns
    -------
    ids : list
        List of the structure_set_ids (each id correspond to a list in which the given structure is part of.

    structure_info : Dataframe
        Datafrme of the given structure, containing general info of it ('acronym', 'graph_id', 'graph_order', 'id', 'name', 'structure_id_path',
       'structure_set_ids', 'rgb_triplet')

    '''
    structure_tree = mcc.get_structure_tree()
    structures = structure_tree.get_structures_by_name([name])
    a=pd.DataFrame(structures)
    ids = a.loc[0]['structure_set_ids']
    structure_info=a.loc[0]
    return (ids,structure_info)
r,a=stru_tree_name('Primary visual area')



def structure_set_id_info (ids) : # for the structure_set_ids of a given structure , we match the properties of those ids , can take one id, but it has to be in list type
    '''
    Returns for each id its associated information list

    Parameters
    ----------
    ids : list
        list of ids corresponding to the ones associated to information lists

    Returns
    -------
    liste : list
        List of all the information lists (containing their name and name description) associated to the given ids

    '''

    oapi = OntologiesApi()
    liste=[]
    structure_tree = mcc.get_structure_tree()
    structure_set_ids = structure_tree.get_structure_sets()
    b = pd.DataFrame(oapi.get_structure_sets(structure_set_ids)) #b gives all the information lists available
    #print(b)
    for i in ids :
        index=b[b['id']==i].index[0]
        liste.append(b.loc[index])
    return (list(liste))
d=structure_set_id_info(r)
o=[396673091]
#print(structure_set_id_info(o))


def stru_tree_id (list_id) : #for a given id refering to a list of info, it returns all the info of this id/list
    '''
    Returns all the structures (and their info) that are contained in the given list (which is referenced here by its id)

    Parameters
    ----------
    list_id : str
        id of the information list

    Returns
    -------
    stru_list : Dataframe
        Dataframe of all the structures contained in the given information list, and general information of those structures ('acronym', 'graph_id', 'graph_order', 'id', 'name', 'structure_id_path',
       'structure_set_ids', 'rgb_triplet')

    '''
    structure_tree = mcc.get_structure_tree()
    summary_structures = structure_tree.get_structures_by_set_id([list_id])
    stru_list = pd.DataFrame(summary_structures)
    return stru_list
q=stru_tree_id(396673091)
#print(q.loc[q[q['acronym']=='VISp2/3'].index[0]]) #gives the line of VISp2/3 by refering to the row number

structure_tree = mcc.get_structure_tree()
isocortex = structure_tree.get_structures_by_name(['Isocortex'])[0]
# find wild-type injections into primary visual area
visp = structure_tree.get_structures_by_acronym(['VISp'])[0]
visp_experiments = mcc.get_experiments(cre=False,
                                       injection_structure_ids=[visp['id']])

iso=structure_tree.get_structures_by_acronym(['SSp-ll1'])[0]
structure_unionizes = mcc.get_structure_unionizes([ e['id'] for e in visp_experiments ], is_injection=False,structure_ids=[isocortex['id']],include_descendants=True)
stru_uni=mcc.get_structure_unionizes([ e['id'] for e in visp_experiments ], is_injection=False,structure_ids=[iso['id']],include_descendants=True)
print("%d VISp non-injection, cortical structure unionizes" % len(structure_unionizes))
print(structure_unionizes.loc[0])

print("%d VISp non-injection, cortical structure unionizes" % len(stru_uni))
print(stru_uni.loc[96])
ee=structure_tree.get_structures_by_id([1030])
print(ee)
all_experiments = mcc.get_experiments(dataframe=True)
#print(all_experiments.loc[307297141])

#307297141 : experiment id

#with unionize we look for each experiment in a given structure (by giving their id experiments) it gives the the density of projecting signal, volume of projecting signal (and where it was mesured ? meaning structure_id)
#may be what we're looking for : from a structure in which the experiment is done, it gives the projection volume in the other structures

#create function : for a given id/structure injection of the experiment, we also give a structure acronym to see if there is projection in this structure and at what proportion/volume/density