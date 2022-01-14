from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi

mcc = MouseConnectivityCache()

# open up a list of all of the experiments
#all_experiments = mcc.get_experiments(dataframe=True) # a dataframe
#structure_tree = mcc.get_structure_tree()
#visp = structure_tree.get_structures_by_name(['Primary visual area'])[0] #[0] if there is more than one structure referenced, but since there is only VISp here the [0] could be removed


def experiments_structure (name_list) : # returns list of tables of all the experiments and their info for the given structures (or a list of structures, in that cas the tables are returns one after an other
    structure_tree = mcc.get_structure_tree()
    all_structures=structure_tree.get_structures_by_name(name_list)
    liste=[]
    for i in all_structures :
        experiments = mcc.get_experiments(dataframe=True, injection_structure_ids=[i['id']])
        liste.append(experiments)
    return liste
s=experiments_structure(['Primary visual area']) #to read a line we have to put within the loc the id number instead of the number of the line

def cre_experiments (name, cre_type) : #returns for a given strcture the experiments having or having not a cre line (can also refers the cre line type (e.g. 'Rbp4-Cre_KL100')
    structure_tree = mcc.get_structure_tree()
    expe = structure_tree.get_structures_by_name([name])[0]
    cre_experiments = mcc.get_experiments(dataframe=True, cre=cre_type, injection_structure_ids=[expe['id']])  # returns only experiments cre line postives
    return cre_experiments
j=cre_experiments('Primary visual area',['Cux2-IRES-Cre'])
#print(j)


def stru_tree_name (name) : #for a given structure it returns the list of structure_set_ids and a short table summarizing the given structure
    # to retrieve the mouse structure tree : list of dictionaries where each dictionary describes a structure.
    structure_tree = mcc.get_structure_tree()
    structures = structure_tree.get_structures_by_name([name])
    a=pd.DataFrame(structures)
    ids = a.loc[0]['structure_set_ids']
    return (ids,a.loc[0])
r,a=stru_tree_name('Primary visual area')


def structure_set_id_info (ids) : # for the structure_set_ids of a given structure , we match the properties of those ids , can take one id, but it has to be in list type
    oapi = OntologiesApi()
    liste=[]
    structure_tree = mcc.get_structure_tree()
    # get the ids of all the structure sets in the tree
    structure_set_ids = structure_tree.get_structure_sets()
    # query the API for information on those structure sets
    b = pd.DataFrame(oapi.get_structure_sets(structure_set_ids))
    for i in ids :
        index=b[b['id']==i].index[0]
        liste.append(b.loc[index])
    return (list(liste))
d=structure_set_id_info(r)
o=[396673091]
#print(structure_set_id_info(o))



def stru_tree_id (stru_id) : #for a given id refering to a list of info, it returns all the info of this id/list
    structure_tree = mcc.get_structure_tree()
    summary_structures = structure_tree.get_structures_by_set_id([stru_id])
    c = pd.DataFrame(summary_structures)
    return c
q=stru_tree_id(396673091)

#print(q.loc[q[q['acronym']=='VISp2/3'].index[0]]) #gives the line of VISp2/3 by refering to the row number





