from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
import webbrowser
import numpy as np
import matplotlib.pyplot as plt


mcc = MouseConnectivityCache()

# open up a list of all of the experiments
#all_experiments = mcc.get_experiments(dataframe=True) # a dataframe
#structure_tree = mcc.get_structure_tree()
#visp = structure_tree.get_structures_by_name(['Primary visual area'])[0] #[0] if there is more than one structure referenced, but since there is only VISp here the [0] could be removed

def experiment_page (experiment_id) :
    '''
    The web page informing about the experiment and the projection volume/density in the Allen Brain Atlas is opened
    Parameters
    ----------
    experiment_id : int
    '''

    link = "https://connectivity.brain-map.org/projection/experiment/" + str(experiment_id)
    webbrowser.open(link)
#experiment_page(113887162)

def injection_cortical_map_page (experiment_id):
    '''
    The web page informing about the experiment and the projection volume/density in the Allen Brain Atlas is opened
    Parameters
    ----------
    experiment_id : int
    '''

    link= "https://connectivity.brain-map.org/projection/experiment/cortical_map/"+str(experiment_id)
    webbrowser.open(link)
#injection_cortical_map_page(113887162)

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
q=stru_tree_id(114512892)
#print(q.loc[q[q['acronym']=='VISp2/3'].index[0]]) #gives the line of VISp2/3 by refering to the row number


def all_structures (): #give all the structures experiments have been done on
    all_experiments = mcc.get_experiments(dataframe=True)
    ids = list(all_experiments['id'])
    structure_dict = dict()
    for i in ids:
        info = all_experiments.loc[i]
        if info['structure_name'] not in structure_dict:
            structure_dict[info['structure_name']] = 1
        else:
            structure_dict[info['structure_name']] = structure_dict[info['structure_name']] + 1
    return structure_dict

a=all_structures()
#print(a)

def plot_x_structures (number) : #plot the first x structures with the most experiments (e.g the structures having the most injection experiments)
    all_experiments = mcc.get_experiments(dataframe=True)
    ids = list(all_experiments['id'])
    structure_dict = dict()
    for i in ids:
        info = all_experiments.loc[i]
        if info['structure_name'] not in structure_dict:
            structure_dict[info['structure_name']] = 1
        else:
            structure_dict[info['structure_name']] = structure_dict[info['structure_name']] + 1
    values=structure_dict.values()
    keys=structure_dict.keys()
    maxi=dict()
    i=0
    while i<number:
        maximum=max(structure_dict,key=structure_dict.get)
        maxi[maximum]=structure_dict[maximum]
        structure_dict.pop(maximum)
        i=i+1
    print(maxi)
    title='the first '+str(number)+' structures with the highest number of experiments'
    plt.bar(maxi.keys(), maxi.values(), width=0.5,align='center')
    plt.xticks(rotation=10)
    plt.title(title)
    plt.show()
#plot_x_structures(5)


def structure_projection (name_structure_injection, experiment_id, cre) : #cre False or True ; 26min for first run and 1min50 for second run
    #for a given id, it returns all the sub structures (takes into account the layers) it projects to
    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_structure_injection])[0]
    structure_injection_experiments = mcc.get_experiments(cre=cre,injection_structure_ids=[structure_injection['id']])
    print('if first run of this function in the script, it may take some time to compute')
    unionize=mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False, include_descendants=True) #long to compute at first
    row_numbers= np.arange(0,len(unionize))
    j=0
    row_list=[]
    for i in row_numbers :
        j=j+1
        if unionize.loc[i]["experiment_id"]==experiment_id :
            row_list.append(j)
    new_dataframe=unionize[row_list[0]:row_list[-1]+1]
    new_row_numbers=np.arange(row_list[0],row_list[-1])
    projection_structures_dict = dict()
    for i in new_row_numbers :
        info=structure_tree.get_structures_by_id([new_dataframe.loc[i]['structure_id']])
        if info[0]['name'] not in projection_structures_dict :
            projection_structures_dict[info[0]['name']] = 1
        else:
            projection_structures_dict[info[0]['name']] = projection_structures_dict[info[0]['name']] + 1
    projection_structures_name=list(projection_structures_dict.keys())
    return (projection_structures_name)

#h=structure_projection('Primary visual area',638314843,None)
#print(h)

#function which takes into account an injection_structure and a projection_structure (but not the id), if target not specified it returns for all the experiment (x3 for the hemisphere) the projections results for each target_structure
def injection_projection_structure (name_structure_injection,cre,target_projection) :
    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_structure_injection])[0]
    structure_injection_experiments = mcc.get_experiments(cre=cre, injection_structure_ids=[structure_injection['id']])
    if target_projection!=None :
        target=structure_tree.get_structures_by_name([target_projection])[0]
        unionize = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False,structure_ids=[target['id']],include_descendants=True)
        print('There are ' + str(len(unionize)/3) + ' experiments in which ' +str(name_structure_injection)+' projects to ' +str(target_projection))
    else :
        unionize = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False, include_descendants=True)
        print("the projection results of the "+str(len(unionize)/282)+" experiments with an injection done within "+str(name_structure_injection))
    return (unionize)

#k=injection_projection_structure('Primary visual area',None,None)
#print(k)


#function which returns, for a given id (and injection_structure) and projection_structure a table with volume_porjection (as well as injection volume) for each hemisphere


def id_projection_structure (name_structure_injection,cre,target_projection,experiment_id) :
    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_structure_injection])[0]
    structure_injection_experiments = mcc.get_experiments(cre=cre, injection_structure_ids=[structure_injection['id']])
    target = structure_tree.get_structures_by_name([target_projection])[0]
    unionize = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False,structure_ids=[target['id']], include_descendants=True)
    row_numbers = np.arange(0, len(unionize))
    row_list = []
    for i in row_numbers:
        if unionize.loc[i]["experiment_id"] == experiment_id:
            row_list.append(i)
    print(row_list)
    print(row_list[-1])
    new_dataframe = unionize[row_list[0]:row_list[-1]+1]
    return (new_dataframe)

dd=id_projection_structure('Primary visual area',None,'Primary visual area, layer 6b',263780729)
print(dd)

#continue the function by creating for a given id a dataframe with in column each hemi (left, right and both) and in row the parameter

#500836840 ; 307297141 ; 503069254 ; 263780729
