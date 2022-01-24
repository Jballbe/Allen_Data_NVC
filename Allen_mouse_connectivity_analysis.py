from allensdk.core.mouse_connectivity_cache import MouseConnectivityCache
import pandas as pd
from allensdk.api.queries.ontologies_api import OntologiesApi
import webbrowser
import numpy as np
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')



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
#experiment_page(263780729)

def cortical_map_page (experiment_id):
    '''
    The web page informing showing the projections of the experiment in a cortical map

    Parameters
    ----------
    experiment_id : int
    '''

    link= "https://connectivity.brain-map.org/projection/experiment/cortical_map/"+str(experiment_id)
    webbrowser.open(link)
#cortical_map_page(263780729)

def connectivity_map (experiment_id) :
    '''
    The web page shows in a 3D volume de connectivity (projections) of an experiment

    Parameters
    ----------
    experiment_id : int
    '''
    link = "https://connectivity.brain-map.org/3d-viewer?v=1&types=STREAMLINE&STREAMLINE=" + str(experiment_id)
    webbrowser.open(link)
#connectivity_map(263780729)

def experiment_info (experiment_id, structure_name) :
    '''
    General info of the given experiment (associated to a specific structure)

    Parameters
    ----------
    experiment_id : int

    structure_name : str
        should be the structure related within the experiment

    Returns
    ----------
    id_info : Dataframe
        Table of the general info of the experiment


    '''

    structure_tree = mcc.get_structure_tree()
    all_structures = structure_tree.get_structures_by_name([structure_name])[0]
    experiments = mcc.get_experiments(dataframe=True, injection_structure_ids=[all_structures['id']])
    c=experiments.loc[experiment_id]
    id_info=pd.DataFrame(c)
    return(id_info)

#print(experiment_info(263780729,'Primary visual area'))

def experiments_structure (structures) :
    '''
    Returns Dataframes containing for each structure all its experiments and its info ('gender', 'injection_structures', 'injection_volume', 'injection_x',
       'injection_y', 'injection_z', 'product_id', 'specimen_name', 'strain',
       'structure_abbrev', 'structure_id', 'structure_name', 'transgenic_line',
       'transgenic_line_id', 'id', 'primary_injection_structure')

    Parameters
    ----------
    structures : list
        List of the structure names we want to take into account (e.g. ['Primary visual area','Hypothalamus'])

    Returns
    -------
    liste: list
        List of Dataframes
    '''

    structure_tree = mcc.get_structure_tree()
    all_structures=structure_tree.get_structures_by_name(structures)
    liste=[]
    for i in all_structures :
        experiments = mcc.get_experiments(dataframe=True, injection_structure_ids=[i['id']]) #injection_structure_ids is here to take into account only the ids that are in the given structure
        liste.append(experiments)
    return liste
#s=experiments_structure(['Primary visual area']) #to read a line we have to put within the loc the id number instead of the number of the line
#print(s[0].loc[156545918])

def cre_lines (cre_line_name) :
    '''
    Returns all the cre lines and the associated experiment ids or, if the name of a specific cre line is referenced, it returns its associated experiment ids

    Parameters
    ----------
    cre_line_name : str
        name of a the cre line

    Returns
    -------
    dico_id : Dictionary
        Dictionary of all the cre lines and their experiment ids

    ids : List
        List of the experiment ids using the specific cre line mentioned

    '''

    all_experiments = mcc.get_experiments(dataframe=True)  # a dataframe
    dico_id = dict()
    dico_number = dict()
    for i in all_experiments.index:
        if all_experiments.loc[i]['transgenic_line'] not in dico_number:
            dico_number[all_experiments.loc[i]['transgenic_line']] = 1
            dico_id[all_experiments.loc[i]['transgenic_line']] = []
            dico_id[all_experiments.loc[i]['transgenic_line']].append(i)

        else:
            dico_number[all_experiments.loc[i]['transgenic_line']] = dico_number[ all_experiments.loc[i]['transgenic_line']] + 1
            dico_id[all_experiments.loc[i]['transgenic_line']].append(i)
    if cre_line_name == None :
        return (dico_id)
    else :
        ids=dico_id[cre_line_name]
        return(ids)

def cre_experiments (structure, cre_type) :
    '''
    Returns Dataframes containing for a given structure and cre line all its experiments and its info ('gender', 'injection_structures', 'injection_volume', 'injection_x',
       'injection_y', 'injection_z', 'product_id', 'specimen_name', 'strain',
       'structure_abbrev', 'structure_id', 'structure_name', 'transgenic_line',
       'transgenic_line_id', 'id', 'primary_injection_structure')

    Parameters
    ----------
    structure : str
        name of the structure (e.g. 'Primary visual area','Hypothalamus')

    cre_type : str, True or False or None
        name of the cre-line. If True it takes into account all the experiments having a cre-line, and False takes into account the one which don't have a cre-line. None takes all the experiments having a cre line or not.

    Returns
    -------
    cre_experiments : Dataframe
        all experiments info for a given cre line and structure
    '''

    structure_tree = mcc.get_structure_tree()
    expe = structure_tree.get_structures_by_name([structure])[0]
    cre_experiments = mcc.get_experiments(dataframe=True, cre=cre_type, injection_structure_ids=[expe['id']])  # returns only experiments cre line postives
    return cre_experiments
j=cre_experiments('Primary visual area',['Cux2-IRES-Cre']) #(e.g. 'Rbp4-Cre_KL100')
#print(j)
#print(j.loc[263780729])


def structure_tree (structure) :
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
    structures = structure_tree.get_structures_by_name([structure])
    a=pd.DataFrame(structures)
    ids = a.loc[0]['structure_set_ids']
    structure_info=a.loc[0]
    return (ids,structure_info)
r,a=structure_tree('Primary visual area')




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
#d=structure_set_id_info(r)

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
#q=stru_tree_id(114512892)

#print(q.loc[q[q['acronym']=='VISp2/3'].index[0]]) #gives the line of VISp2/3 by refering to the row number


def all_structures (): #give all the structures experiments have been done on
    '''
    Returns all the structures that have received injections and their associated experiment ids

    Parameters
    ----------
    None

    Returns
    -------
    structure_dict : Dictionary
        Dictionary of all the structures with the number (as values) of experiments that have been done on them

    structure_id_dict : Dictionary
        Dictionary of all the structures (as keys) with the list of experiment ids (as values) corresponding to the experiments that have been done on those structures

    '''

    all_experiments = mcc.get_experiments(dataframe=True)
    ids = list(all_experiments['id'])
    structure_dict = dict()
    structure_id_dict= dict()
    for i in ids:
        info = all_experiments.loc[i]
        if info['structure_name'] not in structure_dict:
            structure_dict[info['structure_name']] = 1
            structure_id_dict[info['structure_name']] = []
            structure_id_dict[info['structure_name']].append(i)
        else:
            structure_dict[info['structure_name']] = structure_dict[info['structure_name']] + 1
            structure_id_dict[info['structure_name']].append(i)
    return (structure_dict, structure_id_dict)

#a,b=all_structures()



def plot_x_structures (number) : #plot the first x structures with the most experiments (e.g the structures having the most injection experiments)
    '''
    Returns an histogram of the structures having the highest number of experiments done in them. The number of structure plotted can be choosen

    Parameters
    ----------
    number : int or None
        Corresponds to the number of structures we want to be plotted

    Returns
    -------
    plot

    '''

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
    if number!=None :
        while i < number:
            maximum = max(structure_dict, key=structure_dict.get)
            maxi[maximum] = structure_dict[maximum]
            structure_dict.pop(maximum)
            i = i + 1
        print(maxi)
        title = 'the first ' + str(number) + ' structures with the highest number of experiments'
        plt.bar(maxi.keys(), maxi.values(), width=0.5, align='center')
        plt.xticks(rotation=10)
        plt.title(title)
        plt.show()
    else :
        title = 'all the structures according to their number of experiments'
        plt.bar(keys, values, width=0.5, align='center')
        plt.xticks(rotation=10)
        plt.title(title)
        plt.show()

#plot_x_structures(5)


def structure_projection (name_injection_structure, experiment_id, cre) : #cre False or True ; 26min for first run and 1min50 for second run
    #add the is_injection parameter
    '''
    For a given experiment it returns all the targeted structures the structure receiving the injection sends its projections to

    Parameters
    ----------
    name_injection_structure : str
        the structure receiving the injection

    experiment_id : int

    cre : str or False/None
        Name of the cre-line (should be None ? since we already have a unique experiment chosen)

    Returns
    -------
    projection_structures_name : List
        List of the (sub)structures (takes into account the layers of the structures) receiving the projections in the given experiment (and from the given structure receiving the injection)

    '''

    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_injection_structure])[0]
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

#h=structure_projection('Primary visual area',263780729,None)
#print(h)
#print(len(h))


def injection_projection_structure (name_injection_structure,cre,target_projection) :
    '''
    For a given structure receiving the injection and the targeted structure receiving the projections from this injection structure, it returns all the experiment ids and their projection info

    Parameters
    ----------
    name_injection_structure : str
        the structure receiving the injection

    cre : str or False/None
        Name of the cre-line (should be None ? since we already have a unique experiment chosen)

    target_projection : str or None
        the targeted structure receiving the projections. No structures can be referenced, in that case we should put None.

    Returns
    -------
    unionize : Dataframe
        Dataframe of all the experiments (and their projection info) having the given structure receiving the injection, and receiving projections in the given targeted structure. Each of those experiments are
        returned 3 times (or less) because the projections are mesured in the left and right hemsiphere, as well as in both at the same time. If target_projection is None, it returns a Dataframe of all the
        experiments sending projections to all the different targeted structures that the given structure receiving the injection can send to (not a specific structure)

    '''

    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_injection_structure])[0]
    structure_injection_experiments = mcc.get_experiments(cre=cre, injection_structure_ids=[structure_injection['id']])
    if target_projection!=None :
        target=structure_tree.get_structures_by_name([target_projection])[0]
        unionize = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False,structure_ids=[target['id']],include_descendants=True)
        print('There are ' + str(len(unionize)/3) + ' experiments in which ' +str(name_injection_structure)+' projects to ' +str(target_projection))
    else :
        unionize = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False, include_descendants=True)
        print("the projection results of the "+str(len(unionize)/282)+" experiments with an injection done within "+str(name_injection_structure))
    return (unionize)

#k=injection_projection_structure('Primary visual area',None,'Ectorhinal area/Layer 6a')
#print(k)
#in order to access the projection values of those targeted structure for the given experiment id, we can use then the id_projection_structure function

def density_volume_minimum (experiment_id,name_injection_structure,parameters,minimal_values) :
    '''
    For a given experiment id and minimal values for some projection parameters it returns a table with the different targeted structures which respond to this selection

    Parameters
    ----------
    experiment_id : int

    name_injection_structure : str

    parameters : List
        List of parameters that will have a minimal value (in order to take into account only the targeted structures having those parameters at an higher value). Can take : "projection_volume", "projection_density",
        "volume" and "normalized_projection_volume". parameters can take 1 or 2 parameters.

    minimal_values : List
        List of values in float, corresponds to the minimal values of the concerned parameters the targeted structures can have in order to be returned

    Returns
    -------
    info : Dataframe
        Dataframe of the targeted structures which have parameters at a certain value

    '''

    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_injection_structure])[0]
    structure_injection_experiments = mcc.get_experiments(injection_structure_ids=[structure_injection['id']])
    unionizee = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False,
                                           include_descendants=True)
    u=unionizee.set_index('experiment_id')
    unionize=u.loc[experiment_id]
    if len(parameters) == 1 :
        if parameters[0] == 'projection_volume' :
            proj_volume=unionize[unionize.projection_volume > minimal_values[0]]
            info=pd.DataFrame(structure_tree.nodes(proj_volume.structure_id))
            return(info)
        else :
            if parameters[0] == 'projection_density' :
                proj_density = unionize[unionize.projection_density > minimal_values[0]]
                info = pd.DataFrame(structure_tree.nodes(proj_density.structure_id))
                return (info)
            else :
                if parameters[0] == "volume" :
                    vol=unionize[unionize.volume > minimal_values[0]]
                    info = pd.DataFrame(structure_tree.nodes(vol.structure_id))
                    return (info)
                else :
                    norm = unionize[unionize.normalized_projection_volume>minimal_values[0]]
                    info = pd.DataFrame(structure_tree.nodes(norm.structure_id))
                    return (info)
    else :
        if len(parameters) == 2 :
            if parameters[0] == 'projection_volume':
                proj_volume = unionize[unionize.projection_volume > minimal_values[0]]
                if parameters[1] == 'projection_density' :
                    proj_vol_dens = proj_volume[proj_volume.projection_density >minimal_values[1]]
                    info = pd.DataFrame(structure_tree.nodes(proj_vol_dens.structure_id))
                    return(info)
                else :
                    if parameters[1] == "volume" :
                        vol_vol = proj_volume[proj_volume.volume > minimal_values[1]]
                        info = pd.DataFrame(structure_tree.nodes(vol_vol.structure_id))
                        return (info)
                    else :
                        proj_vol_norm = proj_volume[proj_volume.normalized_projection_volume > minimal_values[1]]
                        info = pd.DataFrame(structure_tree.nodes(proj_vol_norm.structure_id))
                        return (info)
            else :
                if parameters[0] == 'projection_density' :
                    proj_density = unionize[unionize.projection_density > minimal_values[0]]
                    if parameters[1] == 'projection_volume' :
                        proj_dens_vol = proj_density[proj_density.projection_volume > minimal_values[1]]
                        info = pd.DataFrame(structure_tree.nodes(proj_dens_vol.structure_id))
                        return (info)
                    else :
                        if parameters[1] == "volume" :
                            proj_dens_vol = proj_density[proj_density.volume > minimal_values[1]]
                            info = pd.DataFrame(structure_tree.nodes(proj_dens_vol.structure_id))
                            return (info)
                        else :
                            proj_dens_norm=proj_density[proj_density.normalized_projection_volume > minimal_values[1]]
                            info = pd.DataFrame(structure_tree.nodes(proj_dens_norm.structure_id))
                            return (info)
                else :
                    if parameters[0] == 'volume' :
                        vol=unionize[unionize.volume > minimal_values[0]]
                        if parameters[1] == 'projection_volume' :
                            vol_proj_vol = vol[vol.projection_volume>minimal_values[1]]
                            info=pd.DataFrame(structure_tree.nodes(vol_proj_vol.structure_id))
                            return (info)
                        else:
                            if parameters[1] == "projection_density" :
                                vol_dens=vol[vol.projection_density>minimal_values[1]]
                                info = pd.DataFrame(structure_tree.nodes(vol_dens.structure_id))
                                return(info)
                            else :
                                vol_norm=vol[vol.normalized_projection_volume>minimal_values[1]]
                                info = pd.DataFrame(structure_tree.nodes(vol_norm.structure_id))
                                return (info)
                    else :
                        if parameters[0] == "normalized_projection_volume" :
                            normalized_vol = unionize[unionize.normalized_projection_volume> minimal_values[0]]
                            if parameters[1] == "projection_volume" :
                                norm_proj_vol=normalized_vol[normalized_vol.projection_volume>minimal_values[1]]
                                info = pd.DataFrame(structure_tree.nodes(norm_proj_vol.structure_id))
                                return (info)
                            else :
                                if parameters[1] == "projection_density" :
                                    norm_proj_dens = normalized_vol[normalized_vol.projection_density>minimal_values[1]]
                                    info = pd.DataFrame(structure_tree.nodes(norm_proj_dens.structure_id))
                                    return (info)
                                else :
                                    norm_vol = normalized_vol[normalized_vol.volume > minimal_values[1]]
                                    info = pd.DataFrame(structure_tree.nodes(norm_vol.structure_id))
                                    return (info)
        else :
            print ("too many info in parameters, this list can only takes 1 or 2 parameters at the time")
            return ([])

#y=density_volume_minimum(263780729,'Primary visual area',['normalized_projection_volume',"projection_density"],[0.05,0.05]).to_string()




def id_projection_structure (name_injection_structure,cre,target_projection,experiment_id) :
    '''
    For a given experiment id, structure receiving the injection, targeted structure receiving the projections and the cre-line, it returns the projection info of this experiment id sending its projection to the given
    targeted structure.

    Parameters
    ----------
    name_injection_structure : str
        the structure receiving the injection

    cre : str or None
        Name of the cre-line (should be None ? since we already have a unique experiment chosen)

    target_projection : str or None
        the targeted structure receiving the projections. No structures can be referenced, in that case we should put None.

    experiment_id : int

    Returns
    -------
    id_projection_info : Dataframe
        Dataframe of the given experiment id and its projection info. It may returns 3 lines corresponding to the projection measurements in each hemisphere and in the two hemisphere at the same time. Without injection info

    id_projection_injection_info : Dataframe
        Dataframe of the given experiment id and its projection info. It may returns 3 lines corresponding to the projection measurements in each hemisphere and in the two hemisphere at the same time. With injection info.

    '''

    structure_tree = mcc.get_structure_tree()
    structure_injection = structure_tree.get_structures_by_name([name_injection_structure])[0]
    structure_injection_experiments = mcc.get_experiments(cre=cre, injection_structure_ids=[structure_injection['id']])
    target = structure_tree.get_structures_by_name([target_projection])[0]
    unionize = mcc.get_structure_unionizes([e['id'] for e in structure_injection_experiments], is_injection=False,structure_ids=[target['id']], include_descendants=True)
    row_numbers = np.arange(0, len(unionize))
    row_list = []
    for i in row_numbers:
        if unionize.loc[i]["experiment_id"] == experiment_id:
            row_list.append(i)
    new_dataframe = unionize[row_list[0]:row_list[-1]+1]
    a = list(new_dataframe.index)
    to_change_id = []
    to_change_hemi = []
    for i in a:
        if new_dataframe.loc[i]['hemisphere_id'] == 1:
            to_change_id.append(i)
            to_change_hemi.append('left hemisphere')
        else:
            if new_dataframe.loc[i]['hemisphere_id'] == 2:
                to_change_id.append(i)
                to_change_hemi.append('right hemisphere')
            else:
                to_change_id.append(i)
                to_change_hemi.append('both hemispheres')
    if len(to_change_id) == 3:
        new = new_dataframe.rename(index={to_change_id[0]: to_change_hemi[0], to_change_id[1]: to_change_hemi[1],
                                  to_change_id[2]: to_change_hemi[2]})
        a = experiment_info(experiment_id, name_injection_structure)
        values = list(a[experiment_id].loc["injection_volume":'injection_z'])
        rows = ["injection_volume", "injection_x", "injection_y", "injection_z"]
        col = [to_change_hemi[0], to_change_hemi[1], to_change_hemi[2]]
        v = [values, values, values]
        new_mini = pd.DataFrame(v, columns=rows, index=col)
    else:
        if len(to_change_id) == 2:
            new = new_dataframe.rename(index={to_change_id[0]: to_change_hemi[0], to_change_id[1]: to_change_hemi[1]})
            a = experiment_info(experiment_id, name_injection_structure)
            values = list(a[experiment_id].loc["injection_volume":'injection_z'])
            rows = ["injection_volume", "injection_x", "injection_y", "injection_z"]
            col = [to_change_hemi[0], to_change_hemi[1]]
            v = [values, values]
            new_mini = pd.DataFrame(v, columns=rows, index=col)
        else:
            new = new_dataframe.rename(index={to_change_id[0]: to_change_hemi[0]})
            a = experiment_info(experiment_id, name_injection_structure)
            values = list(a[experiment_id].loc["injection_volume":'injection_z'])
            rows = ["injection_volume", "injection_x", "injection_y", "injection_z"]
            col = [to_change_hemi[0]]
            v = [values]
            new_mini = pd.DataFrame(v, columns=rows, index=col)
    #print('Experiment id ' + str(experiment_id) + ' projection info:' + '\n', new.transpose())
    id_projection_info=new.transpose() #without injection info
    id_projection_injection_info=id_projection_info.append(new_mini.transpose()) #with injection info
    return (id_projection_info,id_projection_injection_info)

#dd,cc=id_projection_structure('Primary visual area',None,'Primary visual area, layer 6b',263780729)
#print(dd)
#print(cc.to_string())

def structure_name_info (structure_name) :
    '''
    Returns name info (id, descendants, acronym) for the given structure

    Parameters
    ----------
    structure_name : str
        name of a specific structure we want name info on

    Returns
    -------
    structure_info : Dataframe
        table of all the structures with their id, acronym, descendants ids and acronyms

    structure : Dataframe
        table for the given name structure with its id, acronym and descendants ids and acronyms

    '''

    all_experiments = mcc.get_experiments(dataframe=True)  # a dataframe
    row = all_experiments.index
    structure_tree = mcc.get_structure_tree()
    ia_map = structure_tree.get_id_acronym_map()
    liste_structures = []
    liste_structure_ids = []
    liste_structure_acronym = []
    liste_descendants_ids = []
    liste_descendants_acro = []
    for i in row:
        if all_experiments.loc[i]["structure_name"] not in liste_structures:
            liste_structures.append(all_experiments.loc[i]["structure_name"])
            liste_structure_ids.append(all_experiments.loc[i]['structure_id'])
            liste_structure_acronym.append(all_experiments.loc[i]["structure_abbrev"])
            desc = structure_tree.descendant_ids([ia_map[all_experiments.loc[i]["structure_abbrev"]]])[0]
            liste_descendants_ids.append(desc)
            desc_acro = []
            for j in desc:
                desc_acronyms = structure_tree.get_structures_by_id([j])[0]
                desc_acro.append(desc_acronyms['acronym'])
            liste_descendants_acro.append(desc_acro)
    structure_info = pd.DataFrame(index=liste_structures,
                                  columns=["id", "structure acronym", "descendant_acronyms", "descendant_ids"])
    l = np.arange(0, len(liste_structure_ids))
    structure_info["id"] = liste_structure_ids
    structure_info["structure acronym"] = liste_structure_acronym
    structure_info["descendant_ids"] = liste_descendants_ids
    structure_info["descendant_acronyms"] = liste_descendants_acro
    structure=structure_info.loc[structure_name]
    return (structure_info,structure)

m,n=structure_name_info("Primary visual area") #see m for all the structures (218)
g=n.loc['descendant_ids'][1:]

def plot_matrix (name_injection_structure, targeted_structure_ids, projection_parameter, experiment_id_list, cre,hemisphere) :
    '''
    Returns a plotted matrix of the level of projection volume or density in a structures according to the experiments/injection structure

    Parameters
    ----------
    name_injection_structure : str

    targeted_structure_ids : List of int
        List of the structure ids for which we want to see the level of projection density or volume they have

    projection_parameter : str
        Can be projection_density, projection_volume, normalized_projection_volume or volume (for now)

    experiment_id_list : List or None
        if List, it's a list of the experiment ids (corresponding to the injection structure)

    cre : str, True, False, None
        name of the cre line

    hemisphere : List
        List which can contain 1 (left), 2 (right), and/or 3 (both hemispheres)

    Returns
    -------
    Plot

    '''

    structure_tree = mcc.get_structure_tree()
    if experiment_id_list == None :
        structure_injection = structure_tree.get_structures_by_name([name_injection_structure])[0]
        structure_injection_experiments = mcc.get_experiments(cre=cre,injection_structure_ids=[structure_injection['id']])
        experiment_ids = [e['id'] for e in structure_injection_experiments]
    else :
        experiment_ids = experiment_id_list
    # matrix
    pm = mcc.get_projection_matrix(experiment_ids=experiment_ids,
                                   projection_structure_ids=targeted_structure_ids,
                                   hemisphere_ids=hemisphere,  # right hemisphere, ipsilateral
                                   parameter=projection_parameter)
    # plot :
    row_labels = pm['rows']  # these are just experiment ids
    column_labels = [c['label'] for c in pm['columns']]
    matrix = pm['matrix']
    title='injection in the '+str(name_injection_structure)
    fig, ax = plt.subplots(figsize=(15, 15))
    heatmap = ax.pcolor(matrix, cmap=plt.cm.afmhot)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(matrix.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(matrix.shape[0]) + 0.5, minor=False)

    ax.set_xlim([0, matrix.shape[1]])
    ax.set_ylim([0, matrix.shape[0]])

    # want a more natural, table-like display
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    ax.set_xticklabels(column_labels, minor=False)
    ax.set_yticklabels(row_labels, minor=False)
    plt.title(title)
    plt.show()

#liste=[263780729]
#plot_matrix('Primary visual area',g,'normalized_projection_volume',None,False,[1,2])



#500836840 ; 307297141 ; 503069254 ; 263780729






