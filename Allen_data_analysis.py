import allensdk
import json
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.cell_types_api import CellTypesApi
import webbrowser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ast import literal_eval
import time
from plotnine import ggplot,geom_line,aes,geom_abline,geom_point

from scipy.stats import linregress
import random
import plotly.express as plotly

from allensdk.core.swc import Marker
import matplotlib.patches as mpatches
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor



ctc = CellTypesCache(manifest_file='cell_types/manifest.json')

def get_cells(file_name,species,reconstruction=False, morphology=False): 
    '''
    Download the data you want for a given species

    Parameters
    ----------
    file_name : str
        The name you want to give to you file

    species : str
        Name of the species you want to retreive data : Mouse, Human or All if you want both species in you data

    Returns
    -------
    None (the file is directly saved)

    '''

    ctc = CellTypesCache(manifest_file='cell_types/manifest.json') #create the manifest file
    if species=="Mouse":
        my_species=[CellTypesApi.MOUSE]
        ctc.get_cells(file_name=file_name,require_reconstruction=reconstruction,require_morphology=morphology,species=my_species)
    elif species=="Human":
        my_species=[CellTypesApi.HUMAN]
        ctc.get_cells(file_name=file_name,require_reconstruction=reconstruction,require_morphology=morphology,species=my_species)
    elif species=="All":
        ctc.get_cells(file_name=file_name,require_reconstruction=reconstruction,require_morphology=morphology)
    
    return("Data saved in file: "+str(file_name))
#data=ctc.get_all_features(dataframe=True, require_reconstruction=True)

def dict_specimen(name_file): 
    '''
    Generate lists of dictionary for all specimen per species

    Parameters
    ----------
    name_file : str
        File name entered in get_cells function as string.

    Returns
    -------
    d_mouse : list
        List of Dictionnary for all mouse specimen .
    d_human : TYPE
        List of dictionnary for all human specimen.

    '''
    
    f = open(name_file)
    data = json.load(f)
    d_mouse=[]
    d_human=[]
    for i in data :
        if i["donor__species"] == "Mus musculus":
            d_mouse.append(i)
        else :
            d_human.append(i)
    return (d_mouse,d_human)

def take_id(name_file) : 
    '''
    Generates lists of all specimen id per species

    Parameters
    ----------
    name_file : str
        file name entered in get_cells function as string.

    Returns
    -------
    id_mouse : list
        list containing all the mouse specimen id.
    id_human : list
        list containing all the human specimen id.

    '''
   
    f=open(name_file)
    data=json.load(f)
    id_mouse=[]
    id_human=[]
    for i in data :
        if i["donor__species"]=="Mus musculus":
            id_mouse.append(i["specimen__id"])
        else :
            id_human.append(i["specimen__id"])
    return (id_mouse,id_human)

def cell_info (cell_id, dict_species) :
    '''
    Returns info caracterizing the cell

    Parameters
    ----------
    cell_id : str
        one str value of the id_list from take_id function

    dict_species : Dictionary
        dictionary from take_id function.

    Returns
    -------
    dict_cell_info : Dictionary
        contains as keys the parameters of the cell and as values the values corresponding to those parameters

    '''

    for i in dict_species :
        if i["specimen__id"]==cell_id :
            dict_cell_info=i
            return dict_cell_info


def ephys_web_page(cell_id) :
    '''
    The electrophysiological web page of the cell id in the Allen Brain Atlas is opened

    Parameters
    ----------
    cell_id : str
        one str value of the id_list from take_id function

    '''

    link="http://celltypes.brain-map.org/experiment/electrophysiology/"+str(cell_id)
    webbrowser.open(link)
    
def morpho_web_page(cell_id) :
    '''
    The morphological web page of the cell id in the Allen Brain Atlas is opened

    Parameters
    ----------
    cell_id : str
        one str value of the id_list from take_id function

    '''

    link="http://celltypes.brain-map.org/experiment/morphology/"+str(cell_id)
    webbrowser.open(link)


def structure_dict (dict_species) :
    '''
    Returns the cell number within each layer of all the structures of the species brain

    Parameters
    ----------
    dict_species : Dictionary
        dictionary from take_id function.

    Returns
    -------
    dict_structure : Dictionary
        Dictionary containing the cell number of the layers of the structures and substructures of a given species (keys are the structures/substructures/layers and the values are the cell number)

    '''
    
    i=0
    dict_structure=dict()
    while i<len(dict_species) :
        full_name=dict_species[i]['structure__name']
        j=0
        for e in full_name :
            if e=="," :
                j=j+1
        if j==1 :
            name1,layer1=full_name.split(",")
            name=name1.replace('"','')
            layer=layer1.replace('"','')
            if name not in dict_structure :
                dict_structure[name] = dict()
                dict_structure[name][layer] = 1
            else :
                if layer not in dict_structure[name] :
                    dict_structure[name][layer]=1
                else :
                    dict_structure[name][layer]=dict_structure[name][layer]+1
        elif j==2 :
            
            name1,subname,layer1=full_name.split(",")
            name = name1.replace('"', '')
            layer = layer1.replace('"', '')
            if name not in dict_structure:
                dict_structure[name] = dict()
                dict_structure[name][subname] =dict()
                dict_structure[name][subname][layer]=1
            else :
                if subname not in dict_structure[name] :
                    dict_structure[name][subname] = dict()
                    dict_structure[name][subname][layer] = 1
                else :
                    if layer not in dict_structure[name][subname] :
                        dict_structure[name][subname][layer] = 1
                    else:
                        dict_structure[name][subname][layer]=dict_structure[name][subname][layer]+1
        elif j==0 :
            if full_name not in dict_structure :
                dict_structure[full_name]=1
            else :
                dict_structure[full_name]=dict_structure[full_name]+1
        i=i+1
    return(dict_structure)

def transgenic_line_name (cell_id,dict_species) :
    '''
    Returns the name of the transgenic line

    Parameters
    ----------
    cell_id : str
        the id of the cell.

    dict_species : Dictionary
        dictionary from take_id function.

    Returns
    -------
    line_name : str
        Name of the transgenic line

    '''
    for i in dict_species :
        if i["specimen__id"]==cell_id :
            line_name=i["line_name"]
            return line_name


def structures_transgenic_lines (dict_species) :
    '''
    Returns the cell ids of each transgenic line being in each structure/substructure/layer

    Parameters
    ----------
    dict_species : Dictionary
        dictionary from take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)

    Returns
    -------
    dic_type : Dictionary
        Dictionary containing for each transgenic line being expressed in each structure/substructure/layer a list of the cell ids expressing those transgenic lines

    '''

    dic_type=dict()
    for i in dict_species :
        j=0
        cre=i["line_name"]
        str_name=i["structure__name"]
        for e in str_name :
            if e==',' :
                j=j+1
        if j==1 :
            name,layer=str_name.split(",")
            name=name.replace('"','')
            layer=layer.replace('"','')
            if name not in dic_type :
                dic_type[name]=dict()
                dic_type[name][layer]=dict()
                dic_type[name][layer][cre]= []
                dic_type[name][layer][cre].append(i["specimen__id"])
            else :
                if layer not in dic_type[name] :
                    dic_type[name][layer]=dict()
                    dic_type[name][layer][cre] =  []
                    dic_type[name][layer][cre].append(i["specimen__id"])
                else :
                    if cre not in dic_type[name][layer] :
                        dic_type[name][layer][cre] = []
                        dic_type[name][layer][cre].append(i["specimen__id"])
                    else :
                        dic_type[name][layer][cre].append(i["specimen__id"])
        else :
            if j==2 :
                name,subname,layer=str_name.split(",")
                name = name.replace('"', '')
                layer = layer.replace('"', '')
                if name not in dic_type :
                    dic_type[name]=dict()
                    dic_type[name][subname]=dict()
                    dic_type[name][subname][layer]=dict()
                    dic_type[name][subname][layer][cre]=[]
                    dic_type[name][subname][layer][cre].append(i["specimen__id"])
                else :
                    if subname not in dic_type[name] :
                        dic_type[name][subname] = dict()
                        dic_type[name][subname][layer] = dict()
                        dic_type[name][subname][layer][cre] = []
                        dic_type[name][subname][layer][cre].append(i["specimen__id"])
                    else :
                        if layer not in dic_type[name][subname] :
                            dic_type[name][subname][layer] = dict()
                            dic_type[name][subname][layer][cre] = []
                            dic_type[name][subname][layer][cre].append(i["specimen__id"])
                        else :
                            if cre not in dic_type[name][subname][layer] :
                                dic_type[name][subname][layer][cre] = []
                                dic_type[name][subname][layer][cre].append(i["specimen__id"])
                            else :
                                dic_type[name][subname][layer][cre].append(i["specimen__id"])
            else :
                if j==0 :
                    if str_name not in dic_type :
                        dic_type[str_name]=dict()
                        dic_type[str_name][cre]=[]
                        dic_type[str_name][cre].append(i["specimen__id"])
                    else :
                        if cre not in dic_type[str_name] :
                            dic_type[str_name][cre] = []
                            dic_type[str_name][cre].append(i["specimen__id"])
                        else :
                            dic_type[str_name][cre].append(i["specimen__id"])

    return (dic_type)

def layers_transgenic_lines (dict_species) : #return in a list the structures in str followed by their dictionary containing for each marker the layer and its id
    '''
    Within each structure it returns for each transgenic line all the its ids in each layer

    Parameters
    ----------
    dict_species : Dictionary
        dictionary from the take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)

    Returns
    -------
    lines_of_each_layer : list
        A list of the structures names in str followed with their dictionary containing for each transgenic line the layers and its id

    '''

    dic=structures_transgenic_lines(dict_species)
    lines_of_each_layer=[]
    for i in dic :
        lines_of_each_layer.append(i)
        ss_str=list(dic[i].keys())
        m = dict()
        for j in ss_str :
            m_or_l=list(dic[i][j].keys())
            t=type(dic[i][j][m_or_l[0]]) is dict
            if t==False :
                for k in m_or_l :
                    if k not in m :
                        m[k]=dict()
                        m[k][j]=dic[i][j][k]
                    else :
                        if j not in m[k] :
                            m[k][j] = dic[i][j][k]
            else :
                if t==True:
                    for d in m_or_l :
                        mm=list(dic[i][j][d].keys())
                        for l in mm :
                            if l not in m :
                                m[l]=dict()
                                m[l][j]=dict()
                                m[l][j][d]=dic[i][j][d][l]
        lines_of_each_layer.append(m)
    return(liste)

def spiny_reconstruction_transgenic_line_state_plot (dict_species) : #within each structure of a given species, for each marker it plots an histogram of the proportion of spiny/aspiny/sparsely spiny with or without reconstruction in each layer
    '''
    Within each structure of a given species, for each transgenic line it plots an histogram of its proportion of spiny/aspiny/sparsely spiny with or without reconstruction in each layer

    Parameters
    ----------
    dict_species : Dictionary
        dictionary from the take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)

    '''

    l = markers(dict_specimen)
    for i in l:

        t = type(i) is dict
        if t == False:
            titre = i
        if t == True:
            for j in i:
                aa = 0
                t = type(j) is dict
                if t == False:
                    ss_titre = j
                    layers = list(i[j].keys())
                    len_spi_no = []
                    len_spi_yes = []
                    len_asp_no = []
                    len_asp_yes = []
                    len_spa_no = []
                    len_spa_yes = []
                    spi_no = []
                    spi_yes = []
                    asp_no = []
                    asp_yes = []
                    spa_no = []
                    spa_yes = []
                    a, b, c, d, e, f = spiny_reconstruction(dict_specimen)
                    for k in layers:
                        f = type(i[j][k]) is dict
                        if f == False:
                            ss_title = " "
                            lab = layers
                            ids = i[j][k]
                            for z in ids:
                                if z in a:
                                    spi_no.append(z)
                                else:
                                    if z in b:
                                        spi_yes.append(z)
                                    else:
                                        if z in c:
                                            asp_no.append(z)
                                        else:
                                            if z in d:
                                                asp_yes.append(z)
                                            else:
                                                if z in e:
                                                    spa_no.append(z)
                                                else:
                                                    spa_yes.append(z)
                            len_spi_no.append(len(spi_no))
                            len_spi_yes.append(len(spi_yes))
                            len_asp_no.append(len(asp_no))
                            len_asp_yes.append(len(asp_yes))
                            len_spa_no.append(len(spa_no))
                            len_spa_yes.append(len(spa_yes))
                        else:
                            if f == True:
                                ss_title = ", " + str(k) + ", "
                                for layer in i[j][k]:
                                    lab = list(i[j][k].keys())
                                    ids = i[j][k][layer]
                                    for z in ids:
                                        if z in a:
                                            spi_no.append(z)
                                        else:
                                            if z in b:
                                                spi_yes.append(z)
                                            else:
                                                if z in c:
                                                    asp_no.append(z)
                                                else:
                                                    if z in d:
                                                        asp_yes.append(z)
                                                    else:
                                                        if z in e:
                                                            spa_no.append(z)
                                                        else:
                                                            spa_yes.append(z)
                                    len_spi_no.append(len(spi_no))
                                    len_spi_yes.append(len(spi_yes))
                                    len_asp_no.append(len(asp_no))
                                    len_asp_yes.append(len(asp_yes))
                                    len_spa_no.append(len(spa_no))
                                    len_spa_yes.append(len(spa_yes))
                    grd_titre = titre + ss_title + str(ss_titre)
                    plt.bar(lab, len_spi_no, label='spiny without reconstruction')
                    plt.bar(lab, len_spi_yes, label='spiny with reconstruction')
                    plt.bar(lab, len_asp_no, label='aspiny without reconstruction')
                    plt.bar(lab, len_asp_yes, label='aspiny with reconstruction')
                    plt.bar(lab, len_spa_no, label='sparsely spiny without reconstruction')
                    plt.bar(lab, len_spa_yes, label='sparsely spiny with reconstruction')
                    plt.title(grd_titre)
                    plt.legend()
                    plt.show()

def transgenic_line_number (dict_species) : #return a dictionary with the marker names and how many they are ; the "total" is only a verification that the dic_type as indeed all the id of the dict_specimen
    '''
    Returns for each transgenic line its cell number within the species brain (can only be mouse)

    Parameters
    ----------
    dict_species : Dictionary
        dictionary from the take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)

    Returns
    -------
    dic_transgenic_line : Dictionary
        dictionary in which each key is a transgenic line and the values the number of cells expressing this trangenic line

    '''

    dico=structures_transgenic_lines(dict_species)
    total = 0
    dic_transgenic_line = dict()
    for i in dico:
        t=type(dico[i]) is dict
        if t==False :
            if i not in dic_transgenic_line :
                dic_transgenic_line[i]=len(dico[i])
            else :
                dic_transgenic_line[i]=dic_transgenic_line[i]+len(dico[i])
        for j in dico[i]:
            t = type(dico[i][j]) is dict
            if t == False:
                total = total + len(dico[i][j])
                if j not in dic_transgenic_line:
                    dic_transgenic_line[j] = len(dico[i][j])
                else:
                    dic_transgenic_line[j] = dic_transgenic_line[j] + len(dico[i][j])
            else:
                for k in dico[i][j]:
                    t = type(dico[i][j][k]) is dict
                    if t == False:
                        total = total + len(dico[i][j][k])
                        if k not in dico_marker:
                            dic_transgenic_line[k] = len(dico[i][j][k])
                        else:
                            dic_transgenic_line[k] = dic_transgenic_line[k] + len(dico[i][j][k])

                    else:
                        for z in dico[i][j][k]:
                            t = type(dico[i][j][k][z]) is dict
                            if t == False:
                                total = total + len(dico[i][j][k][z])
                                if z not in dic_transgenic_line :
                                    dic_transgenic_line[z] = len(dico[i][j][k][z])
                                else:
                                    dic_transgenic_line[z] = dic_transgenic_line[z] + len(dico[i][j][k][z])
    return (dic_transgenic_line)

def plot_transgenic_line_number(dict_species) :
    '''
    Plot an histogram of all the transgenic lines in the given species

    Parameters
    ----------
    dict_species : Dictionary
        dictionary from the take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)

    '''

    d=transgenic_line_number(dict_species)
    height=d.values()
    bars=d.keys()
    x_pos=np.arange(1,370,10)
    plt.bar(x_pos, height,width=5)
    plt.xticks(x_pos,bars,rotation=90,fontsize=10)
    specie=dict_species[0]["donor__species"]
    title="markers of the "+specie
    plt.tight_layout()
    plt.title(title)
    plt.show()


def plot_number_structures (dict_structure) : 
    '''
    Plots of the number of cells within the structure/substructure/layer for a given species

    Parameters
    ----------
    dict_structure : Dictionary
        Dictionary from structure function

    Returns
    -------
    None.

    '''
    
    structures=list(dict_structure.keys())
    test=type(dict_structure[structures[0]]) is dict
    if test==False :
        li=list(dict_structure.keys())
        vall=list(dict_structure.values())
        plt.bar(li,vall,0.35)
        plt.xticks(rotation=15)
        plt.title("human")
        plt.show()
    else :
        for i in structures :
            ss_structure=list(dict_structure[i].keys())
            t=type(dict_structure[i][ss_structure[0]]) is dict
            if t==True :
                j=0
                r=0
                while j<len(ss_structure) :
                    if type(dict_structure[i][ss_structure[j]])is dict:
                        labels=list(dict_structure[i][ss_structure[j]].keys())
                        val=list(dict_structure[i][ss_structure[j]].values())
                        width=0.35+r
                        plt.bar(labels, val, width, label=str(ss_structure[j]))
                        plt.title(str(i))
                        plt.legend()
                    else :
                        label = ss_structure[j]
                        val=dict_structure[i][ss_structure[j]]
                        width = 0.35
                        plt.bar(label, val, width, label=str(ss_structure[j]))
                        plt.title(str(i))
                        plt.legend()
                    j=j+1
                    r=r-0.1
                plt.show()

            else :
                labels = ss_structure
                val=list(dict_structure[i].values())
                width = 0.35
                fig, ax = plt.subplots()
                ax.bar(labels, val, width, label='layers')
                ax.set_title(str(i))
            plt.show()

def number_of_sweeps (cell_id) :
    '''
    Returns the total number of sweeps for a given cell id

    Parameters
    ----------
    cell_id : str
        the id of the cell.

    Returns
    -------
    number : integer
        the number of sweeps

    '''

    all_sweeps = ctc.get_ephys_sweeps(cell_id)
    number=len(all_sweeps)
    return number

    
def cell_ephys_info (cell_id) : #return the cell electrophysiological info for a given id
    '''
    Returns the general electrophysiological info of a given cell

    Parameters
    ----------
    cell_id : str
        the id of the cell.

    Returns
    -------
    i : Dictionary
        dictionary in which the keys are the electrophysiological parameters and the values are the values

    '''

    all_cell_features = ctc.get_ephys_features()
    for i in all_cell_features :
        if i["specimen_id"]==cell_id :
            return i
        
        
def sweep_info (cell_id,sweep_nber): #return data info of a given sweep of a given id specimen
    '''
    For a given cell and its given sweep number it returns all the info of the sweep


    Parameters
    ----------
    cell_id : str
        the id of the cell.
    sweep_nbre : integer

    Returns
    -------
    i : dictionary
        Dictionary in which the keys are the parameters of the sweep

    '''

    data_set = ctc.get_ephys_sweeps(cell_id)
    for i in data_set :
        if i['sweep_number']==sweep_nber:
            return i
        
        
def plot_stimulus (cell_id,sweep) :
    '''
    Plot the stimulus (in pA) and the cell response (in mV)

    Parameters
    ----------
    cell_id : str
        cell_id found in the id_list from take_id function..
    sweep : str
        sweep number found in stimulus_type function 

    Returns
    -------
    None.

    '''
    
    data_set = ctc.get_ephys_data(cell_id)
    all_sweeps = ctc.get_ephys_sweeps(cell_id)
    for j in all_sweeps :
        if j["sweep_number"]==sweep :
            if j["stimulus_name"]=="Test":
                print ("no plot")
            else:
                sweep_data = data_set.get_sweep(sweep)
                index_range = sweep_data["index_range"]
                i = sweep_data["stimulus"][0:index_range[1] + 1]  # in A
                v = sweep_data["response"][0:index_range[1] + 1]  # in V
                i *= 1e12  # to pA
                v *= 1e3  # to mV
                sampling_rate = sweep_data["sampling_rate"]  # in Hz
                t = np.arange(0, len(v)) * (1.0 / sampling_rate)
                plt.style.use('ggplot')
                fig, axes = plt.subplots(2, 1, sharex=True)
                axes[0].plot(t, v, color='black')
                axes[1].plot(t, i, color='gray')
                axes[0].set_ylabel("mV")
                axes[1].set_ylabel("pA")
                axes[1].set_xlabel("seconds")
                plt.title(j['stimulus_name'])
                plt.show()

def stimulus_type (cell_id):
    '''
    Returns for the different stimuli the sweeps using those stimuli

    Parameters
    ----------
    cell_id : str
        cell_id found in the id_list from take_id function.

    Returns
    -------
    dict_stimuli : dictionary
        dictionary with the keys being the different stimulus types and the values are the lists of sweep numbers.

    '''
    all_sweeps = ctc.get_ephys_sweeps(cell_id)
    dict_stimuli=dict()
    for i in all_sweeps :
        number = i["sweep_number"]
        if i["stimulus_name"] not in dict_stimuli :
            dict_stimuli[i["stimulus_name"]]=[number]
        else :
            dict_stimuli[i["stimulus_name"]].append(number)
    return (dict_stimuli)

def get_structure_name(species_dict,speciment_index): 
    '''
    Return structure substructure and layer for create_species_table function

    Parameters
    ----------
    species_dict : Dictionnary
        From dict_specimen function.
    speciment_index : int
        Index of the specimen of interest.

    Returns
    -------
    name : str

    substructure : str

    layer : str


    '''
    
    full_name=species_dict[speciment_index]['structure__name'] #full_name=dict_mouse[i]["structure__name"]
    j=0
    for e in full_name :
        if e=="," :
            j=j+1
    if j==0:
        name=full_name.replace('"','')
        substructure="No substructure"
        layer="No layer"
    elif j==1:
        name1,layer1=full_name.split(",")
        name=name1.replace('"','')
        layer=layer1.replace('"','')
        substructure="No substructure"
        
    elif j==2:
        name1,substructure,layer1=full_name.split(",")
        name = name1.replace('"', '')
        layer = layer1.replace('"', '')
    return name,substructure,layer



def create_species_sweeps_stim_table(species_dict,from_csv=False):
    '''

    Either from species dict (default) or from already saved table(from_csv=True)
    Create a dataframe containing for a given species in each row
    the specimen_id, 
    the structure acronyme, 
    the structure name, 
    the substructure name, 
    the layer, 
    the dendrite_type,
    the total number of sweeps and 
    the sweeps_id per stimulus kind 

    Parameters
    ----------
    species_dict : Dictionnary or str
        From dict_specimen function
        Or path to csv file
    from_csv : Bool
        Default False: Create table from a csv file

    Returns
    -------
    full_dataframe : DataFrame
        

    '''
    if from_csv==True:
        full_dataframe=pd.read_csv(str(species_dict))
        full_dataframe=full_dataframe.iloc[:,1:]
        
        for row in range(full_dataframe.shape[0]):
            for col in range(7,20):
                if type(full_dataframe.iloc[row,col])== str:
                    full_dataframe.iat[row,col]=literal_eval(full_dataframe.iat[row,col])
                 
                else:
                    full_dataframe.iloc[row,col]=-1
        return full_dataframe
    
    start=time.time()
    full_dataframe=pd.DataFrame()

    for current_specimen in range(len(species_dict)):
        name,substructure,layer=get_structure_name(species_dict, current_specimen)
        number_of_sweeps=0
        
        specimen_stimulus_dict=stimulus_type(species_dict[current_specimen]["specimen__id"])
        

        for key in stimulus_type(species_dict[current_specimen]["specimen__id"]).keys():
            number_of_sweeps+=len(stimulus_type(species_dict[current_specimen]["specimen__id"])[key])
            
        current_row=pd.Series([species_dict[current_specimen]["specimen__id"],
                               species_dict[current_specimen]["structure__acronym"],
                               name,
                               substructure,
                               layer,
                               species_dict[current_specimen]["tag__dendrite_type"],
                               number_of_sweeps],index=['specimen_id','structure_acronyme','structure_name','structure_substructure','layer','dendrite_type','number_of_sweeps'])
        current_row=pd.DataFrame(current_row).T
       

        stim_dict=pd.Series(specimen_stimulus_dict)
        stim_dict=pd.DataFrame(stim_dict).T

        full_row=pd.concat([current_row,stim_dict],axis=1)

        
        full_dataframe=pd.concat([full_dataframe,full_row],ignore_index=True,axis=0,join='outer')
    end=time.time()
    print('Running time= ',end-start)
    return full_dataframe



def pie_plotly(full_species_dict,count_what,based_on):
    #fig=plotly.sunburst(mouse_full_table,path=['structure_name','layer','number of sweeps'])
    my_table=pd.DataFrame(columns=[str(based_on),str(count_what)])
    
    for current_elt in list(set(full_species_dict[based_on])):  
        number=(full_species_dict.loc[full_species_dict[based_on]==current_elt,[count_what]].sum(axis=0))[0]
        


        new_line=pd.Series([current_elt,number],index=[str(based_on),str(count_what)])

        my_table=my_table.append(new_line,ignore_index=True)
       
    figu=plotly.pie(my_table,values=len(count_what),names=based_on)
    figu.update_traces(textposition='inside',textinfo='percent+label')
    figu.show()




            
def create_subtable(full_species_table,which_factor,which_values):
    '''
    Create a subtable from full_species_table by selecting the rows by values

    Parameters
    ----------
    full_species_table : DataFrame
        Full table coming from create_species_sweeps_stim_table function. Each sublist in which_factor must have a corresponding sublist in which_values
    which_factor : list of list
        A list containing lists of str of factor to sort the table by: [["structure name"],["layer"],["Ramp"]].
    which_values : List of list
        A list of list of str of value to keep in each factor. Must be -1 if corresponding to a stimulus [["Primary visual area","Anterolateral visual area"],["layer 5","layer 4"],[-1]]

    Returns
    -------
    new_table : DataFrame
        DESCRIPTION.

    '''
    stim_list=full_species_table.columns[7:]
    if len(which_factor)==1:
        if which_factor[0][0] in stim_list:
            new_table=full_species_table[~full_species_table[which_factor[0][0]].isin (which_values[0])]
            return(new_table)
        else:
            new_table=full_species_table[full_species_table[which_factor[0][0]].isin (which_values[0])]
            return(new_table)
    else:
        
        if which_factor[0][0] in stim_list:
            new_table=full_species_table[~full_species_table[which_factor[0][0]].isin (which_values[0])]
            which_factor=which_factor[1:]
            which_values=which_values[1:]
            new_table=create_subtable(new_table, which_factor,which_values)
            return(new_table)
            
        else :
            new_table=full_species_table[full_species_table[which_factor[0][0]].isin (which_values[0])]
            which_factor=which_factor[1:]
            which_values=which_values[1:]
            new_table=create_subtable(new_table, which_factor,which_values)
            return (new_table)
    
        
def select_specimen_sweep_stimulus(species_table,stimulus):
    '''
    

    Parameters
    ----------
    species_table : DataFrame
        DESCRIPTION.
    stimulus : str
        Select the stimulus you want to keep .

    Returns
    -------
    DataFrame
        Sub DataFrame from the original one containing only the columns of the specimen_id and the column of the selected stimulus.

    '''
    return species_table[["specimen_id",stimulus]]
    
def extract_stim_freq(specimen_id,species_table,stimulus,mydef=True):
    '''
    Function to extract for each specified specimen_id and the corresponding stimulus the frequency of the response
    

    Parameters
    ----------
    specimen_id : List of specimen id
    i.e.[623960880,623960824,...]
    species_table : DataFrame
        DataFrame containing at least a column with the specimen_id and the stimulus of interest.
    stimulus : str
        the stimulus from which we want to extract the frequency.

    Returns
    -------
    f_I_table : DataFrame
        DataFrame with a column "specimen_id"(factor),the sweep number (int),the stimulus amplitude in pA(float),and the computed frequency of the response (float).

    '''
    f_I_table=pd.DataFrame(columns=['specimen','sweep','stim_amplitude_pA','frequence_Hz'])
    index_stim=species_table.columns.get_loc(stimulus)
    for current_specimen in specimen_id:
        index_specimen=species_table.index[species_table["specimen_id"]==current_specimen][0]
        
        my_specimen_data=ctc.get_ephys_data(current_specimen)
        sweep_numbers=species_table.iloc[index_specimen,index_stim]
        
        for current_sweep in sweep_numbers:
            
            nb_spike=len(my_specimen_data.get_spike_times(current_sweep))
            if nb_spike <2:
                freq=nb_spike
            elif mydef==True:
                t_first_spike=my_specimen_data.get_spike_times(current_sweep)[0]
                t_last_spike=my_specimen_data.get_spike_times(current_sweep)[-1]
                freq=nb_spike/(t_last_spike-t_first_spike)
            elif mydef==False:
                freq=nb_spike
            new_line=pd.Series([int(current_specimen),current_sweep,my_specimen_data.get_sweep_metadata(current_sweep)['aibs_stimulus_amplitude_pa'],freq],
                                index=['specimen','sweep','stim_amplitude_pA','frequence_Hz'])
            f_I_table=f_I_table.append(new_line,ignore_index=True)
            
    f_I_table=f_I_table.sort_values(by=["specimen",'stim_amplitude_pA'])
    f_I_table['specimen']=pd.Categorical(f_I_table['specimen'])
    return f_I_table

def average_rate(t, spikes, start, end):
    """Calculate average firing rate during interval between `start` and `end`.

    Parameters
    ----------
    t : numpy array of times in seconds
    spikes : numpy array of spike indexes
    start : start of time window for spike detection
    end : end of time window for spike detection

    Returns
    -------
    avg_rate : average firing rate in spikes/sec
    """

    if start is None:
        start = t[0]

    if end is None:
        end = t[-1]

    spikes_in_interval = [spk for spk in spikes if t[spk] >= start and t[spk] <= end]
    avg_rate = len(spikes_in_interval) / (end - start)
    return avg_rate


def fit_specimen_fi_slope(stim_amps, avg_rates):
    """
    Fit the rate and stimulus amplitude to a line and return the slope of the fit.
    
    Parameters
    ----------
    stim_amps: array of sweeps amplitude in mA
    avg_rates: array of sweeps avergae firing rate in Hz

    Returns
    -------
    m: f-I curve slope for the specimen
    c:f-I curve intercept for the specimen
        
    """

    

    x = stim_amps
    y = avg_rates

    A = np.vstack([x, np.ones_like(x)]).T
    m, c = np.linalg.lstsq(A, y,rcond=None)[0]

    return m,c

def plot_f_I_curve(f_I_table,single_specimen=None):
    '''
    Render a plot of the f-I curve with the gain estimation. Also render a table with in each line the specimen, the slope and the intercept
    Parameters
    ----------
    f_I_table : DataFrame
        Comming from the extract_stim_freq function.
    single_specimen: int
        If you want to display only one specimen data, specify its id

    Returns
    -------
    Plot of the response with a fitted linear trend .

    '''
    if single_specimen != None:
        f_I_table=f_I_table[f_I_table['specimen']==single_specimen]
    linear_fit_param=pd.DataFrame(columns=['specimen','slope','intercept'])
    myplot=ggplot()
    for specimen in f_I_table['specimen'].cat.categories:
       
        sub_specimen_table=f_I_table.loc[f_I_table["specimen"]==specimen]

        sub_table=sub_specimen_table.loc[sub_specimen_table["frequence_Hz"]>0]

        slope,intercept=fit_specimen_fi_slope(np.array(sub_table["stim_amplitude_pA"]),np.array(sub_table['frequence_Hz']))
        new_line=pd.Series([specimen,slope,intercept],index=['specimen','slope','intercept'])
        
        linear_fit_param=linear_fit_param.append(new_line,ignore_index=True)
        
    linear_fit_param['specimen']=pd.Categorical(linear_fit_param['specimen'])
    myplot=ggplot(f_I_table,aes(x="stim_amplitude_pA",y="frequence_Hz",color="specimen"))+geom_point()+geom_abline(linear_fit_param,aes(intercept="intercept",slope="slope",colour="specimen"))
        
    return(myplot,linear_fit_param)
   # my_sweep=my_specimen_data.get_sweep(sweep_nb)
   
def random_id_list(id_list,n):
    
   max_value= len(id_list)
   new_id_list=[]
   for elt in range(n):
       index=random.randint(1,max_value-1)

       new_id_list.append(id_list[index])
   return new_id_list

def ephys_spike_key_extractor(cell_id,sweep,t_start,t_end): #return the keys/feature names caracterizing the spike ; if don't want to put t_start/end, make them equal to None
    """
    For a given cell and sweep it returns all the parameters (feature names) caracterizing the spikes

    Parameters
    ----------
    cell_id : str
    sweep : sweep number
    t_start : start of time window for spike detection (if no value, put None)
    t_end : end of time window for spike detection (if no value, put None)

    Returns
    -------
    liste : list of the parameters caracterizing the spikes
    """

    dict_stimulus = stimulus_type(cell_id)
    if sweep in dict_stimulus['Test']:
        return ('no data to return')
    else:
        data_set = ctc.get_ephys_data(cell_id)
        sweep_data = data_set.get_sweep(sweep)
        print(sweep_data["index_range"])
        index_range = sweep_data["index_range"]
        i = sweep_data["stimulus"][0:index_range[1] + 1]  # in A = values of the stimulus
        v = sweep_data["response"][0:index_range[1] + 1]  # in V = values of the response
        i *= 1e12  # conversion to pA
        v *= 1e3  # conversion to mV

        sampling_rate = sweep_data["sampling_rate"]  # in Hz
        t = np.arange(0, len(v)) * (1.0 / sampling_rate)  # scale of the time ?

        sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, start=t_start, end=t_end)
        sweep_ext.process_spikes()  # à mettre avant d'utiliser les autres fonctions related to the spikes
        liste = sweep_ext.spike_feature_keys()
        return (liste)





def spike_info (cell_id,sweep_number,key) : #with a given feature name, and for a given sweep of a given id specimen, this function returns the data of the feature
    """
    For a given cell, sweep and parameter of the spikes, it returns the info of this parameter

    Parameters
    ----------
    cell_id : str
    sweep_number : str
    key : str
        the spike parameter

    Returns
    -------
    info : list of the value of the spike parameter
    """

    data_set = ctc.get_ephys_data(cell_id)
    sweep_data = data_set.get_sweep(sweep_number)
    index_range = sweep_data["index_range"]
    i = sweep_data["stimulus"][0:index_range[1] + 1]  # in A = values of the stimulus
    v = sweep_data["response"][0:index_range[1] + 1]  # in V = values of the response
    i *= 1e12  # conversion to pA
    v *= 1e3  # conversion to mV

    sampling_rate = sweep_data["sampling_rate"]  # in Hz
    t = np.arange(0, len(v)) * (1.0 / sampling_rate)  # scale of the time ?

    sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i)  # start & end give the interval in which we extract
    sweep_ext.process_spikes()
    info=sweep_ext.spike_feature(key)
    return info

#Morphology :

def spiny_state (dict_species) :
    """
    Returns the cell ids being spiny, sparsely spiny or aspiny.

    Parameters
    ----------
    dict_species : Dictionary

    Returns
    -------
    spiny : list
        list of cell ids bein spiny

    aspiny : list
        list of cell ids being aspiny

    sparsely_spiny : list
        list of cell ids being sparsely spiny
    """

    spiny=[]
    aspiny=[]
    sparsely_spiny=[]
    for i in dict_species :
        if i['tag__dendrite_type']=='spiny' :
            spiny.append(i['specimen__id'])
        else :
            if i['tag__dendrite_type'] == 'aspiny':
                aspiny.append(i['specimen__id'])
            else :
                sparsely_spiny.append(i['specimen__id'])
    return(spiny,aspiny,sparsely_spiny)

def reconstruction(dict_species) :
    """
    Returns the cell ids having and having not a morphological reconstruction

    Parameters
    ----------
    dict_species : Dictionary

    Returns
    -------
    no : list
        list of cell ids having not a reconstruction

    yes : list
        list of cell ids having a reconstruction
    """

    no=[]
    yes=[]
    for i in dict_species :
        if i["nr__reconstruction_type"]==None :
            no.append(i["specimen__id"])
        else :
            yes.append(i["specimen__id"])
    return (no,yes)

def spiny_reconstruction (dict_species) :
    """
    Returns the cell ids according to if they are spiny/aspiny/sparsely spiny and have or not a reconstruction.

    Parameters
    ----------
    dict_species : Dictionary

    Returns
    -------
    spiny_no_reconstruction : list
        list of cell ids being spiny without a reconstruction

    spiny_reconstruction : list
        list of cell ids being spiny with a reconstruction

    aspiny_no_reconstruction : list
        list of cell ids being aspiny without a reconstruction

    aspiny_reconstruction :list
        list cell ids being aspiny with a reconstruction

    sparsely_no_reconstruction : list
        list of cell ids being sparsely spiny without a reconstruction

    sparsely_reconstruction :list
        list of cell ids being sparsely spiny with a reconstruction

    """

    spiny_no_reconstruction=[]
    spiny_reconstruction=[]
    aspiny_no_reconstruction=[]
    aspiny_reconstruction=[]
    sparsely_no_reconstruction=[]
    sparsely_reconstruction=[]
    no,yes=reconstruction(dict_specimen)
    spiny,aspiny,sparsely=spiny_state(dict_species)
    for i in yes :
        if i in spiny :
            spiny_reconstruction.append(i)
        else :
            if i in aspiny :
                aspiny_reconstruction.append(i)
            else :
                sparsely_reconstruction.append(i)
    for i in no :
        if i in spiny :
            spiny_no_reconstruction.append(i)
        else :
            if i in aspiny :
                aspiny_no_reconstruction.append(i)
            else :
                sparsely_no_reconstruction.append(i)
    return(spiny_no_reconstruction,spiny_reconstruction,aspiny_no_reconstruction,aspiny_reconstruction,sparsely_no_reconstruction,sparsely_reconstruction)

def morpho_info(cell_id): #for a given id, returns the morpho info of the cell
    """
    Returns the morphological info of this cell

    Parameters
    ----------
    cell_id : str

    Returns
    -------
    i : Dictionary

    """

    all_morpho = ctc.get_morphology_features()
    for i in all_morpho:
        if i["specimen_id"] == cell_id:
            return i

def swc_morpho(cell_id,name) : #return a link to download the morphology swc file of a given cell
    """
    Returns the swc file of the cell

    Parameters
    ----------
    cell_id : str
    name : str
        name we want to give the file ?

    Returns
    -------
    swc_file : a file + a link

    """
    swc_file=ctc.get_reconstruction(cell_id,name)
    return(swc_file)

def compartment_nber (cell_id) :
    """
    For a given cell it returns the number of all its compartments

    Parameters
    ----------
    cell_id : str

    Returns
    -------
    compartment_number : int

    """

    all_morpho=ctc.get_reconstruction(cell_id)
    compartment_number=len(all_morpho.compartment_list)
    return(compartment_number)


def twoD_morpho(cell_id):
    """
    For a given cell it plots it morphology in different spaces (but always in 2D)

    Parameters
    ----------
    cell_id : str

    Returns
    -------
    None

    """

    morphology = ctc.get_reconstruction(cell_id)
    markers = ctc.get_reconstruction_markers(cell_id)
    fig, axes = plt.subplots(1, 2, sharey=True, sharex=True)
    axes[0].set_aspect('equal', 'box')
    axes[1].set_aspect('equal', 'box')

    # Make a line drawing of x-y and y-z views
    for n in morphology.compartment_list:
        for c in morphology.children_of(n):#print the different children compartements
            if n['type']==2:
                axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='red')
                axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='red')
            else :
                if n['type']==1 : #soma
                    axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='green')
                    axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='green')
                else :
                    axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='black')
                    axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='black')

    # cut dendrite markers
    dm = [m for m in markers if m['name'] == Marker.CUT_DENDRITE] #if m['name'] ==10

    axes[0].scatter([m['x'] for m in dm], [m['y'] for m in dm], color='#3333ff')
    axes[1].scatter([m['z'] for m in dm], [m['y'] for m in dm], color='#3333ff')

    # no reconstruction markers
    nm = [m for m in markers if m['name'] == Marker.NO_RECONSTRUCTION] #if m['name'] ==20

    axes[0].scatter([m['x'] for m in nm], [m['y'] for m in nm], color='#333333')
    axes[1].scatter([m['z'] for m in nm], [m['y'] for m in nm], color='#333333')

    axes[0].set_ylabel('y')
    axes[0].set_xlabel('x')
    axes[1].set_xlabel('z')
    red_patch = mpatches.Patch(color='red', label='axons')
    black_patch=mpatches.Patch(color='black', label='dendrites')
    green_patch=mpatches.Patch(color='green', label='soma')
    blue_patch=mpatches.Patch(color='blue', label='location of the dendrite truncations')
    grey_patch=mpatches.Patch(color='grey', label='location of the axon truncations')
    axes[1].legend(handles=[red_patch,black_patch,green_patch,blue_patch,grey_patch],bbox_to_anchor=(0,1.06,1,0.2))
    plt.show()


    
# ctc = CellTypesCache(manifest_file='cell_types/manifest.json') #create the manifest file

# cells = ctc.get_cells(file_name="all_cells",require_reconstruction=True)
# id_mouse,id_human=take_id("all_cells")
#dict_mouse,dict_human=dict_specimen("all_cells")
#a=structure(dict_human)
#print(a)
#plot_number_structures(a)
#a=stimulus_type(464212183)
#print(a)
