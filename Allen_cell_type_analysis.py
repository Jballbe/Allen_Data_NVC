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
from plotnine import ggplot, geom_line, aes, geom_abline, geom_point, geom_text

from scipy.stats import linregress
from scipy import optimize
from scipy.optimize import curve_fit
import random
import plotly.express as plotly

from allensdk.core.swc import Marker
import matplotlib.patches as mpatches
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor
import warnings

#change '.ix' to '.loc' in allensdk\ephys\ephys_extractor.py (to avoid error for the .process_spikes() function

#cell types data analysis

ctc = CellTypesCache(manifest_file='cell_types/manifest.json')


def get_cells(file_name, species, reconstruction=False, morphology=False):
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

    ctc = CellTypesCache(manifest_file='cell_types/manifest.json')  # create the manifest file
    if species == "Mouse":
        my_species = [CellTypesApi.MOUSE]
        ctc.get_cells(file_name=file_name, require_reconstruction=reconstruction, require_morphology=morphology,
                      species=my_species)
    elif species == "Human":
        my_species = [CellTypesApi.HUMAN]
        ctc.get_cells(file_name=file_name, require_reconstruction=reconstruction, require_morphology=morphology,
                      species=my_species)
    elif species == "All":
        ctc.get_cells(file_name=file_name, require_reconstruction=reconstruction, require_morphology=morphology)

    return ("Data saved in file: " + str(file_name))
    #return (ctc)


# data=ctc.get_all_features(dataframe=True, require_reconstruction=True)

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
        List of Dictionnary for all mouse specimen
    d_human : list
        List of dictionnary for all human specimen.
    '''

    f = open(name_file)
    data = json.load(f)
    d_mouse = []
    d_human = []
    for i in data:
        if i["donor__species"] == "Mus musculus":
            d_mouse.append(i)
        else:
            d_human.append(i)
    return (d_mouse, d_human)


def take_id(name_file):
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

    f = open(name_file)
    data = json.load(f)
    id_mouse = []
    id_human = []
    for i in data:
        if i["donor__species"] == "Mus musculus":
            id_mouse.append(i["specimen__id"])
        else:
            id_human.append(i["specimen__id"])
    return (id_mouse, id_human)


def cell_info(cell_id, dict_species):
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

    for i in dict_species:
        if i["specimen__id"] == cell_id:
            dict_cell_info = i
            return dict_cell_info


def ephys_web_page(cell_id):
    '''
    The electrophysiological web page of the cell id in the Allen Brain Atlas is opened
    Parameters
    ----------
    cell_id : str
        one str value of the id_list from take_id function
    '''

    link = "http://celltypes.brain-map.org/experiment/electrophysiology/" + str(cell_id)
    webbrowser.open(link)


def morpho_web_page(cell_id):
    '''
    The morphological web page of the cell id in the Allen Brain Atlas is opened
    Parameters
    ----------
    cell_id : str
        one str value of the id_list from take_id function
    '''

    link = "http://celltypes.brain-map.org/experiment/morphology/" + str(cell_id)
    webbrowser.open(link)


def structure_dict(dict_species):
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

    i = 0
    dict_structure = dict()
    while i < len(dict_species):
        full_name = dict_species[i]['structure__name']
        j = 0
        for e in full_name:
            if e == ",":
                j = j + 1
        if j == 1:
            name1, layer1 = full_name.split(",")
            name = name1.replace('"', '')
            layer = layer1.replace('"', '')
            if name not in dict_structure:
                dict_structure[name] = dict()
                dict_structure[name][layer] = 1
            else:
                if layer not in dict_structure[name]:
                    dict_structure[name][layer] = 1
                else:
                    dict_structure[name][layer] = dict_structure[name][layer] + 1
        elif j == 2:

            name1, subname, layer1 = full_name.split(",")
            name = name1.replace('"', '')
            layer = layer1.replace('"', '')
            if name not in dict_structure:
                dict_structure[name] = dict()
                dict_structure[name][subname] = dict()
                dict_structure[name][subname][layer] = 1
            else:
                if subname not in dict_structure[name]:
                    dict_structure[name][subname] = dict()
                    dict_structure[name][subname][layer] = 1
                else:
                    if layer not in dict_structure[name][subname]:
                        dict_structure[name][subname][layer] = 1
                    else:
                        dict_structure[name][subname][layer] = dict_structure[name][subname][layer] + 1
        elif j == 0:
            if full_name not in dict_structure:
                dict_structure[full_name] = 1
            else:
                dict_structure[full_name] = dict_structure[full_name] + 1
        i = i + 1
    return (dict_structure)


def transgenic_line_name(cell_id, dict_species):
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
    for i in dict_species:
        if i["specimen__id"] == cell_id:
            line_name = i["line_name"]
            return line_name


def structures_transgenic_lines(dict_species):
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

    dic_type = dict()
    for i in dict_species:
        j = 0
        cre = i["line_name"]
        str_name = i["structure__name"]
        for e in str_name:
            if e == ',':
                j = j + 1
        if j == 1:
            name, layer = str_name.split(",")
            name = name.replace('"', '')
            layer = layer.replace('"', '')
            if name not in dic_type:
                dic_type[name] = dict()
                dic_type[name][layer] = dict()
                dic_type[name][layer][cre] = []
                dic_type[name][layer][cre].append(i["specimen__id"])
            else:
                if layer not in dic_type[name]:
                    dic_type[name][layer] = dict()
                    dic_type[name][layer][cre] = []
                    dic_type[name][layer][cre].append(i["specimen__id"])
                else:
                    if cre not in dic_type[name][layer]:
                        dic_type[name][layer][cre] = []
                        dic_type[name][layer][cre].append(i["specimen__id"])
                    else:
                        dic_type[name][layer][cre].append(i["specimen__id"])
        else:
            if j == 2:
                name, subname, layer = str_name.split(",")
                name = name.replace('"', '')
                layer = layer.replace('"', '')
                if name not in dic_type:
                    dic_type[name] = dict()
                    dic_type[name][subname] = dict()
                    dic_type[name][subname][layer] = dict()
                    dic_type[name][subname][layer][cre] = []
                    dic_type[name][subname][layer][cre].append(i["specimen__id"])
                else:
                    if subname not in dic_type[name]:
                        dic_type[name][subname] = dict()
                        dic_type[name][subname][layer] = dict()
                        dic_type[name][subname][layer][cre] = []
                        dic_type[name][subname][layer][cre].append(i["specimen__id"])
                    else:
                        if layer not in dic_type[name][subname]:
                            dic_type[name][subname][layer] = dict()
                            dic_type[name][subname][layer][cre] = []
                            dic_type[name][subname][layer][cre].append(i["specimen__id"])
                        else:
                            if cre not in dic_type[name][subname][layer]:
                                dic_type[name][subname][layer][cre] = []
                                dic_type[name][subname][layer][cre].append(i["specimen__id"])
                            else:
                                dic_type[name][subname][layer][cre].append(i["specimen__id"])
            else:
                if j == 0:
                    if str_name not in dic_type:
                        dic_type[str_name] = dict()
                        dic_type[str_name][cre] = []
                        dic_type[str_name][cre].append(i["specimen__id"])
                    else:
                        if cre not in dic_type[str_name]:
                            dic_type[str_name][cre] = []
                            dic_type[str_name][cre].append(i["specimen__id"])
                        else:
                            dic_type[str_name][cre].append(i["specimen__id"])

    return (dic_type)


def layers_transgenic_lines(
        dict_species):  # return in a list the structures in str followed by their dictionary containing for each marker the layer and its id
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

    dic = structures_transgenic_lines(dict_species)
    lines_of_each_layer = []
    for i in dic:
        lines_of_each_layer.append(i)
        ss_str = list(dic[i].keys())
        m = dict()
        for j in ss_str:
            m_or_l = list(dic[i][j].keys())
            t = type(dic[i][j][m_or_l[0]]) is dict
            if t == False:
                for k in m_or_l:
                    if k not in m:
                        m[k] = dict()
                        m[k][j] = dic[i][j][k]
                    else:
                        if j not in m[k]:
                            m[k][j] = dic[i][j][k]
            else:
                if t == True:
                    for d in m_or_l:
                        mm = list(dic[i][j][d].keys())
                        for l in mm:
                            if l not in m:
                                m[l] = dict()
                                m[l][j] = dict()
                                m[l][j][d] = dic[i][j][d][l]
        lines_of_each_layer.append(m)
    return (lines_of_each_layer)


def spiny_reconstruction_transgenic_line_state_plot(
        dict_species):  # within each structure of a given species, for each marker it plots an histogram of the proportion of spiny/aspiny/sparsely spiny with or without reconstruction in each layer
    '''
    Within each structure of a given species, for each transgenic line it plots an histogram of its proportion of spiny/aspiny/sparsely spiny with or without reconstruction in each layer
    Parameters
    ----------
    dict_species : Dictionary
        dictionary from the take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)
    '''

    l = layers_transgenic_lines(dict_species)
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
                    a, b, c, d, e, f = spiny_reconstruction(dict_species)
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


def transgenic_line_number(dict_species):
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

    dico = structures_transgenic_lines(dict_species)
    total = 0
    dic_transgenic_line = dict()
    for i in dico:
        t = type(dico[i]) is dict
        if t == False:
            if i not in dic_transgenic_line:
                dic_transgenic_line[i] = len(dico[i])
            else:
                dic_transgenic_line[i] = dic_transgenic_line[i] + len(dico[i])
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
                        if k not in dic_transgenic_line:
                            dic_transgenic_line[k] = len(dico[i][j][k])
                        else:
                            dic_transgenic_line[k] = dic_transgenic_line[k] + len(dico[i][j][k])

                    else:
                        for z in dico[i][j][k]:
                            t = type(dico[i][j][k][z]) is dict
                            if t == False:
                                total = total + len(dico[i][j][k][z])
                                if z not in dic_transgenic_line:
                                    dic_transgenic_line[z] = len(dico[i][j][k][z])
                                else:
                                    dic_transgenic_line[z] = dic_transgenic_line[z] + len(dico[i][j][k][z])
    return (dic_transgenic_line)


def plot_transgenic_line_number(dict_species):
    '''
    Plot an histogram of all the transgenic lines in the given species
    Parameters
    ----------
    dict_species : Dictionary
        dictionary from the take_id function, can only be a dictionary containing mouse cells (no human cells since they don't have any line name)
    '''

    d = transgenic_line_number(dict_species)
    height = d.values()
    bars = d.keys()
    x_pos = np.arange(1, 370, 10)
    plt.bar(x_pos, height, width=5)
    plt.xticks(x_pos, bars, rotation=90, fontsize=10)
    specie = dict_species[0]["donor__species"]
    title = "markers of the " + specie
    plt.tight_layout()
    plt.title(title)
    plt.show()


def plot_number_structures(dict_structure):
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

    structures = list(dict_structure.keys())
    test = type(dict_structure[structures[0]]) is dict
    if test == False:
        li = list(dict_structure.keys())
        vall = list(dict_structure.values())
        plt.bar(li, vall, 0.35)
        plt.xticks(rotation=15)
        plt.title("human")
        plt.show()
    else:
        for i in structures:
            ss_structure = list(dict_structure[i].keys())
            t = type(dict_structure[i][ss_structure[0]]) is dict
            if t == True:
                j = 0
                r = 0
                while j < len(ss_structure):
                    if type(dict_structure[i][ss_structure[j]]) is dict:
                        labels = list(dict_structure[i][ss_structure[j]].keys())
                        val = list(dict_structure[i][ss_structure[j]].values())
                        width = 0.35 + r
                        plt.bar(labels, val, width, label=str(ss_structure[j]))
                        plt.title(str(i))
                        plt.legend()
                    else:
                        label = ss_structure[j]
                        val = dict_structure[i][ss_structure[j]]
                        width = 0.35
                        plt.bar(label, val, width, label=str(ss_structure[j]))
                        plt.title(str(i))
                        plt.legend()
                    j = j + 1
                    r = r - 0.1
                plt.show()

            else:
                labels = ss_structure
                val = list(dict_structure[i].values())
                width = 0.35
                fig, ax = plt.subplots()
                ax.bar(labels, val, width, label='layers')
                ax.set_title(str(i))
            plt.show()


def number_of_sweeps(cell_id):
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
    number = len(all_sweeps)
    return number


def cell_ephys_info(cell_id):  # return the cell electrophysiological info for a given id
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
    for i in all_cell_features:
        if i["specimen_id"] == cell_id:
            return i


def sweep_info(cell_id, sweep_nber):  # return data info of a given sweep of a given id specimen
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
    for i in data_set:
        if i['sweep_number'] == sweep_nber:
            return i


def plot_stimulus(cell_id, sweep):
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
    for j in all_sweeps:
        if j["sweep_number"] == sweep:
            if j["stimulus_name"] == "Test":
                print("no plot")
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





def get_structure_name(species_dict, speciment_index):
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

    full_name = species_dict[speciment_index]['structure__name']  # full_name=dict_mouse[i]["structure__name"]
    j = 0
    for e in full_name:
        if e == ",":
            j = j + 1
    if j == 0:
        name = full_name.replace('"', '')
        substructure = "No substructure"
        layer = "No layer"
    elif j == 1:
        name1, layer1 = full_name.split(",")
        name = name1.replace('"', '')
        layer = layer1.replace('"', '')
        substructure = "No substructure"

    elif j == 2:
        name1, substructure, layer1 = full_name.split(",")
        name = name1.replace('"', '')
        layer = layer1.replace('"', '')
    return name, substructure, layer


def create_species_sweeps_stim_table(species_dict, from_csv=False):
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
    if from_csv == True:
        full_dataframe = pd.read_csv(str(species_dict))
        full_dataframe = full_dataframe.iloc[:, 1:]

        for row in range(full_dataframe.shape[0]):
            for col in range(7, 20):
                if type(full_dataframe.iloc[row, col]) == str:
                    full_dataframe.iat[row, col] = literal_eval(full_dataframe.iat[row, col])

                else:
                    full_dataframe.iloc[row, col] = -1
        return full_dataframe

    start = time.time()
    full_dataframe = pd.DataFrame()

    for current_specimen in range(len(species_dict)):
        name, substructure, layer = get_structure_name(species_dict, current_specimen)
        number_of_sweeps = 0

        specimen_stimulus_dict = stimulus_type(species_dict[current_specimen]["specimen__id"])

        for key in stimulus_type(species_dict[current_specimen]["specimen__id"]).keys():
            number_of_sweeps += len(stimulus_type(species_dict[current_specimen]["specimen__id"])[key])

        current_row = pd.Series([species_dict[current_specimen]["specimen__id"],
                                 species_dict[current_specimen]["structure__acronym"],
                                 name,
                                 substructure,
                                 layer,
                                 species_dict[current_specimen]["tag__dendrite_type"],
                                 number_of_sweeps],
                                index=['specimen_id', 'structure_acronyme', 'structure_name', 'structure_substructure',
                                       'layer', 'dendrite_type', 'number_of_sweeps'])
        current_row = pd.DataFrame(current_row).T

        stim_dict = pd.Series(specimen_stimulus_dict)
        stim_dict = pd.DataFrame(stim_dict).T

        full_row = pd.concat([current_row, stim_dict], axis=1)

        full_dataframe = pd.concat([full_dataframe, full_row], ignore_index=True, axis=0, join='outer')
    end = time.time()
    print('Running time= ', end - start)
    return full_dataframe


def pie_plotly(full_species_dict, count_what, based_on):
    # fig=plotly.sunburst(mouse_full_table,path=['structure_name','layer','number of sweeps'])
    my_table = pd.DataFrame(columns=[str(based_on), str(count_what)])

    for current_elt in list(set(full_species_dict[based_on])):
        number = (full_species_dict.loc[full_species_dict[based_on] == current_elt, [count_what]].sum(axis=0))[0]

        new_line = pd.Series([current_elt, number], index=[str(based_on), str(count_what)])

        my_table = my_table.append(new_line, ignore_index=True)

    figu = plotly.pie(my_table, values=len(count_what), names=based_on)
    figu.update_traces(textposition='inside', textinfo='percent+label')
    figu.show()


def create_subtable(full_species_table, which_factor, which_values):
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
    stim_list = full_species_table.columns[7:]
    if len(which_factor) == 1:
        if which_factor[0][0] in stim_list:
            new_table = full_species_table[~full_species_table[which_factor[0][0]].isin(which_values[0])]
            return (new_table)
        else:
            new_table = full_species_table[full_species_table[which_factor[0][0]].isin(which_values[0])]
            return (new_table)
    else:

        if which_factor[0][0] in stim_list:
            new_table = full_species_table[~full_species_table[which_factor[0][0]].isin(which_values[0])]
            which_factor = which_factor[1:]
            which_values = which_values[1:]
            new_table = create_subtable(new_table, which_factor, which_values)
            return (new_table)

        else:
            new_table = full_species_table[full_species_table[which_factor[0][0]].isin(which_values[0])]
            which_factor = which_factor[1:]
            which_values = which_values[1:]
            new_table = create_subtable(new_table, which_factor, which_values)
            return (new_table)


def select_specimen_sweep_stimulus(species_table, stimulus):
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
    return species_table[["specimen_id", stimulus]]




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

#%%
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
    m, c = np.linalg.lstsq(A, y, rcond=None)[0]

    return m, c


def plot_f_I_curve(f_I_table, single_specimen=None):
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
        f_I_table = f_I_table[f_I_table['specimen'] == single_specimen]
    linear_fit_param = pd.DataFrame(columns=['specimen', 'slope', 'intercept'])
    myplot = ggplot()
    for specimen in f_I_table['specimen'].cat.categories:
        sub_specimen_table = f_I_table.loc[f_I_table["specimen"] == specimen]

        sub_table = sub_specimen_table.loc[sub_specimen_table["frequence_Hz"] > 0]

        slope, intercept = fit_specimen_fi_slope(np.array(sub_table["stim_amplitude_pA"]),
                                                 np.array(sub_table['frequence_Hz']))
        new_line = pd.Series([specimen, slope, intercept], index=['specimen', 'slope', 'intercept'])

        linear_fit_param = linear_fit_param.append(new_line, ignore_index=True)

    linear_fit_param['specimen'] = pd.Categorical(linear_fit_param['specimen'])
    myplot = ggplot(f_I_table,
                    aes(x="stim_amplitude_pA", y="frequence_Hz", color="specimen")) + geom_point() + geom_abline(
        linear_fit_param, aes(intercept="intercept", slope="slope", colour="specimen"))

    return (myplot, linear_fit_param)


# my_sweep=my_specimen_data.get_sweep(sweep_nb)

def random_id_list(id_list, n):
    max_value = len(id_list)
    new_id_list = []
    for elt in range(n):
        index = random.randint(1, max_value - 1)

        new_id_list.append(id_list[index])
    return new_id_list


def ephys_spike_key_extractor(cell_id, sweep, t_start,
                              t_end):  # return the keys/feature names caracterizing the spike ; if don't want to put t_start/end, make them equal to None
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
        index_range = sweep_data["index_range"]
        i = sweep_data["stimulus"][0:index_range[1] + 1]  # in A = values of the stimulus
        v = sweep_data["response"][0:index_range[1] + 1]  # in V = values of the response
        i *= 1e12  # conversion to pA
        v *= 1e3  # conversion to mV

        sampling_rate = sweep_data["sampling_rate"]  # in Hz
        t = np.arange(0, len(v)) * (1.0 / sampling_rate)  # scale of the time ?

        sweep_ext = EphysSweepFeatureExtractor(t=t, v=v, i=i, start=t_start, end=t_end)
        sweep_ext.process_spikes()
        liste = sweep_ext.spike_feature_keys()
        return (liste)


def spike_info(cell_id, sweep_number,key):
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
    info = sweep_ext.spike_feature(key)
    return info


def norm_diff(a):
    """Calculate average of (a[i] - a[i+1]) / (a[i] + a[i+1])."""

    if len(a) <= 1:
        return np.nan

    a = a.astype(float)
    if np.allclose((a[1:] + a[:-1]), 0.):
        return 0.
    norm_diffs = (a[1:] - a[:-1]) / (a[1:] + a[:-1])

    norm_diffs[(a[1:] == 0) & (a[:-1] == 0)] = 0.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning, module="numpy")
        avg = np.nanmean(norm_diffs)
    return avg

def adaptation_index(isis):
    if len(isis) == 0:
        return np.nan

    return norm_diff(isis)

def get_cell_stim_response_table(cell_id,sweep_nb):
    '''
    Create table for a given cell and sweep number with stimulus, response, time arrays, spikes times, sampling rate and input resistance 

    Parameters
    ----------
    cell_id : int
        Cell id.
    sweep_nb : int
        Sweep number.

    Returns
    -------
    Cell_dict:dict
        Dictionnary containing the stimulus array, the response array, the time array, the spike time array, the sampling rate and the input resistance

    '''
    current_data=ctc.get_ephys_data(cell_id).get_sweep(sweep_nb)
    
    index_range=current_data["index_range"]
    
    response=np.array(current_data["response"][0:index_range[1]+1])* 1e12  # to pA
    stimulus=np.array(current_data["stimulus"][0:index_range[1]+1])* 1e3# to mV
    sampling_rate=current_data["sampling_rate"]
    time=np.arange(0, len(response)) * (1.0 / sampling_rate)
    cell_dict={"stimulus":stimulus,
               "response":response,
               "time":time,
               "spike_times":ctc.get_ephys_data(cell_id).get_spike_times(sweep_nb),
               "sampling_rate":current_data["sampling_rate"],
               "input_resistance_mohm":cell_ephys_info(cell_id)['input_resistance_mohm']}
    return(cell_dict)
#%%
def extract_inst_freq_table(specimen_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_spike_nb=False,first_nth_spikes=0):
    '''
    Compute the instananous frequency in each interspike interval per sweep for a cell

    Parameters
    ----------
    specimen_id : int
        specimencell id.
    species_sweep_stim_table : DataFrame
        Coming from create_species_sweeps_stim_table function.

    Returns
    -------
    inst_freq_table: DataFrame
        Table containing for a given cell for each sweep the stimulus amplitude and the instantanous frequency per interspike interval.

    '''
    index_stim = species_sweep_stim_table.columns.get_loc('Long Square')
    index_specimen = species_sweep_stim_table.index[species_sweep_stim_table["specimen_id"] == specimen_id][0]
    
    my_specimen_data = ctc.get_ephys_data(specimen_id)
    sweep_numbers = species_sweep_stim_table.iloc[index_specimen, index_stim]
    maximum_nb_interval =0
 
    for current_sweep in sweep_numbers:
        if per_time==True:
            index_range=my_specimen_data.get_sweep(current_sweep)["index_range"]
            sampling_rate=my_specimen_data.get_sweep(current_sweep)["sampling_rate"]
            current_stim_array=(my_specimen_data.get_sweep(current_sweep)["stimulus"][0:index_range[1]+1])* 1e12 #to pA
            stim_start_index=index_range[0]+next(x for x, val in enumerate(current_stim_array[index_range[0]:]) if val != 0 )
            current_time_array=np.arange(0, len(current_stim_array)) * (1.0 / sampling_rate)
            stim_start_time=current_time_array[stim_start_index]
            end_time=stim_start_time+(first_x_ms*1e-3)
            
            spike_times=my_specimen_data.get_spike_times(current_sweep)[my_specimen_data.get_spike_times(current_sweep) <= end_time ]
            if len(spike_times)>maximum_nb_interval:
                maximum_nb_interval=len(spike_times)
        elif per_spike_nb==True:
            spike_times=my_specimen_data.get_spike_times(current_sweep)[:first_nth_spikes]
            if len(spike_times)>maximum_nb_interval:
                maximum_nb_interval=len(spike_times)
        else:
            spike_times=my_specimen_data.get_spike_times(current_sweep)
            if len(spike_times)>maximum_nb_interval:
                maximum_nb_interval=len(spike_times)
       
    mycolumns=["specimen","sweep","stim_amplitude_pA"]+["interval_"+str(i) for i in range(1,(maximum_nb_interval))]
    inst_freq_table=pd.DataFrame(index=np.arange(len(sweep_numbers)),columns=mycolumns)
    for col in range(inst_freq_table.shape[1]):
        inst_freq_table.iloc[:,col]=np.nan
        
    for line in range(len(sweep_numbers)):
        current_sweep=sweep_numbers[line]
        stim_amplitude=my_specimen_data.get_sweep_metadata(current_sweep)['aibs_stimulus_amplitude_pa']
        if per_time==True:
            index_range=my_specimen_data.get_sweep(current_sweep)["index_range"]
            sampling_rate=my_specimen_data.get_sweep(current_sweep)["sampling_rate"]
            current_stim_array=(my_specimen_data.get_sweep(current_sweep)["stimulus"][0:index_range[1]+1])* 1e12 #to pA
            stim_start_index=index_range[0]+next(x for x, val in enumerate(current_stim_array[index_range[0]:]) if val != 0 )
            current_time_array=np.arange(0, len(current_stim_array)) * (1.0 / sampling_rate)
            stim_start_time=current_time_array[stim_start_index]
            end_time=stim_start_time+(first_x_ms*1e-3)
            
            spike_times=my_specimen_data.get_spike_times(current_sweep)[my_specimen_data.get_spike_times(current_sweep) <= end_time ]
           
        elif per_spike_nb==True:
            spike_times=my_specimen_data.get_spike_times(current_sweep)[:first_nth_spikes]
            
        else:
            spike_times=my_specimen_data.get_spike_times(current_sweep)
            
        
        inst_freq_table.iloc[line,0]=specimen_id
        inst_freq_table.iloc[line,1]=current_sweep   
        inst_freq_table.iloc[line,2]=stim_amplitude

        if len(spike_times) >7:
            for current_spike_time_index in range(1,len(spike_times)):
                current_inst_frequency=1/(spike_times[current_spike_time_index]-spike_times[current_spike_time_index-1])
                
                inst_freq_table.iloc[line,(current_spike_time_index+2)]=current_inst_frequency
        
            inst_freq_table.iloc[line,3:]/=max(inst_freq_table.iloc[line,3:])
    inst_freq_table = inst_freq_table.sort_values(by=["specimen", 'stim_amplitude_pA'])
    inst_freq_table['specimen']=pd.Categorical(inst_freq_table['specimen'])
    
    return(inst_freq_table)

def table_to_fit(inst_freq_table):
    '''
    Create table of interspike index-instantanous frequency- stimulus amplitude

    Parameters
    ----------
    inst_freq_table : DataFrame
        Coming from extract_inst_freq_table.

    Returns
    -------
    interval_freq_table : DataFrame
        Reorganized table to fit exponential function in fit_exponential_decay function.

    '''
    interval_freq_table=pd.DataFrame(columns=['interval','inst_frequency','stimulus_amplitude'])
    isnull_table=inst_freq_table.isnull()
    for col in range(3,(inst_freq_table.shape[1])):
        for line in range(inst_freq_table.shape[0]):
            if isnull_table.iloc[line,col] == False:
                new_line=pd.Series([int(col-2),inst_freq_table.iloc[line,col],inst_freq_table.iloc[line,2]],
                                   index=['interval','inst_frequency','stimulus_amplitude'])
                interval_freq_table=interval_freq_table.append(new_line,ignore_index=True)
   
    specimen=pd.Series(np.array([inst_freq_table.iloc[0,0]]*interval_freq_table.shape[0]))
    interval_freq_table=pd.concat([specimen,interval_freq_table],axis=1)
    interval_freq_table.columns=["specimen",'interval','inst_frequency','stimulus_amplitude']
    interval_freq_table['specimen']=pd.Categorical(interval_freq_table['specimen'])

    return interval_freq_table

def my_exponential_decay(x,Y0,tau,b):
    '''
    Parameters
    ----------
    x : Array
        interspike interval index array.
    Y0 : flt
        initial instantanous frequency .
    tau : flt
        Adaptation index constant.
    b : flt
        intantaneous frequency limit.

    Returns
    -------
    y : array
        Modelled instantanous frequency.

    '''
    y=Y0*np.exp(-tau*x)+b
    
    return y

def fit_exponential_decay(interval_freq_table):
    '''
    Parameters
    ----------
    interval_freq_table : DataFrame
        Comming from table_to_fit function.

    Returns
    -------
    my_plot : ggplot
        
    starting_freq : flt
        estimated initial instantanous frequency.
    adapt_cst : flt
        Adaptation index constant.
    limit_freq : flt
        intantaneous frequency limit.
    pcov_overall : 2-D array
        The estimated covariance of popt

    '''
    
    try:
        
        x_data=interval_freq_table.iloc[:,1]
        y_data=interval_freq_table.iloc[:,2]
        
        initial_amount=np.mean(interval_freq_table[interval_freq_table['interval']==1]['inst_frequency'])
        initial_tau=1
        initial_limit=np.mean(interval_freq_table[interval_freq_table['interval']==max(interval_freq_table['interval'])]['inst_frequency'])
        initial_estimate=[initial_amount,initial_tau,initial_limit]
        parameters_boundaries=([0,0,0],[max(interval_freq_table[interval_freq_table['interval']==1]['inst_frequency']),np.inf,max(interval_freq_table[interval_freq_table['interval']==1]['inst_frequency'])])

    
        popt_overall,pcov_overall=curve_fit(my_exponential_decay,x_data,y_data,p0=initial_estimate,bounds=parameters_boundaries,check_finite=False)
        
        sim_interval=np.arange(1,(max(interval_freq_table['interval'])+1))
        
        sim_freq=np.array(my_exponential_decay(sim_interval, *popt_overall))
        sim_interval=pd.Series(sim_interval)
        sim_freq=pd.Series(sim_freq)
        specimen=pd.Series(np.array([interval_freq_table.iloc[0,0]]*len(sim_freq)))
        new_data=pd.concat([specimen,sim_interval,sim_freq],axis=1)
        new_data.columns=['specimen','interval','inst_freq']
    
        starting_freq,adapt_cst,limit_freq=popt_overall
        starting_freq=my_exponential_decay(1.0,*popt_overall)
        limit_freq=my_exponential_decay(max(interval_freq_table['interval']),*popt_overall)
        starting_freq,adapt_cst,limit_freq=round(starting_freq,2),round(adapt_cst,2),round(limit_freq,2)
        
        steady_state_frequency=limit_freq/starting_freq
        my_plot=np.nan
        # my_plot=ggplot(interval_freq_table,aes(x=interval_freq_table.columns[1],y=interval_freq_table.columns[2],color="stimulus_amplitude"))+geom_point()+geom_line(new_data,aes(x=new_data.columns[1],y=new_data.columns[2]),color='black')
        # my_plot= my_plot+geom_text(x=(max(interval_freq_table['interval'])+1)/2,y=max(interval_freq_table[interval_freq_table['interval']==1]['inst_frequency']),
        #                            label="Tau ="+str(round(adapt_cst,2))+"+/-"+str(round(pcov_overall[1,1],2))+
        #                            ' , f_in='+str(round(starting_freq,2))+'+/-'+str(round(pcov_overall[0,0],2))+
        #                            ', f_lim='+str(round(limit_freq,2))+'+/-'+str(round(pcov_overall[2,2],2)),color="black",size=10)
        
        return my_plot,starting_freq,adapt_cst,limit_freq,steady_state_frequency,pcov_overall,popt_overall
    except (StopIteration):
        print("Stopped Iteration")
        my_plot=np.nan
        starting_freq=np.nan
        adapt_cst=np.nan
        limit_freq=np.nan
        steady_state_frequency=np.nan
        pcov_overall=np.nan
        popt_overall=np.nan
        return my_plot,starting_freq,adapt_cst,limit_freq,steady_state_frequency,pcov_overall,popt_overall
    except (ValueError):
        print("stopped_valueError")
        my_plot=np.nan
        starting_freq=np.nan
        adapt_cst=np.nan
        limit_freq=np.nan
        steady_state_frequency=np.nan
        pcov_overall=np.nan
        popt_overall=np.nan
        return my_plot,starting_freq,adapt_cst,limit_freq,steady_state_frequency,pcov_overall,popt_overall
    except(RuntimeError):
        print("Can't fit exponential, least-square optimization failed")
        my_plot=np.nan
        starting_freq=np.nan
        adapt_cst=np.nan
        limit_freq=np.nan
        steady_state_frequency=np.nan
        pcov_overall=np.nan
        popt_overall=np.nan
        return my_plot,starting_freq,adapt_cst,limit_freq,steady_state_frequency,pcov_overall,popt_overall

def goodness_of_fit(y_data,simulated_y_data):
    '''
    Compute determination coefficient R2

    Parameters
    ----------
    y_data : Array
        Actual values.
    simulated_y_data : Array
        Estimated value.

    Returns
    -------
    r2 : flt
        Goodness of fit.

    '''
    # residual sum of squares
    ss_res = np.sum((y_data - simulated_y_data) ** 2)
    
    # total sum of squares
    ss_tot = np.sum((y_data - np.mean(y_data)) ** 2)
    
    # r-squared
    r2 = 1 - (ss_res / ss_tot)
    r2=round(r2,2)
    return r2

def extract_ISI(spike_times):
    '''
    Extract the ISI from the spike times array

    Parameters
    ----------
    spike_times : array
        Numpy array containing the spike times

    Returns
    -------
    Array
        Numpy array containing the ISI of a sweep.

    '''
    if len(spike_times) <= 1:
        return np.array([])
    return spike_times[1:] - spike_times[:-1]
#%%
def extract_stim_freq(specimen_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0):
    '''
    Function to extract for each specified specimen_id and the corresponding stimulus the frequency of the response
    Frequency is defined as the number of spikes divided by the time between the stimulus start and the time of the specified index
    Parameters
    ----------
    specimen_id : int

    
    Returns
    -------
    f_I_table : DataFrame
        DataFrame with a column "specimen_id"(factor),the sweep number (int),the stimulus amplitude in pA(float),and the computed frequency of the response (float).
    '''
    print(specimen_id)
    # import os
    # is_file_ok=False
    # while is_file_ok ==False:
    #     try:
    #         index_stim = species_sweep_stim_table.columns.get_loc('Long Square')
    #         index_specimen = species_sweep_stim_table.index[species_sweep_stim_table["specimen_id"] == specimen_id][0]
            
    #         my_specimen_data = ctc.get_ephys_data(specimen_id)
    #         sweep_numbers = species_sweep_stim_table.iloc[index_specimen, index_stim]
    #         index_range=my_specimen_data.get_sweep(sweep_numbers[0])["index_range"]
    #         is_file_ok=True
    #         print("ok")
    #     except (OSError):
    #         try:
    #             os.remove(str("/Users/julienballbe/My_Work/Allen_Data/Common_Script/Full_analysis_cell_types/specimen_"+str(specimen_id)+"/ephys.nwb"))
    #         except(FileNotFoundError):
    #             os.remove(str("/Users/julienballbe/My_Work/Allen_Data/Common_Script/Full_analysis_cell_types/specimen_"+str(specimen_id)+"/ephys_sweeps.json"))
    #             os.rmdir("/Users/julienballbe/My_Work/Allen_Data/Common_Script/Full_analysis_cell_types/specimen_"+str(specimen_id))
    #         print("Truncated file for specimen ",str(specimen_id) ,". File removed to be redownloaded")
    #         continue

    
            
    # print(specimen_id)
    f_I_table = pd.DataFrame(columns=['specimen', 'sweep', 'stim_amplitude_pA', 'frequence_Hz'])
    index_stim = species_sweep_stim_table.columns.get_loc('Long Square')
    index_specimen = species_sweep_stim_table.index[species_sweep_stim_table["specimen_id"] == specimen_id][0]
    
    my_specimen_data = ctc.get_ephys_data(specimen_id)
    sweep_numbers = species_sweep_stim_table.iloc[index_specimen, index_stim]
    
    
    for current_sweep in sweep_numbers:
        
        index_range=my_specimen_data.get_sweep(current_sweep)["index_range"]
        
        sampling_rate=my_specimen_data.get_sweep(current_sweep)["sampling_rate"]
        current_stim_array=(my_specimen_data.get_sweep(current_sweep)["stimulus"][0:index_range[1]+1])* 1e12 #to pA
    
        
        stim_start_index=index_range[0]+next(x for x, val in enumerate(current_stim_array[index_range[0]:]) if val != 0 )
        current_time_array=np.arange(0, len(current_stim_array)) * (1.0 / sampling_rate)
        
        stim_start_time=current_time_array[stim_start_index]
       
       
        
        if len(my_specimen_data.get_spike_times(current_sweep)) <2:
            freq = 0

        else :
            if per_nth_spike==True:
                reshaped_spike_times=my_specimen_data.get_spike_times(current_sweep)[:first_nth_spike]
    
                nb_spike = len(reshaped_spike_times)
                t_last_spike = reshaped_spike_times[-1]
                freq = nb_spike / (t_last_spike - stim_start_time)
                
            elif per_time==True:
                end_time=stim_start_time+(first_x_ms*1e-3)

                reshaped_spike_times=my_specimen_data.get_spike_times(current_sweep)[my_specimen_data.get_spike_times(current_sweep) <= end_time ]
                nb_spike = len(reshaped_spike_times)
                if nb_spike !=0:
                    t_last_spike = reshaped_spike_times[-1]
                    freq = nb_spike / (t_last_spike - stim_start_time)
                else:
                    freq=0
        new_line = pd.Series([int(specimen_id), current_sweep,
                              my_specimen_data.get_sweep_metadata(current_sweep)['aibs_stimulus_amplitude_pa'],
                              freq],
                             index=['specimen', 'sweep', 'stim_amplitude_pA', 'frequence_Hz'])
        f_I_table = f_I_table.append(new_line, ignore_index=True)
    
    f_I_table = f_I_table.sort_values(by=["specimen", 'stim_amplitude_pA'])
    f_I_table['specimen'] = pd.Categorical(f_I_table['specimen'])
    return f_I_table
    

def mysigmoid(x,maxi,x0,slope):
    y=maxi/(1+np.exp((x0-x)/slope))
    return y

def fit_sigmoid(f_I_table):
    '''
    Fit a sigmoid curve to the stim/amplitude data points to extract several I/O metrcis : the threshold, the saturation and the gain

    Parameters
    ----------
    f_I_table : DataFrame
        Stimulus_frequency table for one cell.

    Returns
    -------
    estimated_gain : float
        estimated gain of the I/O.
    estimated_threshold : float
        estimated neuron threshold.
    estimated_saturation : float
        estimated neuron saturation firing rate.
    my_plot : ggplot
        plot of the data point with the sigmoid fit and the linear fit to the linear part of the sigmoid.
    pcov: 2-D array
        The estimated covariance of popt
    popt: 1D array
        Estimated parameters of function fit
    
    '''
    try:
        x_data=f_I_table.iloc[:,2]
        y_data=f_I_table.iloc[:,3]
        
        
        ##Get the initial estimate for the fit of sigmoid
        #Get the maximum firing rate of the data
        maxi=max(x_data)
        
        #Get the index corresponding to the median non-zero firing rate
        

        without_zero_index=next(x for x, val in enumerate(y_data) if val >0 )
        median_firing_rate_index=next(x for x, val in enumerate(y_data) if val >= np.median(y_data.iloc[without_zero_index:]))
        #Get the stimulus amplitude correspondingto the median non-zero firing rate
    
        x0=x_data.iloc[median_firing_rate_index]
    
        #Get the slope from the linear fit of the firing rate
        slope=fit_specimen_fi_slope(x_data,y_data)[0]
    
    
        initial_estimate=[maxi,x0,slope]
        parameters_boundaries=([0,0,0],[np.inf,np.inf,np.inf])
        

        popt,pcov=curve_fit(mysigmoid,x_data,y_data,p0=initial_estimate,bounds=parameters_boundaries,check_finite=False)
        
        new_x_data=pd.Series(np.arange(min(x_data),max(x_data),1))
        new_y_data=pd.Series(mysigmoid(new_x_data,*popt))
        new_data=pd.concat([new_x_data,new_y_data],axis=1,ignore_index=True)
        new_data.columns=["stim_amplitude_pA","frequence_Hz"]
        slope_confidence_threshold=4
        if trust_sigmoid(new_x_data, *popt, slope_confidence_threshold)==False:
            estimated_gain=np.nan

            estimated_saturation=np.nan
            estimated_threshold=np.nan
            my_plot=np.nan
            pcov=np.nan
            popt=np.nan
            return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
        else:

            #get index 25% and 75% of max firing rate
            twentyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.25*popt[0]))
            seventyfive_index=next(x for x, val in enumerate(new_y_data) if val >(0.75*popt[0]))
            #fit linear line to linear sigmoid portion
            linear_estimated_slope,linear_estimated_intercept=fit_specimen_fi_slope(new_x_data.iloc[twentyfive_index:seventyfive_index],mysigmoid(new_x_data.iloc[twentyfive_index:seventyfive_index],*popt))
            estimated_threshold=(0-linear_estimated_intercept)/linear_estimated_slope
            my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(popt[0],popt[1],popt[2])))
           
            if my_derivative[-1]<0.001:
                
                estimated_saturation=popt[0]
            else:
                estimated_saturation=np.nan
            estimated_gain=linear_estimated_slope
            my_plot=np.nan
            # my_plot=ggplot(f_I_table,aes(x=f_I_table.columns[2],y=f_I_table.columns[3]))+geom_point()+geom_line(new_data,aes(x=new_data.columns[0],y=new_data.columns[1]),color='blue')+geom_abline(aes(intercept=linear_estimated_intercept,slope=linear_estimated_slope))
            # my_plot+=geom_text(x=10,y=estimated_saturation,label="Gain="+str(round(estimated_gain,2))+
            #                    ', thresh='+str(round(estimated_threshold,2))+
            #                    ', sat='+str(round(estimated_saturation,2))+'+/-'+str(round(pcov[0,0],2)),size=10,color="black")
            return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (StopIteration):
        print("Stopped Iteration")
        print(f_I_table.iloc[1,1])
        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (ValueError):
        print("stopped_valueError")
        print(f_I_table.iloc[1,1])
        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    except (RuntimeError):
        print("Can't fit sigmoid, least-square optimization failed")
        print(f_I_table.iloc[1,1])
        estimated_gain=np.nan
        estimated_saturation=np.nan
        estimated_threshold=np.nan
        my_plot=np.nan
        pcov=np.nan
        popt=np.nan
        return estimated_gain,estimated_threshold,estimated_saturation,my_plot,pcov,popt
    

#%%General function to extract parameters
def extract_feature(specimen_id,species_sweep_stim_table,per_time=False,first_x_ms=0,per_nth_spike=False,first_nth_spike=0):
    '''
    

    Parameters
    ----------
    specimen_id : list
        list of specimen id.
    species_sweep_stim_table : DataFrame
        DataFrame coming from create_species_sweeps_stim_table function.
    per_time : Boolean, optional
        Extract feature per time . The default is False.
    first_x_ms : int, optional
        If per_time==True , indicate the time to take into account to extract features in ms (x ms after the start of the stimulus). The default is 0.
    per_nth_spike : Boolean, optional
        Extract feature per spike number. The default is False.
    first_nth_spike : int, optional
        If per_nth_spike==True, indicate the n first spike to take into account for feature extraction. The default is 0.

    Returns
    -------
    full_table : DataFrame
        DataFrame containing the different feature values for each cell .

    '''
    mycolumns=['Cell_id','Gain','Threshold','Saturation','Adapt_cst','Starting_frequency','Limit_frequency','Steady_state_frequency']
    full_table=pd.DataFrame(columns=mycolumns)

    for current_specimen_id in specimen_id:

        f_I_table=extract_stim_freq(current_specimen_id,species_sweep_stim_table,per_time,first_x_ms,per_nth_spike,first_nth_spike)
        f_I_gain,f_I_threshold,f_I_saturation,f_I_plot,f_I_pcov,f_I_popt=fit_sigmoid(f_I_table)
        
        inst_freq_table=table_to_fit(extract_inst_freq_table(current_specimen_id,species_sweep_stim_table,per_time,first_x_ms,per_spike_nb=per_nth_spike,first_nth_spikes=first_nth_spike))
        inst_plot,inst_starting_freq,inst_adapt_cst,inst_limit_freq,inst_steady_state_frequency,inst_pcov,inst_popt=fit_exponential_decay(inst_freq_table)
        
        
        new_line=pd.Series([int(current_specimen_id),
                            round(f_I_gain,2),
                            round(f_I_threshold,2),
                            round(f_I_saturation,2),
                            round(inst_adapt_cst,2),
                            round(inst_starting_freq,2),
                            round(inst_limit_freq,2),
                            round(inst_steady_state_frequency,2)],
                           index=mycolumns)
        
        full_table=full_table.append(new_line, ignore_index=True)
    full_table['Cell_id']=pd.Categorical(full_table['Cell_id'])
    return full_table
#%%


#%%My tables
# mycol=['cell_id',
#        'Gain',
#        'Thresold',
#        'Saturation',
#        'Adapt_cst',
#        'Init_freq',
#        'Limit_freq',
#        'Steady_state_freq']
# units=['--',
#        'Hz/pA',
#        'pA',
#        'Hz',
#        'spike_index',
#        'WU',
#        'WU',
#        'WU']
# first_line=pd.Series(mycol,index=mycol)
# second_line=pd.Series(units,index=mycol)
# first_two_lines=pd.DataFrame([first_line,second_line])

# main_start_time=time.time()
# start_time=time.time()

# features_5ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=5)
# features_5ms=features_5ms.reindex(columns=mycol)
# File_5ms=pd.concat([first_two_lines,features_5ms])
# File_5ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_5ms)+".csv"),na_rep="nan",index=False)
# end_5ms=time.time()
# print("time for 5ms=",end_5ms-start_time,"s")

# start_time=time.time()
# features_10ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=10)
# features_10ms=features_10ms.reindex(columns=mycol)
# File_10ms=pd.concat([first_two_lines,features_10ms])
# File_10ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_10ms)+".csv"),na_rep="nan",index=False)
# end_10ms=time.time()

# print("time for 10ms=",end_10ms-start_time,"s")
# start_time=time.time()
# features_25ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=25)
# features_25ms=features_25ms.reindex(columns=mycol)
# File_25ms=pd.concat([first_two_lines,features_25ms])
# File_25ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_25ms)+".csv"),na_rep="nan",index=False)
# end_25ms=time.time()
# print("time for 25ms=",end_25ms-start_time,"s")

# start_time=time.time()
# features_50ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=50)
# features_50ms=features_50ms.reindex(columns=mycol)
# File_50ms=pd.concat([first_two_lines,features_50ms])
# File_50ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_50ms)+".csv"),na_rep="nan",index=False)
# end_50ms=time.time()
# print("time for 50ms=",end_50ms-start_time,"s")

# start_time=time.time()
# features_100ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=100)
# features_100ms=features_100ms.reindex(columns=mycol)
# File_100ms=pd.concat([first_two_lines,features_100ms])
# File_100ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_100ms)+".csv"),na_rep="nan",index=False)
# end_100ms=time.time()
# print("time for 100ms=",end_100ms-start_time,"s")

# start_time=time.time()
# features_250ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=250)
# features_250ms=features_250ms.reindex(columns=mycol)
# File_250ms=pd.concat([first_two_lines,features_250ms])
# File_250ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_250ms)+".csv"),na_rep="nan",index=False)
# end_250ms=time.time()
# print("time for 250ms=",end_250ms-start_time,"s")

# start_time=time.time()
# features_500ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=500)
# features_500ms=features_500ms.reindex(columns=mycol)
# File_500ms=pd.concat([first_two_lines,features_500ms])
# File_500ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_500ms)+".csv"),na_rep="nan",index=False)
# end_500ms=time.time()
# print("time for 500ms=",end_500ms-start_time,"s")

# start_time=time.time()
# features_1000ms=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_time=True,first_x_ms=1000)
# features_1000ms=features_1000ms.reindex(columns=mycol)
# File_1000ms=pd.concat([first_two_lines,features_1000ms])
# File_1000ms.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_1000ms)+".csv"),na_rep="nan",index=False)
# end_1000ms=time.time()
# print("time for 1000ms=",end_1000ms-start_time,"s")

# start_time=time.time()
# features_4spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=4)
# features_4spikes=features_4spikes.reindex(columns=mycol)
# File_4spikes=pd.concat([first_two_lines,features_4spikes])
# File_4spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_4spikes)+".csv"),na_rep="nan",index=False)
# end_4spike=time.time()
# print("time for 4spikes=",end_4spike-start_time,"s")

# start_time=time.time()
# features_5spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=5)
# features_5spikes=features_5spikes.reindex(columns=mycol)
# File_5spikes=pd.concat([first_two_lines,features_5spikes])
# File_5spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_5spikes)+".csv"),na_rep="nan",index=False)
# end_5spike=time.time()
# print("time for 5spikes=",end_5spike-start_time,"s")

# start_time=time.time()
# features_6spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=6)
# features_6spikes=features_6spikes.reindex(columns=mycol)
# File_6spikes=pd.concat([first_two_lines,features_6spikes])
# File_6spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_6spikes)+".csv"),na_rep="nan",index=False)
# end_6spike=time.time()
# print("time for 6spikes=",end_6spike-start_time,"s")

# start_time=time.time()
# features_7spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=7)
# features_7spikes=features_7spikes.reindex(columns=mycol)
# File_7spikes=pd.concat([first_two_lines,features_7spikes])
# File_7spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_7spikes)+".csv"),na_rep="nan",index=False)
# end_7spike=time.time()
# print("time for 7spikes=",end_7spike-start_time,"s")

# start_time=time.time()
# features_8spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=8)
# features_8spikes=features_8spikes.reindex(columns=mycol)
# File_8spikes=pd.concat([first_two_lines,features_8spikes])
# File_8spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_8spikes)+".csv"),na_rep="nan",index=False)
# end_8spike=time.time()
# print("time for 8spikes=",end_8spike-start_time,"s")

# start_time=time.time()
# features_9spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=9)
# features_9spikes=features_9spikes.reindex(columns=mycol)
# File_9spikes=pd.concat([first_two_lines,features_9spikes])
# File_9spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_9spikes)+".csv"),na_rep="nan",index=False)
# end_9spike=time.time()
# print("time for 9spikes=",end_9spike-start_time,"s")

# start_time=time.time()
# features_10spikes=extract_feature(mouse_id_list,mouse_sweep_stim_table,per_nth_spike=True,first_nth_spike=10)
# features_10spikes=features_10spikes.reindex(columns=mycol)
# File_10spikes=pd.concat([first_two_lines,features_10spikes])
# File_10spikes.to_csv(path_or_buf=str("/Users/julienballbe/My_Work/Allen_Data/"+str(features_10spikes)+".csv"),na_rep="nan",index=False)
# end_10spike=time.time()
# print("time for 10spikes=",end_10spike-start_time,"s")
# start_time=time.time()

# main_end_time=time.time()
# print("Total time=",main_end_time-main_start_time,'s')
    
#%%Check slope of sigmoid
from scipy.misc import derivative
def trust_sigmoid(new_x_data,maxi,x0,slope,slope_confidence_threshold):
    
    my_derivative=np.array(derivative(mysigmoid,new_x_data,dx=1e-1,args=(maxi,x0,slope)))
   
    if max(my_derivative)>slope_confidence_threshold:
        print("slope to high")
        return False
    elif my_derivative[-1]==0:
        print("end slope too low")
        return False
    
    
    return True
# Morphology :

def spiny_state(dict_species):
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

    spiny = []
    aspiny = []
    sparsely_spiny = []
    for i in dict_species:
        if i['tag__dendrite_type'] == 'spiny':
            spiny.append(i['specimen__id'])
        else:
            if i['tag__dendrite_type'] == 'aspiny':
                aspiny.append(i['specimen__id'])
            else:
                sparsely_spiny.append(i['specimen__id'])
    return (spiny, aspiny, sparsely_spiny)


def reconstruction(dict_species):
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

    no = []
    yes = []
    for i in dict_species:
        if i['nr__reconstruction_type'] == None :
            no.append(i["specimen__id"])
        else:
            yes.append(i["specimen__id"])
    return (no, yes)


def spiny_reconstruction(dict_species):
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

    spiny_no_reconstruction = []
    spiny_reconstruction = []
    aspiny_no_reconstruction = []
    aspiny_reconstruction = []
    sparsely_no_reconstruction = []
    sparsely_reconstruction = []
    no, yes = reconstruction(dict_species)
    spiny, aspiny, sparsely = spiny_state(dict_species)
    for i in yes:
        if i in spiny:
            spiny_reconstruction.append(i)
        else:
            if i in aspiny:
                aspiny_reconstruction.append(i)
            else:
                sparsely_reconstruction.append(i)
    for i in no:
        if i in spiny:
            spiny_no_reconstruction.append(i)
        else:
            if i in aspiny:
                aspiny_no_reconstruction.append(i)
            else:
                sparsely_no_reconstruction.append(i)
    return (spiny_no_reconstruction, spiny_reconstruction, aspiny_no_reconstruction, aspiny_reconstruction,
            sparsely_no_reconstruction, sparsely_reconstruction)


def morpho_info(cell_id):  # for a given id, returns the morpho info of the cell
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


def swc_morpho(cell_id, name):  # return a link to download the morphology swc file of a given cell
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
    swc_file = ctc.get_reconstruction(cell_id, name)
    return (swc_file)


def compartment_nber(cell_id):
    """
    For a given cell it returns the number of all its compartments
    Parameters
    ----------
    cell_id : str
    Returns
    -------
    compartment_number : int
    """

    all_morpho = ctc.get_reconstruction(cell_id)
    compartment_number = len(all_morpho.compartment_list)
    return (compartment_number)


def twoD_morpho(cell_id):
    """
    For a given cell it plots its morphology in different spaces (but always in 2D)
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
        for c in morphology.children_of(n):  # print the different children compartements
            if n['type'] == 2:
                axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='red')
                axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='red')
            else:
                if n['type'] == 1:  # soma
                    axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='green')
                    axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='green')
                else:
                    axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='black')
                    axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='black')

    # cut dendrite markers
    dm = [m for m in markers if m['name'] == Marker.CUT_DENDRITE]  # if m['name'] ==10

    axes[0].scatter([m['x'] for m in dm], [m['y'] for m in dm], color='#3333ff')
    axes[1].scatter([m['z'] for m in dm], [m['y'] for m in dm], color='#3333ff')

    # no reconstruction markers
    nm = [m for m in markers if m['name'] == Marker.NO_RECONSTRUCTION]  # if m['name'] ==20

    axes[0].scatter([m['x'] for m in nm], [m['y'] for m in nm], color='#333333')
    axes[1].scatter([m['z'] for m in nm], [m['y'] for m in nm], color='#333333')

    axes[0].set_ylabel('y')
    axes[0].set_xlabel('x')
    axes[1].set_xlabel('z')
    red_patch = mpatches.Patch(color='red', label='axons')
    black_patch = mpatches.Patch(color='black', label='dendrites')
    green_patch = mpatches.Patch(color='green', label='soma')
    blue_patch = mpatches.Patch(color='blue', label='location of the dendrite truncations')
    grey_patch = mpatches.Patch(color='grey', label='location of the axon truncations')
    axes[1].legend(handles=[red_patch, black_patch, green_patch, blue_patch, grey_patch],
                   bbox_to_anchor=(0, 1.06, 1, 0.2))
    plt.show()

def twoD_morpho_without_truncation (cell_id) :
    """
    For a given cell it plots its morphology in different spaces (but always in 2D), but without the indications of axon/dendrite truncations
    Parameters
    ----------
    cell_id : str

    Returns
    -------
    None
    """
    morphology = ctc.get_reconstruction(cell_id)
    fig, axes = plt.subplots(1, 2, sharey=True, sharex=True)
    axes[0].set_aspect('equal', 'box')
    axes[1].set_aspect('equal', 'box')

    # Make a line drawing of x-y and y-z views
    for n in morphology.compartment_list:
        for c in morphology.children_of(n):  # print the different children compartments
            if n['type'] == 2:
                axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='red')
                axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='red')
            else:
                if n['type'] == 1:  # soma
                    axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='green')
                    axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='green')
                else:
                    axes[0].plot([n['x'], c['x']], [n['y'], c['y']], color='black')
                    axes[1].plot([n['z'], c['z']], [n['y'], c['y']], color='black')
    axes[0].set_ylabel('y')
    axes[0].set_xlabel('x')
    axes[1].set_xlabel('z')
    red_patch = mpatches.Patch(color='red', label='axons')
    black_patch = mpatches.Patch(color='black', label='dendrites')
    green_patch = mpatches.Patch(color='green', label='soma')
    axes[1].legend(handles=[red_patch, black_patch, green_patch],bbox_to_anchor=(0, 1.06, 1, 0.2))
    plt.show()

def id_reconstruction_given_structure (dic_species,structure_parent_acronym, reconstructed_yes_or_no) :
    """
    For a given species, structure and if there's a reconstruction or not, the function returns the cell_id having those caracteristics
    Parameters
    ----------
    dic_species : Dictionary
        Dictionary of all the cells of a species

    structure_parent_acronym : str
        the acronym of the structure of interest (ex : VISp)

    reconstructed_yes_or_no : list
        List of ids for cell having a reconstruction or having not a reconstruction

    Returns
    -------
    structure_with_reconstruction : list
        List of ids being in the given structure and having or not a reconstruction
    """

    structure_with_reconstruction = []
    for i in reconstructed_yes_or_no :
        info = cell_info(i, dic_species)
        if info['structure_parent__acronym'] == structure_parent_acronym :
            structure_with_reconstruction.append(i)
    return(structure_with_reconstruction)

def sorted_layers (id_list,dic_species) :
    """
    For a structure, it returns for each layer its ids

    Parameters
    ----------
    id_list : list
        list of ids of a structure

    dic_species : Dictionary
        Dictionary of all the cells of a species

    Returns
    -------
    layer_1, layer_23, layer_4, layer_5, layer_6a, layer_6b : lists
        Lists of ids of each layer
    """

    layer_1 = []
    layer_23 = []
    layer_4 = []
    layer_5 = []
    layer_6a = []
    layer_6b = []
    for i in id_list:
        if cell_info(i, dic_species)['structure__layer'] == '6b':
            layer_6b.append(cell_info(i,dic_species))
        else:
            if cell_info(i, dic_species)['structure__layer'] == '6a':
                layer_6a.append(cell_info(i,dic_species))
            else:
                if cell_info(i, dic_species)['structure__layer'] == '5':
                    layer_5.append(cell_info(i,dic_species))
                else:
                    if cell_info(i, dic_species)['structure__layer'] == '4':
                        layer_4.append(cell_info(i,dic_species))
                    else:
                        if cell_info(i, dic_species)['structure__layer'] == '2/3':
                            layer_23.append(cell_info(i,dic_species))
                        else:
                            layer_1.append(cell_info(i,dic_species))
    return (layer_1,layer_23,layer_4,layer_5,layer_6a,layer_6b)


# ctc = CellTypesCache(manifest_file='cell_types/manifest.json') #create the manifest file

# cells = ctc.get_cells(file_name="all_cells",require_reconstruction=True)
# id_mouse,id_human=take_id("all_cells")
# dict_mouse,dict_human=dict_specimen("all_cells")
# a=structure(dict_human)
# print(a)
# plot_number_structures(a)
# a=stimulus_type(464212183)
# print(a)