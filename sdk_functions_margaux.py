import allensdk
import json
from allensdk.core.cell_types_cache import CellTypesCache
from allensdk.api.queries.cell_types_api import CellTypesApi
import webbrowser
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from allensdk.core.swc import Marker
import matplotlib.patches as mpatches
from allensdk.ephys.ephys_extractor import EphysSweepFeatureExtractor

ctc = CellTypesCache(manifest_file='cell_types/manifest.json') #create the manifest file
#this is a new test
cells = ctc.get_cells(file_name="all_cells",require_reconstruction=True)


#data=ctc.get_all_features(dataframe=True, require_reconstruction=True)

def dict_specimen(name_file): #return all the cell info for all cell and for each specimen ; name_file=all_cells
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

def cell_info (id_specimen, dict_specimen) : #return the cell info for a given id and type of specie
    for i in dict_specimen :
        if i["specimen__id"]==id_specimen :
            return i
def cell_marker (id_specimen, dict_specimen): #only works for mouse, the human species doesn't have a line name
    for i in dict_specimen :
        if i["specimen__id"]==id_specimen :
            return i["line_name"]


def all_cell_markers (dict_specimen) :
    dic_type=dict()
    for i in dict_specimen :
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

def markers (dict_specimen) : #return in a list the structures in str followed by their dictionary containing for each marker the layer and its id
    dic=all_cell_markers(dict_specimen)
    liste=[]
    for i in dic :
        liste.append(i)
        ss_str=list(dic[i].keys())
        #print(i)
        marker = list(marker_number(dict_specimen).keys())
        m = dict()
        for j in ss_str :
            m_or_l=list(dic[i][j].keys()) #les markers si pas de substructures, sinon c'est les layers de la substructure
            t=type(dic[i][j][m_or_l[0]]) is dict # pour savoir si c'est layer ou marker
            if t==False : #c'est des markers
                for k in m_or_l : #pour chaque marker
                    if k not in m :
                        m[k]=dict()
                        m[k][j]=dic[i][j][k]
                    else :
                        if j not in m[k] :
                            m[k][j] = dic[i][j][k]

            else :
                if t==True:
                    #print(j) #ss-structure
                    for d in m_or_l : #d est la layer
                        mm=list(dic[i][j][d].keys())
                        for l in mm :
                            if l not in m :
                                m[l]=dict()
                                m[l][j]=dict()
                                m[l][j][d]=dic[i][j][d][l]


        liste.append(m)
    return(liste)


def plot_marker (dict_specimen) : # within each structure of a given species, for each marker it plots an histogram of its proportion of spiny/aspiny/sparsely spiny with or without reconstruction in each layer

    l=markers(dict_specimen)
    for i in l :

        t=type(i) is dict
        if t==False :
            titre=i
        if t==True :
            #print(i)
            for j in i :
                aa=0
                t=type(j) is dict
                #print(j)
                if t==False:
                    ss_titre=j
                    layers=list(i[j].keys())
                    len_spi_no=[]
                    len_spi_yes=[]
                    len_asp_no=[]
                    len_asp_yes=[]
                    len_spa_no=[]
                    len_spa_yes=[]

                    spi_no = []
                    spi_yes = []
                    asp_no = []
                    asp_yes = []
                    spa_no = []
                    spa_yes = []
                    a, b, c, d, e, f = spiny_reconstruction(dict_specimen)
                    for k in layers :
                        f=type(i[j][k]) is dict
                        if f==False:
                            ss_title=" "
                            lab=layers
                            #print(lab)
                            ids = i[j][k]
                            # print(a)
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
                        else :
                            if f==True :
                                sub_structure=i[j][k]
                                ss_title=", "+str(k)+", "
                                for layer in i[j][k] :
                                    lab=list(i[j][k].keys())
                                    ids=i[j][k][layer]
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
                    #fig, ax = plt.subplots()
                    grd_titre=titre+ss_title+str(ss_titre)
                    plt.bar(lab,len_spi_no,label='spiny without reconstruction')
                    plt.bar(lab,len_spi_yes,label='spiny with reconstruction')
                    plt.bar(lab,len_asp_no,label='aspiny without reconstruction')
                    plt.bar(lab,len_asp_yes,label='aspiny with reconstruction')
                    plt.bar(lab,len_spa_no,label='sparsely spiny without reconstruction')
                    plt.bar(lab,len_spa_yes,label='sparsely spiny with reconstruction')
                    plt.title(grd_titre)
                    plt.legend()
                    plt.show()


#'Retrosplenial area' : au lieu de layer ça affiche le nom de la substructure et ça met 1 sparsely pour tous les markers
#regarder au niveau de la boucle for k in layer

def marker_number (dict_specimen) : #return a dictionary with the marker names and how many they are ; the "total" is only a verification that the dic_type as indeed all the id of the dict_specimen
    dico=all_cell_markers(dict_specimen)
    total = 0
    dico_marker = dict()
    for i in dico:
        t=type(dico[i]) is dict
        if t==False :
            if i not in dico_marker :
                dico_marker[i]=len(dico[i])
            else :
                dico_marker[i]=dico_marker[i]+len(dico[i])
        for j in dico[i]:
            t = type(dico[i][j]) is dict
            if t == False:
                total = total + len(dico[i][j])
                if j not in dico_marker:
                    dico_marker[j] = len(dico[i][j])
                else:
                    dico_marker[j] = dico_marker[j] + len(dico[i][j])
            else:
                for k in dico[i][j]:
                    t = type(dico[i][j][k]) is dict
                    if t == False:
                        total = total + len(dico[i][j][k])
                        if k not in dico_marker:
                            dico_marker[k] = len(dico[i][j][k])
                        else:
                            dico_marker[k] = dico_marker[k] + len(dico[i][j][k])

                    else:
                        for z in dico[i][j][k]:
                            t = type(dico[i][j][k][z]) is dict
                            if t == False:
                                total = total + len(dico[i][j][k][z])
                                if z not in dico_marker:
                                    dico_marker[z] = len(dico[i][j][k][z])
                                else:
                                    dico_marker[z] = dico_marker[z] + len(dico[i][j][k][z])
    return (dico_marker)


def plot_marker_number(dict_specimen) : #plot an histogram of the all the markers of the specimen and how many cells have those markers
    d=marker_number(dict_specimen)
    height=d.values()
    bars=d.keys()
    x_pos=np.arange(1,370,10)
    plt.bar(x_pos, height,width=5)
    plt.xticks(x_pos,bars,rotation=90,fontsize=10)
    specie=dict_specimen[0]["donor__species"]
    title="markers of the "+specie
    plt.tight_layout()
    plt.title(title)
    plt.show()

def take_id(name_file) : #return two lists : the mouse id list and human id list
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


def ephys_web_page(specimen_id) : #the specimen id web page in the Allen Brain Atlas is opened
    link="http://celltypes.brain-map.org/experiment/electrophysiology/"+str(specimen_id)
    webbrowser.open(link)

def morpho_web_page(specimen_id) : #the specimen id web page in the Allen Brain Atlas is opened
    link="http://celltypes.brain-map.org/experiment/morphology/"+str(specimen_id)
    webbrowser.open(link)



def structure (dico_specimen) : #return all the structures, as well as its substructures and layers, of a given specimen as a dictionary ; need to have run the dict_specimen in first place
    i=0
    liste=[]
    dico_structure=dict()
    while i<len(dico_specimen) : #dict_mouse
        full_name=dico_specimen[i]['structure__name'] #full_name=dict_mouse[i]["structure__name"]
        j=0
        for e in full_name :
            if e=="," :
                j=j+1
        if j==1 :
            name1,layer1=full_name.split(",")
            name=name1.replace('"','')
            layer=layer1.replace('"','')
            if name not in dico_structure :
                dico_structure[name] = dict()
                dico_structure[name][layer] = 1
            else :
                if layer not in dico_structure[name] :
                    dico_structure[name][layer]=1
                else :
                    dico_structure[name][layer]=dico_structure[name][layer]+1
        else :
            if j==2 :
                name1,subname,layer1=full_name.split(",")
                name = name1.replace('"', '')
                layer = layer1.replace('"', '')
                if name not in dico_structure:
                    dico_structure[name] = dict()
                    dico_structure[name][subname] =dict()
                    dico_structure[name][subname][layer]=1
                else :
                    if subname not in dico_structure[name] :
                        dico_structure[name][subname] = dict()
                        dico_structure[name][subname][layer] = 1
                    else :
                        if layer not in dico_structure[name][subname] :
                            dico_structure[name][subname][layer] = 1
                        else:
                            dico_structure[name][subname][layer]=dico_structure[name][subname][layer]+1
        if j==0 :
            if full_name not in dico_structure :
                dico_structure[full_name]=1
            else :
                dico_structure[full_name]=dico_structure[full_name]+1
        i=i+1
    return(dico_structure)


def plot_number_structures (dico_structure) : #returns plots of the number of cells within the structure/substructure/layer for a given specimen ; need of the dictionary of the structure of the specimen
    structures=list(dico_structure.keys()) #liste des noms de structures : primary visual area, etc. =21 structures
    test=type(dico_structure[structures[0]]) is dict
    if test==False :
        li=list(dico_structure.keys())
        vall=list(dico_structure.values())
        plt.bar(li,vall,0.35)
        plt.xticks(rotation=15)
        plt.title("human")
        plt.show()
    else :
        for i in structures :
            ss_structure=list(dico_structure[i].keys()) #pour un nom de structures, on prend ses keys : peut être les layers ou des sous-structures
            t=type(dico_structure[i][ss_structure[0]]) is dict #on regarde si la première key de la structure est un dico ou non : si est un dico c'est une sous structure, si non c'est une layer
            if t==True :
                j=0
                r=0
                while j<len(ss_structure) : #pour chaque sous structure
                    if type(dico_structure[i][ss_structure[j]])is dict:
                        labels=list(dico_structure[i][ss_structure[j]].keys()) #liste des layers présentes dans la sous structure
                        val=list(dico_structure[i][ss_structure[j]].values()) #liste du nbre de cell pour chaque layer
                        width=0.35+r
                        plt.bar(labels, val, width, label=str(ss_structure[j]))
                        plt.title(str(i))
                        plt.legend()
                    else :
                        label = ss_structure[j]
                        val=dico_structure[i][ss_structure[j]]
                        width = 0.35
                        plt.bar(label, val, width, label=str(ss_structure[j]))
                        plt.title(str(i))
                        plt.legend()
                    j=j+1
                    r=r-0.1
                plt.show()

            else :
                labels = ss_structure
                val=list(dico_structure[i].values())
                width = 0.35
                fig, ax = plt.subplots()
                ax.bar(labels, val, width, label='layers')
                ax.set_title(str(i))
            plt.show()

#electrophysiology

def cell_ephys_info (id_specimen) : #return the cell electrophysiological info for a given id
    all_cell_features = ctc.get_ephys_features()
    for i in all_cell_features :
        if i["specimen_id"]==id_specimen :
            return i

def number_of_sweeps (id_specimen) :
    all_sweeps = ctc.get_ephys_sweeps(id_specimen)
    print(len(all_sweeps))

def sweep_info (id_specimen,sweep_nber): #return data info of a given sweep of a given id specimen
    data_set = ctc.get_ephys_sweeps(id_specimen)
    for i in data_set :
        if i['sweep_number']==sweep_nber:
            return i



def plot_stimulus (id_specimen,sweep) : #plot the stimulus and the cell response (in mV) ; to choose the sweep number, you can use the stimulus_type function to know which sweep (by their numbers) does which type of stimulus
    data_set = ctc.get_ephys_data(id_specimen)
    all_sweeps = ctc.get_ephys_sweeps(id_specimen)
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

def stimulus_type (id_specimen): #return a dictionry with the keys being the different stimulus types and the values are the sweep numbers for a given id
    all_sweeps = ctc.get_ephys_sweeps(id_specimen)
    dico_stimuli=dict()
    for i in all_sweeps :
        number = i["sweep_number"]
        if i["stimulus_name"] not in dico_stimuli :
            dico_stimuli[i["stimulus_name"]]=[number]
        else :
            dico_stimuli[i["stimulus_name"]].append(number)
    return (dico_stimuli)




def ephys_spike_key_extractor(id_specimen,sweep,t_start,t_end): #return the keys/feature names caracterizing the spike ; if don't want to put t_start/end, make them equal to None
    dict_stimulus = stimulus_type(id_specimen)
    if sweep in dict_stimulus['Test']:
        return ('no data to return')
    else:
        data_set = ctc.get_ephys_data(id_specimen)
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





def spike_info (id_specimen,sweep_number,key) : #with a given feature name, and for a given sweep of a given id specimen, this function returns the data of the feature
    data_set = ctc.get_ephys_data(id_specimen)
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

def spiny_state (dict_specimen) :
    spiny=[]
    aspiny=[]
    sparsely_spiny=[]
    for i in dict_specimen :
        if i['tag__dendrite_type']=='spiny' :
            spiny.append(i['specimen__id'])
        else :
            if i['tag__dendrite_type'] == 'aspiny':
                aspiny.append(i['specimen__id'])
            else :
                sparsely_spiny.append(i['specimen__id'])
    return(spiny,aspiny,sparsely_spiny)

def reconstruction(dict_specimen) :
    no=[]
    yes=[]
    for i in dict_specimen :
        if i["nr__reconstruction_type"]==None :
            no.append(i["specimen__id"])
        else :
            yes.append(i["specimen__id"])
    return (no,yes)

def spiny_reconstruction (dict_specimen) :
    spiny_no_reconstruction=[]
    spiny_reconstruction=[]
    aspiny_no_reconstruction=[]
    aspiny_reconstruction=[]
    sparsely_no_reconstruction=[]
    sparsely_reconstruction=[]
    no,yes=reconstruction(dict_specimen)
    spiny,aspiny,sparsely=spiny_state(dict_specimen)
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

def morpho_info(id_specimen): #for a given id, returns the morpho info of the cell
    all_morpho = ctc.get_morphology_features()
    for i in all_morpho:
        if i["specimen_id"] == id_specimen:
            return i

def swc_morpho(id_specimen,name) : #return a link to download the morphology swc file of a given cell
    swc_file=ctc.get_reconstruction(id_specimen,name)
    return(swc_file)

def compartment_nber (id_specimen) :
    all_morpho=ctc.get_reconstruction(id_specimen)
    compartment_number=len(all_morpho.compartment_list)
    return(compartment_number)


def twoD_morpho(id_specimen):
    morphology = ctc.get_reconstruction(id_specimen)
    markers = ctc.get_reconstruction_markers(id_specimen)
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



ctc = CellTypesCache(manifest_file='cell_types/manifest.json') #create the manifest file

cells = ctc.get_cells(file_name="all_cells",require_reconstruction=True)
id_mouse,id_human=take_id("all_cells")
#print(id_mouse)
dict_mouse,dict_human=dict_specimen("all_cells")

#web_page(565871768)

#dico=all_cell_markers(dict_mouse)

#print(dico)

#print(cell_marker(565871768,dict_mouse))

struct_mouse=structure(dict_human)

#print(sweep_info(565871768,1))

cell=cell_info(52501190,dict_human)
#print(cell)

#plot_marker(dict_mouse)
#print(stimulus_type(565871768))

#print(stimulus_type(464212183))
#print(ephys_spike_key_extractor(565871768,94,0.4,0.6))

#print(cell_info(565871768,dict_mouse))
#twoD_morpho(565871768)


#print(list(marker_number(dict_mouse).keys()))

#print(reconstruction(dict_mouse)[1])
#twoD_morpho(475057898)

#print(all_cell_markers(dict_mouse).keys())

#a,b,c,d,e,f=spiny_reconstruction(dict_human)
#print(len(a),len(b),len(c),len(d),len(e),len(f))





#plot_marker(dict_mouse)

#'tag__dendrite_type'

#a=structure(dict_human)
#print(a)
#plot_number_structures(a)

#web_page(529878215)
#print(a)
plot_stimulus(464212183,102)

#print(morphology.compartment_list[0])
#morphology = ctc.get_reconstruction(480114344)
#web_page(464212183)
#plot_stimulus(464212183,35)

#plot_marker_number(dict_mouse)

#data_set = ctc.get_ephys_data(464212183)
#all_sweeps = ctc.get_ephys_sweeps(464212183)


#keys=ephys_spike_key_extractor(464212183,35)
#print(keys)


#a=ephys_spike_key_extractor(565871768,94)
#print(a)

#s_info=spike_info(464212183,35,'threshold_t')
#print(s_info)

#a=cell_info(464212183,dict_mouse)
#e=cell_type(474626527,dict_mouse)
#print(e)
#print(cell_info(464212183,dict_mouse))

#print(structure(dict_mouse))
#print(dict_mouse[0])








d=marker_number(dict_human)
#print(d)




#plot_marker_number(dict_mouse)

#print(dict_human)
